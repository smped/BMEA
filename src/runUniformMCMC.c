/* The idea of this script is to write a wrapper that can receive data from R
 * then stick it into the updateValues function & run the entire MCMC process.
 * Begin by just running a single chain. */

/* Written in August 2010 */

/* As a process try structuring it like this:
 * Accept the data as an SEXP
 * Define the structure of the data
 * Run the checks
 * Generate initial values
 * Run the chain */

/* The data from R should include:
 * 1 - The observed PM data as an I*K vector (PM_hijk)
 * 2 - A vector defining which chip belongs to which condition (conditions)
 * 3 - A vector defining which probe belongs to which exon (exons)
 * 4 - A vector with the lambda_ik values
 * 5 - A vector with the delta_ik values
 * 6 - A list with the MCMC input parameters. This must be in the order: nChains, nIter, nBurnin, nThin
*/

/* All other data structures can be extracted from this information */

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

/* The necessary headers */
#include<R.h>
#include<Rmath.h>
#include<Rdefines.h>

#include "bmea.h" // Prototypes many sampling functions

/* The common structural variables for each gene will be: */
int H; // # of conditions
int I; // # of chips
int J; // # of exons
int K; // # of probes
int totParam; // This is the total number of parameters in the model

/* Square Root constants required of the structural variables */
//double SQRT_I;

/* The other variable vectors */
SEXP cond_count, exon_count; // Able to be defined as of variable size. Pointed to by *I_h & *K_j
int *I_h, *K_j; // Used to point to the SEXPs above as required. These are the counts of chips/probes for each h/j respectively
double *PM_hijk; // Points at the actual data
int *p_conditions, *p_exons; // Global pointers to the initial SEXP vectors

// The background parameters are easier if handled globally
double *lambda_ik, *delta_ik; // Define them as pointers to the vectors 'bg_means' & 'bg_sds' which will be sent from R

/* The MCMC parameters can also be defined globally */
int nChains, nThin, nKeep;
long nIter, nBurnin; // Just in case nIter > 32,768

/* The variances of the proposal generating distributions can also be declared globally for use in updateValues */
SEXP sigmaJump_mu, sigmaJump_p; // Can be of variable size to account for changes to nChains. Use the following pointers to access the values
double *tau_mu; // Set this initially to 0.5
double *tau_p; // Set this initially to 1

/* The common variables used for iterating will be: */
int h, i, j, k;
int cur_h, cur_i, cur_j, cur_ik, cur_hj;
int chain, iter, samp;

/* An object to store the inits */
SEXP inits;

/* Define the functions:
 * 1 - setDimensions (done)
 * 2 - initialiseValues (done)
 * 3 - updateValues (done)
 * 4 - runUniformMCMC */


/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/

static void setDimensions (SEXP conditions, SEXP exons){

    // Set I as the length of conditions & K as the length of exons:
    I = LENGTH(conditions);
    K = LENGTH(exons);

    // Set H as the maximum of the conditions vector:
    H=0; // Make sure that H starts at zero.
    for (i=0; i<I; i++) {
	if (p_conditions[i]>H) H=p_conditions[i];
    }

    // Set J as the maximum of the exons vector:
    J=0; // Make sure that J begins at zero for this process!
    for (k=0; k<K; k++) {
	if (p_exons[k]>J) J=p_exons[k];
    }

    // The total number of parameters
    totParam = I*K + I + H + K + H*J + 3;

    // Set the square roots. This will speed up the iterations considerably
    //    SQRT_I = sqrt(I);

    // Define cond_count as an integer SEXP of length H to help define the vector I_h
    PROTECT(cond_count=NEW_INTEGER(H));
    I_h = INTEGER_POINTER(cond_count); // Set I_h to point at the SEXP
    for (h=0; h<H; h++) I_h[h] = 0; // Make sure each value is intially a zero
    for (i=0; i<I; i++) {
	cur_h = p_conditions[i]-1;
	I_h[cur_h]++;
    }

    // Use the same process for the exons & K_j
    PROTECT(exon_count=NEW_INTEGER(J));
    K_j = INTEGER_POINTER(exon_count); // Set K_j to point at the SEXP
    for (j=0; j<J; j++) K_j[j] = 0; // Make sure each value is intially a zero
    for (k=0; k<K; k++) {
	cur_j = p_exons[k]-1;
	K_j[cur_j]++;
    }

    UNPROTECT(2);

    /*******************
     * End of Function *
     ******************/
}

/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/

static void initialiseValues (double *S_hijk, double *sigma_S, double *c_hi, double *mu_h, double *sigma_mu, double *p_k, double *sigma_p, double *phi_hj) {

    /* This will provide initial values for all the above parameter vectors.
     * The vectors are declared in the main MCMC function & the correct size will be defined */

    double sumSq, mu_p;
    double *p_inits;
    int cur_chainKeep; // Used for indexing the chain & iteration below

    /* Do each chain separately */
    for (chain=0; chain<nChains; chain++) {

	// Reset the cumulative variables
	for (h=0; h<H; h++) mu_h[chain*H + h] = 0;
	for (k=0; k<K; k++) p_k[chain*K + k] = 0;
	sumSq = 0;
	mu_p = 0;

	// Get the random.seed from R separately for each chain
	GetRNGstate();

	/* First do S_hijk */
	for (i=0; i<I; i++) {
	    for (k=0; k<K; k++) {
		cur_ik = i*K + k;
		S_hijk[chain*I*K + cur_ik] = sample_S(PM_hijk[cur_ik], lambda_ik[cur_ik], delta_ik[cur_ik]);
	    }
	}

	/* c_hi & mu_h */
	for (i=0; i<I; i++) {
	    samp = (int)runif(0, K); // Choose a random integer between 0 & K for the current I
	    cur_i = chain*I;
	    c_hi[cur_i + i ] = log(S_hijk[cur_i*K + i*K + samp]); // Take the log of the signal estimate at this probe on chip 'i' as the chip effect
	    cur_h = p_conditions[i]-1;
	    mu_h[chain*H + cur_h] += c_hi[cur_i + i] / I_h[cur_h]; // Take the average of the sampled c_hi values as the mu_h estimate
	}

	/* sigma_mu */
	for (i=0; i<I; i++){
	    cur_h = p_conditions[i]-1;
	    sumSq += pow(c_hi[chain*I + i] - mu_h[chain*H + cur_h], 2) / (I - 1); // The variance of the random sampled c_hi around the mu_h values
	}
	if (sumSq < 25) { // make sure that sigma_mu < 5 as defined by the prior
	    sigma_mu[chain] = sqrt(sumSq);
	}
	else {
	    sigma_mu[chain] = runif(0, 5);
	}

	/* phi_hj */
	if (J == 1) {
	    for (h=0; h<H; h++ ) {
		phi_hj[chain*H + h] = 1; // If there is only one exon, phi_hj=1 for all conditions
	    }
	}
	else {
	    for (h=0; h<H; h++ ) {
		for (j=0; j<J; j++) {
		    phi_hj[chain*H*J + h*J + j] = runif(0, 1); // A random sample from the prior
		}
	    }
	}

	/* p_k */
	for (k=0; k<K; k++) {
	    for (i=0; i<I; i++) {
		cur_ik = i*K + k;
		cur_hj = (p_conditions[i] - 1)*J + p_exons[k] - 1;
		// Find each probe effect as the mean of log(S_hijk) - c_hi - log(phi_hj)
		p_k[chain*K + k] += (log(S_hijk[chain*I*K + cur_ik]) - c_hi[chain*I + i] - log(phi_hj[chain*H*J + cur_hj])) / I;
	    }
	    mu_p += p_k[chain*K + k]/K; // Calculate the mean as we go along...
	}

	/*  sigma_p */
	sumSq = 0;
	for (k=0; k<K; k++) {
	  sumSq += pow(p_k[chain*K + k] - mu_p, 2) / (K-1); // Calculate the initial variance
	}
	if (sumSq < 100){ // Ensure the intial value < 10 as defined by the prior
	  sigma_p[chain] = sqrt(sumSq);
	}
	else{
	  sigma_p[chain] = runif(0, 10);
	}

	/* sigma_S */
	// calculate a frequentist estimate
	sumSq = 0; // Reset the dummy variable
	for (i=0; i<I; i++){
	  for (k=0; k<K; k++){
	    cur_ik = i*K + k;
	    cur_hj = (p_conditions[i] - 1)*J + p_exons[k] - 1;
	    sumSq += pow(log(S_hijk[chain*I*K + cur_ik]) - c_hi[chain*I + i] - p_k[chain*K + k] - log(phi_hj[chain*H*J + cur_hj]) , 2);
	  }
	}
	sigma_S[chain] = sqrt(sumSq/(I*K - 1));

	PutRNGstate();
    }

    /* Set up & save the inits */
    PROTECT(inits = alloc3DArray(REALSXP, totParam, 1, nChains)); 
    p_inits = NUMERIC_POINTER(inits); // Assign the pointer    

    // Use the SEXP 'inits' for returning the data via the pointer p_inits
    for (chain=0; chain<nChains; chain++) {

	cur_chainKeep = totParam*chain;
	p_inits[cur_chainKeep + I*K] = sigma_S[chain]; //sigma_S
	p_inits[cur_chainKeep + I*K + 1 + I + H] = sigma_mu[chain]; // sigma_mu
	p_inits[cur_chainKeep + I*K + 2 + I + H + K] = sigma_p[chain]; // sigma_p

	for (i=0; i<I; i++) {
	    for(k=0; k<K; k++) {
		cur_ik = i*K + k;
		p_inits[cur_chainKeep + cur_ik] = S_hijk[chain*I*K + cur_ik]; // Store the S values first
	    }
	    p_inits[cur_chainKeep + I*K + 1 + i] = c_hi[chain*I + i]; // c_hi
	}

	for (h=0; h<H; h++) {
	    p_inits[cur_chainKeep + I*K + 1 + I + h] = mu_h[chain*H + h]; // mu_h
	    for (j=0; j<J; j++) {
		cur_hj = h*J + j;
		p_inits[cur_chainKeep + I*K + 3 + I + H + K + cur_hj] = phi_hj[chain*H*J + cur_hj]; // phi_hj
	    }
	}

	for (k=0; k<K; k++) {
	    p_inits[cur_chainKeep + I*K + 2 + I + H + k] = p_k[chain*K + k]; //p_k
	}

    }


    /*******************
     * End of Function *
     ******************/
}

/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/

static void updateValues (double *S_hijk, double *sigma_S, double *c_hi, double *mu_h, double *sigma_mu, double *p_k, double *sigma_p, double *phi_hj) {

    ////////////////////
    /* Updating order */
    ////////////////////

    /* 1 - S_hijk & sigma_S
     * 2 - sigma_mu, mu_h, c_hi
     * 3 - p_k, sigma_p
     * 4 - phi_hj */

    // These are used for each M-H step. These can be recycled for each new parameter
    double newValue, r, qOld, qNew, pOld, pNew;

    // Those specific to step 1
    int degf = I*K; // the degrees of freedom for the rchisq variable...
    double sumSq = 0; // holds the cumulative sum of squares

    // Those for step 2
    double xBar_hi[I];
    double sumSq_h[H]; // Holds the sum of squares
    double sigmaSq_h, sigmaSq_mu, sigmaSq_S;  
    double V_h[H], muHat_h[H]; // Used for updating mu_h
    double V_hi, thetaHat_hi[I]; // Used for updating c_hi

    // Those for step 3
    double xBar_k[K];
    double sigmaSq_k, V_k, theta_k, sigmaSq_p;

    // And those for step 4
    double xBar_hj[H*J];
    double sigma_hj[H*J];
    double theta_hj; // As this is only used temporarily in each loop, just declare it once.

/******************************
* Step 1: S_hijk then sigma_S *
*******************************/

    for (chain=0; chain<nChains; chain++) { // Do each chain separately

	GetRNGstate(); // Set the random.seed separately for each chain

	for (i=0; i<I; i++) {
	  
	  cur_i = chain*I;

	    for (k=0; k<K; k++){

		cur_ik = K*i + k; // The calculated index for the current i & k;
		cur_hj = J*(p_conditions[i] - 1) + p_exons[k] - 1; // The calculated index for the current h & j

		/* S */
		newValue = sample_S(PM_hijk[cur_ik], lambda_ik[cur_ik], delta_ik[cur_ik]); 

		// Calculate the values of the pdf for the old & new values
		pOld = exp(-0.5*( pow( (log( PM_hijk[cur_ik] - S_hijk[cur_i*K + cur_ik] ) - lambda_ik[cur_ik] ) / delta_ik[cur_ik], 2) + pow( (log( S_hijk[cur_i*K + cur_ik] / phi_hj[chain*H*J + cur_hj ] ) - c_hi[cur_i + i] - p_k[chain*K + k]) / sigma_S[chain], 2) ) ) / (S_hijk[cur_i*K + cur_ik]*(PM_hijk[cur_ik] - S_hijk[cur_i*K + cur_ik]));
		pNew =  exp(-0.5*( pow( (log( PM_hijk[cur_ik] - newValue ) - lambda_ik[cur_ik]) / delta_ik[cur_ik], 2) + pow( (log( newValue / phi_hj[chain*H*J + cur_hj ] ) - c_hi[cur_i + i] - p_k[chain*K + k]) / sigma_S[chain], 2) ) ) / (newValue*(PM_hijk[cur_ik] - newValue));

		// The probability of acceptance
		r = pNew / pOld ;
		if (r >1) r=1;

		// Sample with probability r & update each S_hijk value if successful
		samp = rbinom((double)1, r);
		if (samp == 1) S_hijk[cur_i*K + cur_ik] = newValue;

	    }
	}

	/* The sum of squares & degrees of freedom for sigma_S */
	sumSq = 0; // Reset the cumulative variables for each chain

	for (i=0; i<I; i++){
	    for (k=0; k<K; k++){

		cur_ik = K*i + k; // The calculated index for the current i & k;
		cur_hj = J*(p_conditions[i] - 1) + p_exons[k] - 1; // The calculated index for the current h & j
		sumSq += pow( log(S_hijk[chain*I*K + cur_ik]) - c_hi[chain*I + i] - p_k[chain*K + k] - log(phi_hj[chain*H*J + cur_hj]), 2); // Increase the cumulative sum of squares

	    }
	}

	/* The new value for sigmaS */
	sigma_S[chain]= sqrt(sumSq / rchisq((double)degf) );
	PutRNGstate();
    }

/********************************
* Step 2: sigma_mu, mu_h & c_hi *
*********************************/

    for (chain=0; chain<nChains; chain++){

	GetRNGstate();

	// Initialise the variables of length H which are incremented during the update at zero for each chain
	for (h=0; h<H; h++) {
	    sumSq_h[h] = 0;
	    muHat_h[h] = 0;
	}
	// Reset the pOld & pNew values for updating sigma_mu
	pOld = 1;
	pNew = 1;
	sigmaSq_S = pow(sigma_S[chain], 2);
	sigmaSq_mu = pow(sigma_mu[chain], 2);
	sigmaSq_h = sigmaSq_S / K; // This is used in all steps & is obtained after sigma_S has been updated

	for (i=0; i<I; i++){

	    xBar_hi[i] = 0; // Initialise this variable at zero. This is also used for mu_h so needs to be set for each 'i'
	    cur_h = p_conditions[i]-1; // Get which condition the chip belongs to

	    for (k=0; k<K; k++) {
		cur_ik = K*i + k; // The calculated index for the current i & k;
		cur_hj = J*(cur_h) + p_exons[k] - 1; // The calculated index for the current h & j
		xBar_hi[i] += ( log(S_hijk[chain*I*K + cur_ik]) - p_k[chain*K + k] - log(phi_hj[chain*H*J + cur_hj]) ) / K;
	    }

	    sumSq_h[ cur_h ] += pow( xBar_hi[i] - mu_h[chain*H + cur_h], 2) ; // Calculate the sum of squares
	}

	/************
	 * sigma_mu *
	 ************/
	newValue = r_truncnorm(0.0, 5.0, sigma_mu[chain], tau_mu[chain]);

	// The q values, i.e. P(new | old) & P(old | new)
	qNew = dnorm(newValue, sigma_mu[chain], tau_mu[chain], 0) / (pnorm(5, sigma_mu[chain], tau_mu[chain], 1, 0) - pnorm(0, sigma_mu[chain], tau_mu[chain], 1, 0));
	qOld = dnorm(sigma_mu[chain], newValue, tau_mu[chain], 0) / (pnorm(5, newValue, tau_mu[chain], 1, 0) - pnorm(0, newValue, tau_mu[chain], 1, 0));

	// Calculate the log p-values as cumulative products.
	for (h=0; h<H; h++) {
	    pOld *= exp(-0.5*((I_h[h]-1)*log( sigmaSq_h + sigmaSq_mu ) + sumSq_h[h] / ( sigmaSq_h + sigmaSq_mu ) ));
	    pNew *= exp(-0.5*((I_h[h]-1)*log( sigmaSq_h + pow(newValue, 2) ) + sumSq_h[h] / ( sigmaSq_h + pow(newValue, 2) ) ));
	}

	// Calculate r
	r = (pNew/pOld)*(qOld/qNew);
	if (r>1) r=1;

	// Sample the new value for sigma_mu with probability r
	samp = rbinom((double)1, r);
	if (samp == 1) sigma_mu[chain] = newValue;

	/********
	 * mu_h *
	 ********/
	for (i=0; i<I; i++) {
	    cur_h = p_conditions[i] - 1; // Get which condition the curent chip belongs to
	    muHat_h[cur_h] += xBar_hi[i] / I_h[cur_h]; // Cumulative posterior mean for the new value
	}

	for (h=0; h<H; h++) {
	    V_h[h] = (sigmaSq_h + sigmaSq_mu) / I_h[h]; // Posterior variance for each condtion
	    mu_h[chain*H + h] = r_truncnorm(0.0, scanLimit, muHat_h[h], sqrt(V_h[h]));
	}

	/********
	 * c_hi *
	 ********/
	V_hi = sigmaSq_mu*sigmaSq_S / ( K*sigmaSq_mu + sigmaSq_S);
	for (i=0; i<I; i++) {

	    cur_h = p_conditions[i] - 1;
	    thetaHat_hi[i] = V_hi*( K*xBar_hi[i] / sigmaSq_S + mu_h[chain*H + cur_h ] / sigmaSq_mu);
	    c_hi[chain*I + i] = r_truncnorm(0.0, scanLimit, thetaHat_hi[i], sqrt(V_hi)); 
	}

	PutRNGstate();
    }

/************************
* Step 3: sigma_p & p_k *
*************************/

    for (chain=0; chain<nChains; chain++) {

        
	GetRNGstate();

	sigmaSq_p = pow(sigma_p[chain], 2);	
	sigmaSq_k = pow(sigma_S[chain], 2) / I; // Calculate once sigma_S has been updated
	V_k = sigmaSq_k * sigmaSq_p / ( sigmaSq_p + sigmaSq_k);

	/*******
	 * p_k *
	 *******/
	for (k=0; k<K; k++) { // for each probe

	    xBar_k[k] = 0; //Set to 0
	    cur_j = p_exons[k] - 1;
	    cur_i = chain*I;

	    for (i=0; i<I; i++) {   // Cumulative mean
		cur_ik = i*K+k;
		cur_hj = (p_conditions[i] - 1)*J + cur_j;
		xBar_k[k] += (log( S_hijk[cur_i*K + cur_ik] ) - c_hi[cur_i + i] - log( phi_hj[chain*H*J + cur_hj ] )) / (double)I;
	    }
	    
	    theta_k = V_k*xBar_k[k] / sigmaSq_k; // The posterior mean
	    p_k[chain*K + k] = rnorm( theta_k, sqrt(V_k) ); // Assign the new value for that probe

	}

	/***********
	 * sigma_p *
	 ***********/
	newValue = r_truncnorm(0.0, 10.0, sigma_p[chain], tau_p[chain]);
	qNew = dnorm(newValue, sigma_p[chain], tau_p[chain], 0) / (pnorm(10, sigma_p[chain], tau_p[chain], 1, 0) - pnorm(0, sigma_p[chain], tau_p[chain], 1, 0)); // P(new | old)
	qOld = dnorm(sigma_p[chain], newValue, tau_p[chain], 0) / (pnorm(10, newValue, tau_p[chain], 1, 0) - pnorm(0, newValue, tau_p[chain], 1, 0)); // P(old | new)

	// Reset then calculate the log p-values as cumulative sums.
	pNew = 1;
	pOld = 1;
	for (k=0; k<K; k++) {
	    pOld *= exp(-0.5*(log(sigmaSq_k + sigmaSq_p) + pow(xBar_k[k], 2)/(sigmaSq_k + sigmaSq_p)));
	    pNew *= exp(-0.5*(log(sigmaSq_k + pow(newValue, 2)) + pow(xBar_k[k], 2)/(sigmaSq_k + pow(newValue, 2))));
	}

	// Calculate r
	r = (pNew/pOld)*(qOld/qNew);
	if (r>1) r=1;

	// Sample the new value with probability r
	samp = rbinom((double)1, r);
	if (samp == 1) sigma_p[chain] = newValue;
	PutRNGstate();

}

/*****************
* Step 4: phi_hj *
******************/
    if (J != 1) { // phi_hj will only be updated if there is more than one exon
	for (chain=0; chain<nChains; chain++) {

	    GetRNGstate();

	    // Set the key variables for this step
	    for (h=0; h<H; h++) {
		for (j=0; j<J; j++) {
		    sigma_hj[h*J +j] = sigma_S[chain] / sqrt(K_j[j]*I_h[h]); // NB: These could be converted to constants!!!
		    xBar_hj[h*J + j] = 0; // Initialise for incrementing in the next step
		}
	    }

	    for (i=0; i<I; i++) {

	        cur_i = chain*I;

		for (k=0; k<K; k++) {
		    cur_ik = i*K + k;
		    cur_h = p_conditions[i]-1;
		    cur_j = p_exons[k]-1;
		    cur_hj = cur_h*J + cur_j;
		    xBar_hj[cur_hj] += (log(S_hijk[cur_i*K + cur_ik]) - c_hi[cur_i + i] - p_k[chain*K + k]) / (double)(K_j[cur_j]*I_h[cur_h]); // Cumulative mean
		}
	    }

/********************************
 * Update the values for phi_hj *
 ********************************/
	    for (h=0; h<H; h++) {
		for (j=0; j<J; j++) {

		    cur_hj = h*J + j;
		    theta_hj = r_righttruncnorm(0.0, xBar_hj[cur_hj] + pow(sigma_hj[cur_hj], 2), sigma_hj[cur_hj]);
		    phi_hj[chain*H*J + cur_hj]  = exp( theta_hj ); // Update each value
		    		    
		}
	    }
	    PutRNGstate();
	}
    }

    /*******************
     * End Of Function *
     *******************/
}


/****************************************************************************************************/
/****************************************************************************************************/
/****************************************************************************************************/


/************************************************/
/* Define the output of the function as an SEXP */
/************************************************/

SEXP runUniformMCMC (SEXP data, SEXP conditions, SEXP exons, SEXP bg_means, SEXP bg_sds, SEXP mcmcParam) {

    /* Description of the input parameters: */
    /* 'data' holds the observed data as an I*K vector
     * 'conditions' defines which chip/array belongs to which condition/treatment
     * 'exons' defines which probes belongs to which exon
     * 'bg_means' & 'bg_sds' are the means & standard deviations for the BG prior for each probe
     * 'mcmcParam' is a list with components nChains, nIter, nBurnin, nThin
     */

    // Protect the main objects & coerce them to doubles & integers
    PROTECT(data=AS_NUMERIC(data));
    PROTECT(conditions=AS_INTEGER(conditions));
    PROTECT(exons=AS_INTEGER(exons));
    PROTECT(bg_means=AS_NUMERIC(bg_means));
    PROTECT(bg_sds=AS_NUMERIC(bg_sds));
    PROTECT(mcmcParam=AS_LIST(mcmcParam));

    // Make sure that PM_hijk points at the observed data & likewise for p_conditions & p_exons
    PM_hijk = NUMERIC_POINTER(data);
    p_conditions = INTEGER_POINTER(conditions);
    p_exons = INTEGER_POINTER(exons);
    // Set the background means
    lambda_ik = NUMERIC_POINTER(bg_means);
    delta_ik = NUMERIC_POINTER(bg_sds);

    // Set the MCMC parameters. These are all declared as global variables
    nChains = INTEGER(VECTOR_ELT(mcmcParam,0))[0]; // This MUST be the first element in the list
    nIter = INTEGER(VECTOR_ELT(mcmcParam,1))[0]; // Again, this MUST be the second element in the list
    nBurnin = INTEGER(VECTOR_ELT(mcmcParam,2))[0]; // The 3rd element
    nThin = INTEGER(VECTOR_ELT(mcmcParam,3))[0]; // The 4th element
    nKeep = (nIter-nBurnin)/nThin; // Set the number of iterations to keep

    // And setup tau_mu & tau_p with initial values of 0.2 & 0.25. These values are close to ideal & this should speed up the adapting process
    PROTECT(sigmaJump_mu=NEW_NUMERIC(nChains));
    PROTECT(sigmaJump_p=NEW_NUMERIC(nChains));
    tau_mu = NUMERIC_POINTER(sigmaJump_mu);
    tau_p = NUMERIC_POINTER(sigmaJump_p);
    for (chain=0; chain<nChains; chain++) {
	tau_mu[chain] = 1.0;
	tau_p[chain] = 0.25;
    }


/**************************************************************
 *Determine the dimensions of the current gene/dataset/matrix *
 **************************************************************/

    setDimensions(conditions, exons); // Call the above function to set H, I, J, K, I_h & K_J

    // Protect the global SEXPs now the dimensions have been defined
    PROTECT(cond_count=AS_INTEGER(cond_count));
    PROTECT(exon_count=AS_INTEGER(exon_count));

/**************************************************************************
 * Declare all the parameter vectors then provide initial values for them *
 **************************************************************************/

    double S_hijk[I*K*nChains], sigma_S[nChains]; // Holds the signal estimates & associated standard deviation
    double c_hi[I*nChains], mu_h[H*nChains], sigma_mu[nChains]; // The expression-level parameters
    double p_k[K*nChains], sigma_p[nChains]; // The probe-level parameters
    double phi_hj[H*J*nChains]; // The exon parameters

    // Initialise values for each chain
    initialiseValues(S_hijk, sigma_S, c_hi, mu_h, sigma_mu, p_k, sigma_p, phi_hj);


/*********************************
 * Perform the burnin iterations *
**********************************/

    // These can be discarded, but adapting needs to be done for tau_mu & tau_p
    // Later modifications could be to assess convergence while running the burnin phase & only exit
    // burnin once all rHat < 1.1
    double temp_mu[nBurnin*nChains], temp_p[nBurnin*nChains]; // Holds the values for sigma_mu & sigma_p
    int mu_done[nChains], p_done[nChains]; // A boolean variable indicating whether adapting is done
    float accRate_mu, accRate_p; // The acceptance rate for each. Use the last 100 samples to estimate this
    float mean_mu, mean_p; // Gets the mean for each normal approximation
    double var_mu=0, var_p; // Get the var for each normal approximation
    int t; // A local variable for iterating through during the adapting phase
    int nextAdapt = 100; // Used to set the next iteration to adapt tau_mu & tau_p

    // Initialise the boolean variables
    for (chain=0; chain<nChains; chain++) {
	mu_done[chain]=0;
	p_done[chain]=0;
    }

    if(nBurnin > nIter) nBurnin = nIter; // Make sure the burnin isn't longer the the number of iterations
    for (iter=0; iter<nBurnin; iter++) {

	// Update the values. All chains are done in this function
      updateValues(S_hijk, sigma_S, c_hi, mu_h, sigma_mu, p_k, sigma_p, phi_hj);

	for (chain=0; chain<nChains; chain++) {

	    // Add the lastest values
	    temp_mu[chain*nBurnin + iter] = sigma_mu[chain];
	    temp_p[chain*nBurnin + iter] = sigma_p[chain];

	    if (iter==nextAdapt) { // Begin adapting every 100 iterations after the first 100 iterations

		// The adapting phase for sigma_mu
		if (mu_done[chain]==0) {
		    accRate_mu=0; // Reset the acceptance rate to 100%
		    mean_mu = 0; // Reset the cumulative mean
		    var_mu = 0; // Reset the cumulative var
		    for (t=iter; t>iter-100; t--) { // Increment backwards through the last 100 samples
			if (temp_mu[chain*nBurnin + t]!=temp_mu[chain*nBurnin+ t-1]) accRate_mu++; // Reduce the acceptance rate if the value is the same as previous
			mean_mu += temp_mu[chain*nBurnin + t]/100; // Also get the mean as a cumulative sum
		    }

		    if (accRate_mu/100>0.42 && accRate_mu/100<0.46) {
			mu_done[chain]=1; // Declare adapting complete if the acceptance rate is ~0.44
		    }
		    else { // Otherwise keep adapting
			for (t=iter; t>iter-100; t--) {
			    var_mu += pow( temp_mu[chain*nBurnin + t] - mean_mu, 2) / 99; // Get the variance of the last 100
			}
			if (var_mu!=0) tau_mu[chain] = sqrt(2.4*var_mu); // Set the global variable to the sd of the last 100, if it is non-zero
		    }
		}

		// The adapting phase for sigma_p
		if (p_done[chain]==0) {
		    accRate_p=0; // Reset the acceptance rate to 100%
		    mean_p = 0; // Reset the cumulative mean
		    var_p = 0; // Reset the cumulative var
		    for (t=iter; t>iter-100; t--) { // Increment backwards through the last 100 samples
			if (temp_p[chain*nBurnin + t]!=temp_p[chain*nBurnin+ t-1]) accRate_p++; // Reduce the acceptance rate if the value is the same as previous
			mean_p += temp_p[chain*nBurnin + t]/100; // Also get the mean as a cumulative sum
		    }

		    if (accRate_p/100>0.42 && accRate_p/100<0.46) {
			p_done[chain]=1; // Declare adapting complete if the acceptance rate is ~0.44
		    }
		    else { // Otherwise keep adapting
			for (t=iter; t>iter-100; t--) {
			    var_p += pow( temp_p[chain*nBurnin + t] - mean_p, 2) / 99; // Get the variance of the last 100
			}
			if (var_p!=0) tau_p[chain] = sqrt(2.4*var_p); // Set the global variable to the sd of the last 100, if it is non-zero
		    }
		}
		////////////////////////////////////////////////////
		/*   Adapting has been checked for nBurnin = 500  */
		////////////////////////////////////////////////////
		/* No monitoring of convergence has been included */
		////////////////////////////////////////////////////
	    }

	}

	if (iter==nextAdapt) nextAdapt += 100; // Set the next adapting iteration

    }

/*********************************************************
 * The actual MCMC process where the iterations are kept *
 *********************************************************/

    // Set up the object for storing the output from the parameter updates:
    /* At this stage the output object should contain all parameters.
     * The total estimable parameters in the model is: totParam = I*K + I + H + K + H*J + 3
     * Signal Level: I*K + 1
     * Expression Level: I + H + 1
     * Probe Level: K + 1
     * Exon Level: H*J

     * Each chain must be recorded separately to calculate rHat
     * Each chain then needs an object of dimension totParam*nKeep

     * The best solution would be to create an SEXP array of dim(totParam, nKeep, nChains)

     * The iterations to be kept must also be indexed with a vector
     * Kept iterations would then be from the region between nBurnin+1:nIter */


    // Set up the vector indexing the iterations to be kept
    int keep[nKeep];
    int nextKeep = 0; // Used for accessing these values during iteration
    int cur_chainKeep; // Used for indexing the chain & iteration below

    t=0; // Reset t after being used in the adapting stage
    keep[t] = nBurnin + nThin - 1; // Set the first value to be kept
    do {
	t++;
	keep[t] = keep[t-1] + nThin; // Add nThin to the previous kept value
	if (keep[t] > nIter) { // If the next kept value is > nIter, abort here
	    keep[t] = nIter;
	    nKeep = t; // This will set a new value for nKeep and will restrict the size of any subsequent object to being correct
	}
    }
    while (t < nKeep || keep[t] < nIter); // Exit once the values exceed nIter

    // The best structure for the output is as a 3DArray in R with dimensions (totParam, nKeep, nChains)
    SEXP output, dimnames, paramNames, chainNames, class;
    double *p_output; // Points at the output SEXP
    char temp[12]; // Holds the temporary names for each parameter/chain

    // Protect & set up the output
    PROTECT(output = alloc3DArray(REALSXP, totParam, nKeep, nChains)); // The possible variability in nKeep means this needs to be placed here
    p_output = NUMERIC_POINTER(output); // Assign the pointer

    // Protect the names objects & define them
    PROTECT(dimnames = allocVector(VECSXP, 3));
    PROTECT(chainNames = allocVector(STRSXP,nChains));
    PROTECT(paramNames = allocVector(STRSXP,totParam));

    // Setup the chain names, i.e. chain1, chain2 etc...
    for (chain=0; chain<nChains; chain++) {
	sprintf(temp,"chain[%i]",chain+1);
	SET_STRING_ELT(chainNames, chain, mkChar(temp));
    }

    // Setup the parameter names
    SET_STRING_ELT(paramNames, I*K, mkChar("sigma_S"));
    SET_STRING_ELT(paramNames, I*K + 1 + I + H, mkChar("sigma_mu"));
    SET_STRING_ELT(paramNames, I*K + 2 + I + H + K, mkChar("sigma_p"));
    for (i=0; i<I; i++) {
	for (k=0; k<K; k++) {
	    sprintf(temp,"S[%i,%i]",i+1,k+1);
	    SET_STRING_ELT(paramNames, i*K + k, mkChar(temp)); // S_hijk
	    if (i==0) { // Just do the K values once
		sprintf(temp,"p[%i]",k+1);
		SET_STRING_ELT(paramNames, I*K + 2 + I + H + k, mkChar(temp)); // p_k
	    }
	}
	sprintf(temp,"c[%i]",i+1);
	SET_STRING_ELT(paramNames, I*K + 1 + i, mkChar(temp)); // c_hi
    }
    for (h=0; h<H; h++) {
	sprintf(temp, "mu[%i]",h+1);
	SET_STRING_ELT(paramNames, I*K + 1 + I + h, mkChar(temp)); // mu_h
	for (j=0; j<J; j++) {
	    sprintf(temp,"phi[%i,%i]",h+1,j+1);
	    SET_STRING_ELT(paramNames, I*K + 3 + I + H + K + h*J + j, mkChar(temp)); // phi_hj
	}
    }

    // Assign the names to the output array
    // The names for the 2nd dimension are left blank as these are just the samples. Indexing these is not important at all
    SET_VECTOR_ELT(dimnames, 0, paramNames); // paramNames must be defined above as an SEXP
    SET_VECTOR_ELT(dimnames, 2, chainNames); // chainNames must be defined above as an SEXP
    setAttrib(output, R_DimNamesSymbol, dimnames);
    // Also do this for the inits object
    setAttrib(inits, R_DimNamesSymbol, dimnames);

    // Define the class as a BMEA.MCMC object
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("BMEA.MCMC"));
    setAttrib(output, R_ClassSymbol, class);

    // Run the process
    for (iter=nBurnin; iter<nIter; iter++) {

        // Update the values. All chains are done automatically in this function using the global variables
      updateValues(S_hijk, sigma_S, c_hi, mu_h, sigma_mu, p_k, sigma_p, phi_hj);

	if (iter==keep[nextKeep]) {

	    // Use the SEXP 'output' for returning the data via the pointer p_output
	    for (chain=0; chain<nChains; chain++) {

		cur_chainKeep = totParam*nKeep*chain + totParam*nextKeep;
		p_output[cur_chainKeep + I*K] = sigma_S[chain]; //sigma_S
		p_output[cur_chainKeep + I*K + 1 + I + H] = sigma_mu[chain]; // sigma_mu
		p_output[cur_chainKeep + I*K + 2 + I + H + K] = sigma_p[chain]; // sigma_p

		for (i=0; i<I; i++) {
		    for(k=0; k<K; k++) {
			cur_ik = i*K + k;
			p_output[cur_chainKeep + cur_ik] = S_hijk[chain*I*K + cur_ik]; // Store the S values first
		    }
		    p_output[cur_chainKeep + I*K + 1 + i] = c_hi[chain*I + i]; // c_hi
		}

		for (h=0; h<H; h++) {
		    p_output[cur_chainKeep + I*K + 1 + I + h] = mu_h[chain*H + h]; // mu_h
		    for (j=0; j<J; j++) {
			cur_hj = h*J + j;
			p_output[cur_chainKeep + I*K + 3 + I + H + K + cur_hj] = phi_hj[chain*H*J + cur_hj]; // phi_hj
		    }
		}

		for (k=0; k<K; k++) {
		    p_output[cur_chainKeep + I*K + 2 + I + H + k] = p_k[chain*K + k]; //p_k
		}

	    }

	    nextKeep++; // Now increment to index the next kept iteration

	}
    }

    /* The protected objects so far are:
     * 1: conditions
     * 2: exons
     * 3: data
     * 4: bg_means
     * 5: bg_sds
     * 6: mcmcParam
     * 7: sigmaJump_mu
     * 8: sigmaJump_p
     * 9: cond_count
     * 10: exon_count
     * 11: output
     * 12: dimnames
     * 13: chainNames
     * 14: paramNames
     * 15: class
     */

    UNPROTECT(15);
    UNPROTECT(1); // The inits object
    
    return(output);
    // return(inits);
    // return(truncP_hijk); // Checked: The correct p_values are returned

    /*******************
     * End of Function *
     ******************/

}

/* Currently just returns NaN for everything... Bummer */
