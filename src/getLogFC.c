/********************************************************************
 * Function for estimating logFC from the sims for each value of mu *
 ********************************************************************/

/* Begun by Steve Pederson Oct 2010 */

/* The supplied input should be a SEXP mu which is an H*(nSims*nChains) matrix and 
 * a SEXP contrasts which is a H*nContrasts matrix. 
 * The chains should have already been merged into one matrix in R */

#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include<Rmath.h>


SEXP getLogFC(SEXP mu, SEXP contrasts, SEXP includeSims) {

    // mu is the matrix of sims for mu. nRows = sims, nCols = H
    // contrasts is the matrix of contrasts

    SEXP dimMu, dimContrasts;
    SEXP tempSims, sims, summary; // Hold the sims & the summary respectively
    SEXP output; // The output
    int H, N, M; // Hold the number of conditions, total sims & total contrasts
    int h, n, m; // Used for the looping
    int *p_incSims; // Determines whether the sims are to be included in the output or not
    SEXP contNames, summaryColNames, summaryNames, listNames, simNames;
    double *p_sims, *p_summary, *p_contrasts, *p_mu; 
    double curMean=0, curVar=0;
    double nGr0=0, nLs0=0, B; // Hold the numbers greater or less than zero
    
    // Protect the input matrices
    PROTECT(mu=AS_NUMERIC(mu));
    PROTECT(contrasts=AS_NUMERIC(contrasts));
    PROTECT(includeSims=AS_INTEGER(includeSims));
    p_contrasts = NUMERIC_POINTER(contrasts);
    p_mu = NUMERIC_POINTER(mu);
    p_incSims = INTEGER_POINTER(includeSims);
    
    // The number of parameters & contrasts
    PROTECT(dimMu = getAttrib(mu,R_DimSymbol));
    H = INTEGER_POINTER(dimMu)[1]; // The number of conditions/cellTypes
    N = INTEGER_POINTER(dimMu)[0]; // The total number of sims
    PROTECT(dimContrasts = getAttrib(contrasts, R_DimSymbol));
    M = INTEGER_POINTER(dimContrasts)[1]; // The number of contrasts


    // Setup the sims & summary
    PROTECT(sims = allocMatrix(REALSXP, N, M)); // Each col is a contrast, each row is a sim
    PROTECT(summary = allocMatrix(REALSXP, M, 9)); // Each row is a contrast, the columns are mean, sd & the quantiles
    PROTECT(tempSims = allocVector(REALSXP, N));
    p_sims = NUMERIC_POINTER(tempSims);
    p_summary = NUMERIC_POINTER(summary);

    // The actual process
    for (m=0; m<M; m++) { // Do it for each contrast

	// Reset the variables
	curMean = 0;
	curVar = 0;

	for (n=0; n<N; n++) { // For each sim

	    p_sims[n] = 0; // Reset the value

	    for (h=0; h<H; h++) { // For each condition
		p_sims[n] += p_mu[h*N + n]*p_contrasts[H*m + h]; // Sample the value of logFC for each sim
	    }

	    curMean += p_sims[n]/N; // Find the cumulative mean
	}

	for (n=0; n<N; n++) {
	    curVar += pow(p_sims[n] - curMean, 2) / (N-1); // Find the cumulative sd
	}

	if (p_incSims[0] != 0 ) { // If the sims matrix is required for the output
	    for (n=0; n<N; n++) {
		NUMERIC_POINTER(sims)[m*N + n] = p_sims[n];
	    }
	}

	R_rsort(p_sims, N); // Sort in ascending order
	
	// Assign the mean & Sd to the output
	p_summary[m] = curMean;
	p_summary[M + m] = sqrt(curVar);
	p_summary[2*M + m] = p_sims[ imax2(0,(int)(0.025*N-1)) ] ; // The 2.5th Percentile
	p_summary[3*M + m] = p_sims[ imax2(0,(int)(0.25*N-1)) ] ; // The 25th Percentile
	p_summary[4*M + m] = p_sims[ imax2(0,(int)(0.5*N-1)) ] ; // The 50th Percentile
	p_summary[5*M + m] = p_sims[ imax2(0,(int)(0.75*N-1)) ] ; // The 75th Percentile
	p_summary[6*M + m] = p_sims[ imax2(0,(int)(0.975*N-1)) ] ; // The 97.5th Percentile
	
	// Calculate B
	nGr0 = 0;
	nLs0 = 0;
	if (p_sims[0] > 0) {
	    nGr0 = N;
	}
	else if (p_sims[N] < 0 ) {
	    nLs0 = N;
	}
	else {
	    for (n=0; n<N; n++) {
		if (p_sims[n] <=0) nLs0 ++;
		if (p_sims[n] >=0) nGr0 ++;
	    }
	}
	B = log(((double)nGr0+0.5)/((double)nLs0+0.5)); // Scaled as per the empirical logit transformation
	// Enter it in the output 
	p_summary[7*M + m] = fmax2((double)nGr0/(double)N, (double)nLs0/(double)N);
	p_summary[8*M + m] = B; // The B statistic

    }

    // Setup the names for the output
    PROTECT(contNames = getAttrib(contrasts, R_DimNamesSymbol));
    // The column names for the summary matrix
    PROTECT(summaryColNames = allocVector(STRSXP, 9));
    SET_STRING_ELT(summaryColNames,0,mkChar("Mean"));
    SET_STRING_ELT(summaryColNames,1,mkChar("Sd"));
    SET_STRING_ELT(summaryColNames,2,mkChar("2.5%"));
    SET_STRING_ELT(summaryColNames,3,mkChar("25%"));
    SET_STRING_ELT(summaryColNames,4,mkChar("50%"));
    SET_STRING_ELT(summaryColNames,5,mkChar("75%"));
    SET_STRING_ELT(summaryColNames,6,mkChar("97.5%"));
    SET_STRING_ELT(summaryColNames,7,mkChar("maxP"));
    SET_STRING_ELT(summaryColNames,8,mkChar("B"));
    // The dimnames for the summary matrix
    PROTECT(summaryNames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(summaryNames, 0, VECTOR_ELT(contNames,1));
    SET_VECTOR_ELT(summaryNames, 1, summaryColNames);
    setAttrib(summary, R_DimNamesSymbol, summaryNames);
    

    // Setup the output as a list with the summary (& sims if required)
    PROTECT(output=allocVector(VECSXP,2)); // A list with a slot for the summary ([[1]]) & the sims ([[2]])
    SET_VECTOR_ELT(output, 0, summary);
    if (p_incSims[0] != 0) { // Otherwise include the sims
	PROTECT(simNames =allocVector(VECSXP, 2));
	SET_VECTOR_ELT(simNames, 1, VECTOR_ELT(contNames, 1));
	SET_VECTOR_ELT(output, 1, sims); 
	setAttrib(sims, R_DimNamesSymbol, simNames);
    }
    // Set the names for the list
    PROTECT(listNames=allocVector(STRSXP,2));
    SET_STRING_ELT(listNames, 0, mkChar("summary"));
    SET_STRING_ELT(listNames, 1, mkChar("sims"));
    setAttrib(output, R_NamesSymbol, listNames);


    /* Protected items:
     * 1 - mu
     * 2 - contrasts
     * 4 - includeSims
     * 5 - dimMu
     * 6 - dimContrasts
     * 7 - sims
     * 8 - summary
     * 9 - output
     * 11 - contNames
     * 12 - summaryColNames
     * 13 - summaryNames
     * 14 - listNames
     * 15 - simNames (done separately)
     */
    if (p_incSims[0] != 0) UNPROTECT(1); // If sims are included, there will be an extra componenet to the outpu
    UNPROTECT(13);
//    UNPROTECT(6);

    return(output);

}
