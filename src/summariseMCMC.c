/**********************************************************************
 * This is a function to summarise the MCMC output for a single gene  *
 * Just the summary statistics will be provided as the output         *
 **********************************************************************/

/* Commenced by Steve Pederson on 8-Oct-2010 */

/* Last updated on 8-Oct-2010 */

/* Written to be accessed using .Call */

/* The object that must be supplied is the MCMC output */

/* The output should be a matrix with all the required summary information (Mean, SD, Quantiles, Rhat, nEff etc) for each parameter */

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

/* The necessary headers */
#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include<Rmath.h>

double quickCov(double *x, double *y, int nChains) {
    // This will just be a quick calculation of covariance without any checks.
    // These can be performed elsewhere.
    double cov=0;
    int i;
    double xmean=0, ymean=0;
    
    for (i=0; i<nChains; i++) {
	xmean += x[i]/nChains;
	ymean += y[i]/nChains;
    }
    
    for (i=0; i<nChains; i++) {
	cov += (x[i] - xmean)*(y[i] - ymean) / (nChains - 1);
    }
    
    return(cov);
}

SEXP summariseMCMC(SEXP data, SEXP transform) {

    /* data must be the output from BMEA.MCMC for a single gene, as a 3D array with nonNull dimnames(data)[[1]]*/
    /* transform is essantially a boolean variable where FALSE(0) indicates no logit transformation */

    // The required variables for the process will be
    SEXP dimNames, paramNames, outNames, colNames;//, class;
    SEXP dimData;
    int nParam, nChains, nSims, nEff=0; // nChains = nChains, nSims = nSims, nParam = nParam
    int curParam, i, j, cur_ij; 
    SEXP tempParam;
    double muhat=0, tempVar=0;
    double *p_data, *p_tempParam;
    int *p_transform;
    SEXP xdot, s2, xdot2;
    double *p_xdot, *p_s2, *p_xdot2;
    double B=0, W=0, rHat=0, sig2hat=0;
    double varB=0, varW=0, covWB=0;
    double postVar=0, varPostVar=0, postDf=0;
    SEXP curName; // Holds the current parameter name
    SEXP curType;
    char c='0', *p_c;

    // Required variables for the output
    SEXP output;

    // Protect the main object
    PROTECT(data=AS_NUMERIC(data));

    // Get the parameter names:
    PROTECT(dimNames=GetArrayDimnames(data));
    PROTECT(paramNames=VECTOR_ELT(dimNames,0));

    // The number of parameters
    PROTECT(dimData = getAttrib(data,R_DimSymbol));
    nParam = INTEGER_POINTER(dimData)[0];
    nSims = INTEGER_POINTER(dimData)[1];
    nChains = INTEGER_POINTER(dimData)[2];

    // Define the output matrix
    if (nChains>1){ // If there is more than one chain
	PROTECT(output=allocMatrix(REALSXP,nParam,9));
	PROTECT(colNames=allocVector(STRSXP,9)); // Holds the colnames
	SET_STRING_ELT(colNames,7,mkChar("rHat"));
	SET_STRING_ELT(colNames,8,mkChar("nEff"));
    }
    else {
	PROTECT(output=allocMatrix(REALSXP,nParam,7));
	PROTECT(colNames=allocVector(STRSXP,7)); // Holds the colnames
    }
    // Set the rest of the column names
    SET_STRING_ELT(colNames,0,mkChar("Mean"));
    SET_STRING_ELT(colNames,1,mkChar("Sd"));
    SET_STRING_ELT(colNames,2,mkChar("2.5%"));
    SET_STRING_ELT(colNames,3,mkChar("25%"));
    SET_STRING_ELT(colNames,4,mkChar("50%"));
    SET_STRING_ELT(colNames,5,mkChar("75%"));
    SET_STRING_ELT(colNames,6,mkChar("97.5%"));

    PROTECT(outNames=allocVector(VECSXP,2)); // Holds the dimnames
    SET_VECTOR_ELT(outNames,0,paramNames);
    SET_VECTOR_ELT(outNames,1,colNames);
    setAttrib(output, R_DimNamesSymbol, outNames);

    // Set up the pointers to the data
    PROTECT(tempParam = allocVector(REALSXP,nSims*nChains));
    p_data = NUMERIC_POINTER(data);
    p_tempParam = NUMERIC_POINTER(tempParam);
    p_c = &c; // Set the pointer to the address of c
    p_transform = INTEGER_POINTER(transform); // Set the boolean pointer 

    // The values requied for estimation of rHat are PsiBar_j (mean for chain j) & SumSq_j (var for chain j * (nSims-1))
    PROTECT(xdot = allocVector(REALSXP,nChains));
    PROTECT(s2 = allocVector(REALSXP,nChains));
    PROTECT(xdot2 = allocVector(REALSXP,nChains));
    p_xdot = NUMERIC_POINTER(xdot);
    p_s2 = NUMERIC_POINTER(s2);
    p_xdot2 = NUMERIC_POINTER(xdot2);

    // Set up the curType variable as a STRSXP of length 3
    PROTECT(curType=allocVector(STRSXP,3));
    SET_STRING_ELT(curType, 0, mkChar("S"));     // For the 'S' parameters
    SET_STRING_ELT(curType, 1, mkChar("sig"));   // For the sigma parameters
    SET_STRING_ELT(curType, 2, mkChar("phi"));   // For the phi parameters
    PROTECT(curName=allocVector(STRSXP,2)); // The values are set in the loop for each parameter

    // Righto, now find the stats for each parameter;
    for (curParam=0; curParam<nParam; curParam++) {

	// Reset the appropriate variables
	muhat = 0;
	tempVar = 0;
	SET_STRING_ELT(curName,0,mkCharLen(CHAR(STRING_ELT(paramNames,curParam)),1)); // This will be the first character in the parameter name
	SET_STRING_ELT(curName,1,mkCharLen(CHAR(STRING_ELT(paramNames,curParam)),3)); // This will be the first 3 characters in the parameter name

	// Decide whether the variable should be log or logit transformed:
	// If it is one of the 'S' or 'sigma' variables, it should be log-transformed (c==1)
	// The phi variables should be logit transformed. (c==2)
	// Otherwise (c==0) no data transformation will be done.
	if (strcmp(CHAR(STRING_ELT(curName,0)), CHAR(STRING_ELT(curType,0))) == 0 && p_transform[0]!=0) *p_c='1'; // Checks for 'S' if it is to be transformed
	else if (strcmp(CHAR(STRING_ELT(curName,1)), CHAR(STRING_ELT(curType,1))) == 0 ) *p_c='1'; // Checks for 'sig'
	else if (strcmp(CHAR(STRING_ELT(curName,1)), CHAR(STRING_ELT(curType,2))) == 0 && p_transform[0]!=0) *p_c='2'; // Checks if phi is to be transformed to the logit scale
	else *p_c='0';

	// The means:
	for (j=0; j<nChains; j++) {
	    // Reset the chain specific variables
	    p_xdot[j] = 0;
	    p_s2[j] = 0;

	    for (i=0; i<nSims; i++) {
		cur_ij = j*nSims + i;
		p_tempParam[cur_ij] = p_data[curParam + cur_ij*nParam]; // Set the vector of the current parameter for later sorting
		p_xdot[j] += p_data[curParam + cur_ij*nParam] / nSims; // The chain specific mean
	    }
	    muhat += p_xdot[j] / nChains; // Calculate the overall mean, as each chain will have an equal number of sims
	}

	R_rsort(p_tempParam, nSims*nChains); // Sort in ascending order

	// The standard deviations:
	for(j=0; j<nChains; j++) {
	    for(i=0; i<nSims; i++) {
		cur_ij = j*nSims + i;
		tempVar += pow(p_tempParam[cur_ij]-muhat,2)/ (nSims*nChains - 1);
	    }
	}

	// Now assign these values to the output
	NUMERIC_POINTER(output)[curParam] = muhat;
	NUMERIC_POINTER(output)[nParam + curParam] = sqrt(tempVar);
	NUMERIC_POINTER(output)[2*nParam + curParam] = p_tempParam[ imax2(0,(int)(0.025*nSims*nChains-1)) ] ; // The 2.5th Percentile
	NUMERIC_POINTER(output)[3*nParam + curParam] = p_tempParam[ imax2(0,(int)(0.25*nSims*nChains-1)) ] ; // The 25th Percentile
	NUMERIC_POINTER(output)[4*nParam + curParam] = p_tempParam[ imax2(0,(int)(0.5*nSims*nChains-1)) ] ; // The 50th Percentile
	NUMERIC_POINTER(output)[5*nParam + curParam] = p_tempParam[ imax2(0,(int)(0.75*nSims*nChains-1)) ] ; // The 75th Percentile
	NUMERIC_POINTER(output)[6*nParam + curParam] = p_tempParam[ imax2(0,(int)(0.975*nSims*nChains-1)) ] ; // The 97.5th Percentile

	// The rest of the calculations for rHat should only be done if nChains>1:
	if (nChains>1) {

	    B = 0;
	    W = 0;
	    varB = 0;
	    varW = 0;
	    covWB = 0;

	    // Calculate the convergence statistics here. If the data doesn't need to be transformed (c==0), it's easy.
	    // Otherwise (c==1 || c==2), some variables will need to be recalculated (psiBar)
	    // These also only can be calculated if tempVar!=0
	    if (tempVar<0.00000001) { // This will skip the following calculations if it is a degenerate chain
		rHat = 1; // Maybe this should be NA
		nEff = 1;
	    }
	    else {
		
		switch(c) {
		    
		case '0': // If no data transformation is required, xdot & muhat have already been done during calculation of the mean

		    // The sum of squares:
		    for(j=0; j<nChains; j++) {
		    
			p_s2[j] = 0; // Reset the variable
		    
			for(i=0; i<nSims; i++) {
			    cur_ij = j*nSims + i;
			    p_s2[j] += pow(p_data[curParam + cur_ij*nParam] - p_xdot[j], 2) / (nSims -1); // The within chain variance
			}
		    }
		
		    for (j=0; j<nChains; j++) {
			B += nSims*pow(p_xdot[j]-muhat,2) / (nChains - 1); // The between chains variance
			W += p_s2[j] / nChains;
		    }
		
		    for (j=0; j<nChains; j++) {
			varW += pow(p_s2[j] - W, 2) / (pow(nChains, 2) - nChains);
			p_xdot2[j] = pow(p_xdot[j], 2); // Required for the covariance calculations
		    }
		    covWB = (nSims/nChains)*(quickCov(p_s2,p_xdot2,nChains) - 2*muhat*quickCov(p_s2, p_xdot, nChains)) ;
		    varB = pow(B, 2) * 2 / (nChains - 1);
		    sig2hat = (nSims-1)*W/nSims + B/nSims;
		    postVar = sig2hat + B/(nChains*nSims);
		    varPostVar = (pow((double)nSims - 1, 2)*varW + pow(1 + 1/(double)nChains, 2)*varB + 2*((double)nSims-1)*(1 + 1/(double)nChains)*covWB)/pow(nSims, 2);
		    if (varPostVar<0) varPostVar = 0;
		    postDf = 2*pow(postVar, 2)/varPostVar;
		    if (postDf >1000) postDf = 1000;
		    rHat = sqrt((postVar/W)*(postDf+3)/(postDf+1));
		    nEff = nChains*nSims*sig2hat/B;
		    if (nEff > nChains*nSims) nEff = nChains*nSims;
		
		    break;
		
		case '1': // For log-transformed data:
		
		    muhat = 0; // Reset for use here
		
		    // Calculate xdot
		    for (j=0; j<nChains; j++) {
			
			// Reset the chain specific variables
			p_xdot[j] = 0;
			p_s2[j] = 0;
			
			for (i=0; i<nSims; i++) {
			    cur_ij = j*nSims + i;
			    p_xdot[j] += log(p_data[curParam + cur_ij*nParam]) / nSims; // The chain specific mean (log_transformed)
			}
		    }
		    
		    // The sum of squares:
		    for(j=0; j<nChains; j++) {
			for(i=0; i<nSims; i++) {
			    cur_ij = j*nSims + i;
			    p_s2[j] += pow(log(p_data[curParam + cur_ij*nParam]) - p_xdot[j], 2) / (nSims -1); // The within chain variance
			}
		    }
		    
		    // Now the convergence statistics:
		    for (j=0; j<nChains; j++) {
			muhat += p_xdot[j] / nChains;
		    }
		    for (j=0; j<nChains; j++) {
			B += nSims*pow(p_xdot[j]-muhat,2) / (nChains - 1); // The between chains variance
			W += p_s2[j] / nChains;
		    }
		    
		    for (j=0; j<nChains; j++) {
			varW += pow(p_s2[j] - W, 2) / (pow(nChains, 2) - nChains);
			p_xdot2[j] = pow(p_xdot[j], 2); // Required for the covariance calculations
		    }
		    covWB = (nSims/nChains)*(quickCov(p_s2,p_xdot2,nChains) - 2*muhat*quickCov(p_s2, p_xdot, nChains)) ;
		    varB = pow(B, 2) * 2 / (nChains - 1);
		    sig2hat = (nSims-1)*W/nSims + B/nSims;
		    postVar = sig2hat + B/(nChains*nSims);
		    varPostVar = (pow((double)nSims - 1, 2)*varW + pow(1 + 1/(double)nChains, 2)*varB + 2*((double)nSims-1)*(1 + 1/(double)nChains)*covWB)/pow(nSims, 2);
		    if (varPostVar<0) varPostVar = 0;
		    postDf = 2*pow(postVar, 2)/varPostVar;
		    if (postDf >1000) postDf = 1000;
		    rHat = sqrt((postVar/W)*(postDf+3)/(postDf+1));
		    nEff = nChains*nSims*sig2hat/B;
		    if (nEff > nChains*nSims) nEff = nChains*nSims;
		    
		    break;
		    
		case '2': // For logit-transformed data:
		    
		    muhat = 0; // Reset for use here
		    
		    // Calculate xdot
		    for (j=0; j<nChains; j++) {
			
			// Reset the chain specific variables
			p_xdot[j] = 0;
			p_s2[j] = 0;
			
			for (i=0; i<nSims; i++) {
			    cur_ij = j*nSims + i;
			    p_xdot[j] += log(p_data[curParam + cur_ij*nParam]/(1-p_data[curParam + cur_ij*nParam])) / nSims; // The chain specific mean (log_transformed)
			}
		    }
		    
		    // The sum of squares:
		    for(j=0; j<nChains; j++) {
			for(i=0; i<nSims; i++) {
			    cur_ij = j*nSims + i;
			    p_s2[j] += pow(log(p_data[curParam + cur_ij*nParam]/(1-p_data[curParam + cur_ij*nParam])) - p_xdot[j], 2) / (nSims -1); // The within chain variance
			}
		    }
		    
		    // Now the convergence statistics:
		    for (j=0; j<nChains; j++) {
			muhat += p_xdot[j] / nChains;
		    }
		    for (j=0; j<nChains; j++) {
			B += nSims*pow(p_xdot[j]-muhat,2) / (nChains - 1); // The between chains variance
			W += p_s2[j] / nChains;
		    }
		    
		    for (j=0; j<nChains; j++) {
			varW += pow(p_s2[j] - W, 2) / (pow(nChains, 2) - nChains);
			p_xdot2[j] = pow(p_xdot[j], 2); // Required for the covariance calculations
		    }
		    covWB = (nSims/nChains)*(quickCov(p_s2,p_xdot2,nChains) - 2*muhat*quickCov(p_s2, p_xdot, nChains)) ;
		    varB = pow(B, 2) * 2 / (nChains - 1);
		    sig2hat = (nSims-1)*W/nSims + B/nSims;
		    postVar = sig2hat + B/(nChains*nSims);
		    varPostVar = (pow((double)nSims - 1, 2)*varW + pow(1 + 1/(double)nChains, 2)*varB + 2*((double)nSims-1)*(1 + 1/(double)nChains)*covWB)/pow(nSims, 2);
		    if (varPostVar<0) varPostVar = 0;
		    postDf = 2*pow(postVar, 2)/varPostVar;
		    if (postDf >1000) postDf = 1000;
		    rHat = sqrt((postVar/W)*(postDf+3)/(postDf+1));
		    nEff = nChains*nSims*sig2hat/B;
		    if (nEff > nChains*nSims) nEff = nChains*nSims;
		    break;
		    
		}    
	    }

	    // Place the values in the output
	    NUMERIC_POINTER(output)[7*nParam + curParam] = rHat;
	    NUMERIC_POINTER(output)[8*nParam + curParam] = nEff;
	}

    }

    /* Protected objects:
     * 1 - data
     * 2 - dimNames
     * 3 - paramNames
     * 4 - dimData
     * 5 - output
     * 6 - outNames
     * 7 - colNames
     * 8 - tempParam
     * 9 - xdot
     *10 - s2
     *11 - curName
     *12 - curType
     *13 - xdot2
     */

    UNPROTECT(13);

    return(output);

    /***********************************
     * End Of Function summariseChains *
     ***********************************/

}
