/*********************************************************************
 * Function for estimating logFC from the sims for each value of phi *
 *********************************************************************/

/* Begun by Steve Pederson Oct 2010 */

/* The supplied input should be a SEXP phi which is an (H x J)*(nSims x nChains) matrix and 
 * a SEXP contrasts which is a H x nContrasts matrix. 
 * The chains should have already been merged into one matrix in R */

#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include<Rmath.h>


SEXP getPhiLogFC(SEXP phi, SEXP contrasts, SEXP exonNames, SEXP includeSims) {

    // 'phi' is the matrix of sims for phi. nRows = sims, nCols = H x J
    // 'contrasts' is the matrix of contrasts
    // 'exonNames' is a vector of the exon names
    // 'includeSims' is an integer where 0==FALSE & any other value==TRUE

    SEXP dimPhi, dimContrasts;
    SEXP tempSims, sims, summary; // Hold the sims & the summary respectively
    SEXP output; // The output
    int H, J, N, M; // Hold the number of conditions, exons, total sims & total contrasts
    int h, j, n, m; // Used for the looping
    int *p_incSims; // Determines whether the sims are to be included in the output or not
    SEXP contNames, summaryColNames, summaryNames, listNames, simNames;
    double *p_sims, *p_summary, *p_contrasts, *p_phi; 
    double curMean=0, curVar=0;
    double nGr0=0, nLs0=0, B; // Hold the numbers greater or less than zero
    int protectCount=0; // Keep track of the number of protected objects
   
    // Protect the input matrices
    PROTECT(phi=AS_NUMERIC(phi)); protectCount++;
    PROTECT(contrasts=AS_NUMERIC(contrasts)); protectCount++;
    PROTECT(exonNames=AS_CHARACTER(exonNames)); protectCount++;
    PROTECT(includeSims=AS_INTEGER(includeSims)); protectCount++;
    p_contrasts = NUMERIC_POINTER(contrasts);
    p_phi = NUMERIC_POINTER(phi);
    p_incSims = INTEGER_POINTER(includeSims);
    
    // The number of parameters & contrasts
    PROTECT(dimPhi = getAttrib(phi,R_DimSymbol)); protectCount++;
    N = INTEGER_POINTER(dimPhi)[0]; // The total number of sims
    PROTECT(dimContrasts = getAttrib(contrasts, R_DimSymbol)); protectCount++;
    H = INTEGER_POINTER(dimContrasts)[0]; // The number of conditions/cellTypes
    J = INTEGER_POINTER(dimPhi)[1] / H;   // The number of exons
    M = INTEGER_POINTER(dimContrasts)[1]; // The number of contrasts

    // Setup the summary object as a list with a separate component for each contrast
    // Each component of the list is a Jx9 matrix
    PROTECT(summary = allocVector(VECSXP, M)); protectCount++;
    for (m=0; m<M; m++) {
	SET_VECTOR_ELT(summary, m, allocMatrix(REALSXP, J, 9)); // Set each list component to a Jx9 matrix
    }

    // Setup the sims object
    // Again each list component is a contrast
    if (p_incSims[0] == 0) { // If sims are not required
	PROTECT(sims = R_NilValue); protectCount++;
    }
    else { // If sims are required
	PROTECT(sims = allocVector(VECSXP, M)); protectCount++;
	PROTECT(simNames =allocVector(VECSXP, 2)); protectCount++; // Setup the colnames for the sims
	SET_VECTOR_ELT(simNames, 0, R_NilValue);
	SET_VECTOR_ELT(simNames, 1, exonNames);
	for (m=0; m<M; m++) {
	    SET_VECTOR_ELT(sims, m, allocMatrix(REALSXP, N, J)); // Set each list component to an NxJ matrix
	    setAttrib(VECTOR_ELT(sims,m), R_DimNamesSymbol, simNames); // Set the dimnames for each matrix
	}
    }

    // Setup the output
    PROTECT(contNames = getAttrib(contrasts, R_DimNamesSymbol)); protectCount++;
    // The column names for the summary matrices
    PROTECT(summaryColNames = allocVector(STRSXP, 9)); protectCount++;
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
    PROTECT(summaryNames = allocVector(VECSXP, 2)); protectCount++;
    SET_VECTOR_ELT(summaryNames, 0, exonNames);
    SET_VECTOR_ELT(summaryNames, 1, summaryColNames);
    // Assign the names
    setAttrib(summary, R_NamesSymbol, VECTOR_ELT(contNames, 1));
    if (p_incSims[0] != 0) {
	setAttrib(sims, R_NamesSymbol, VECTOR_ELT(contNames, 1));
    }
    for (m=0; m<M; m++) {
	setAttrib(VECTOR_ELT(summary,m), R_DimNamesSymbol, summaryNames);
    }
        
    // The actual output
    PROTECT(output = allocVector(VECSXP, 2)); protectCount++;
    SET_VECTOR_ELT(output, 0, summary);
    SET_VECTOR_ELT(output, 1, sims);
    // Set the names for the list
    PROTECT(listNames=allocVector(STRSXP,2)); protectCount++;
    SET_STRING_ELT(listNames, 0, mkChar("summary"));
    SET_STRING_ELT(listNames, 1, mkChar("sims"));
    setAttrib(output, R_NamesSymbol, listNames);

    // Set the sims pointer for the process
    PROTECT(tempSims = allocVector(REALSXP, N)); protectCount++;
    p_sims = NUMERIC_POINTER(tempSims);


    // The process:
    for (m=0; m<M; m++) { // For each contrast
	
	p_summary = NUMERIC_POINTER(VECTOR_ELT(summary,m));

	for (j=0; j<J; j++) { // For each exon

	    // Reset the variables
	    curMean = 0;
	    curVar = 0;
	    
	    for (n=0; n<N; n++) { // For each sim
		
		p_sims[n] = 0; // Reset the value for the tempSims object

		for (h=0; h<H; h++) { // For each condition
		    //p_sims[n] += log(p_phi[h*N + n])*p_contrasts[H*m + h]; // Sample the value for each sim
		    p_sims[n] += log(p_phi[h*J*N + j*N + n])*p_contrasts[H*m + h]; // Sample the value for each sim
		}

		curMean += p_sims[n]/N; // Find the cumulative mean
	    }

	    for (n=0; n<N; n++) {
		curVar += pow(p_sims[n] - curMean, 2) / (N-1); // Find the cumulative sd
	    }

	    // If the sims matrices are required for the output
	    // Each contrast will be a list component with a NxJ matrix
	    if (p_incSims[0] != 0 ) { 
		for (n=0; n<N; n++) {
		    NUMERIC_POINTER(VECTOR_ELT(sims,m))[j*N + n] = p_sims[n];
		}

	    }

	    R_rsort(p_sims, N); // Sort in ascending order

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
		    if (p_sims[n] <=0) nLs0++;
		    if (p_sims[n] >=0) nGr0++;
		}
	    }
	    B = log(((double)nGr0+0.5)/((double)nLs0+0.5)); // Scaled as per the empirical logit transformation

	    // Enter it all in the output 
	    p_summary[j] = curMean; 
	    p_summary[J + j] = sqrt(curVar);
	    p_summary[2*J + j] = p_sims[ imax2(0,(int)(0.025*N-1)) ] ; // The 2.5th Percentile
	    p_summary[3*J + j] = p_sims[ imax2(0,(int)(0.25*N-1)) ] ; // The 25th Percentile
	    p_summary[4*J + j] = p_sims[ imax2(0,(int)(0.5*N-1)) ] ; // The 50th Percentile
	    p_summary[5*J + j] = p_sims[ imax2(0,(int)(0.75*N-1)) ] ; // The 75th Percentile
	    p_summary[6*J + j] = p_sims[ imax2(0,(int)(0.975*N-1)) ] ; // The 97.5th Percentile
	    p_summary[7*J + j] = fmax2((double)nGr0, (double)nLs0)/(double)N;
	    p_summary[8*J + j] = B; // The B statistic

	}
    }

    UNPROTECT(protectCount); // Unprotect all the SEXP objects

    return(output);

}
