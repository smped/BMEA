#' @title Fit the BMEA model for a single gene.
#' 
#' @description Fits the BMEA model using the MCMC process for a single gene.
#' 
#' @details 
#' This function runs the complete BMEA process for a single gene. 
#' Currently the only implemented model is the uniform prior for the exon proportions term (phi). 
#' 
#' Computation times will vary to a large degree based on the processor involved. 
#' It should be noted that the process for a single gene can take in the order of minutes. 
#' 
#' The MCMC process itself is wrtten in C and this function is a wrapper to this function.
#' 
#' @param PM the raw intensity data for a single gene.
#' @param conditions a vector of factors defining which array belongs to which cell-type, or condition.
#' @param exons a vector of factors defining which probe belongs to which exon (or group).
#' @param lambda a matrix of the same dimensions as the PM matrix. 
#' Contains the probe-specific means for the log-normal distributions used as the priors for the background signal.
#' @param delta a matrix of the smae dimensions as the PM matrix. 
#' Contains the standard deviations for the log-normal distributions used as the priors for the background signal.
#' @param mcmcParam a list defining the parameters for the MCMC process. 
#' Requires the arguments: 
#'  \itemize{
#'  \item{$nChains}{ the number of chains to run. Defaults to 3 if not specified.}
#'  \item{$nIter}{ the number of iterations to run. Defaults to 12,000 if not specified.}
#'  \item{$nBurnin}{ the number of burnin iterations. Defaults to the first half of the iterations.}
#'  \item{$nThin}{ the number of iterations used to thin the stored data, 
#'  e.g. if \code{$nThin} = 10, every 10th iteration will be kept. Defaults to 6 if not specified.}
#'  }
#' @param model a character string defining which BMEA model to run. 
#'  Currently only the full model (\code{model="Full"}) or one with no probeEffects 
#'  (\code{model="noProbes"}) are implemented.
#' @param ... not used
#' 
#' @return 
#' Returns a 3D-array of class BMEA.MCMC. 
#' The first dimension contains the model parameters for the given gene. 
#' The second dimension stores the iterations for each chain, 
#' whilst the third dimension contains the number of chains used for the process.
#' 
#' The names for the first dimension are vital for downstream analysis &  
#' will specify which parameters are mu & phi.
#' 
#' @useDynLib BMEA runSingleExonUniformMCMC
#' @useDynLib BMEA runUniformMCMC
#' @useDynLib BMEA runUniformMCMCnoProbe
#' 
#' @export
runMCMC.BMEA <- function(PM, conditions, exons, lambda, delta, mcmcParam=NULL, model="Full",...) {

    # Wrapper for the C function that runs the BMEA process

    # The output is a 3D array of dimension nParam x nSims x nChains. The class of BMEA.MCMC will be assigned

    # Checks on the expression matrix
    I <- length(conditions)
    PM <- matrix(PM, ncol=I)
    lambda <- matrix(lambda, ncol=I)
    delta <- matrix(delta, ncol=I)
    K <- nrow(PM)
    if (min(PM)<0) stop("The PM matrix cannot contain negative values\n")
    if (max(PM)>2^16) stop("The PM matrix contains a value greater than the Affymetrix scanner resolution\n")
    if (length(which(is.na(PM)))!=0) stop("The PM matrix cannot contain NA values")

    # Checks on the conditions vector
    if (!is.factor(conditions)) stop("The conditions vector must be supplied as factors\n")

    # Checks on the exons vector
    if (length(exons)!=K) stop("The exons vector does not match the number of probes\n")
    if (!is.factor(exons)) stop("The exons vector must be supplied as factors\n")
    J <- length(levels(exons)) # The number of exons

    # Checks on the lambda matrix
    if (nrow(lambda)!=K) stop("The matrix of background means (lambda) has an incorrect number of probes\n")

    # Checks on the delta matrix
    if (nrow(delta)!=K) stop("The matrix of background means (delta) has an incorrect number of probes\n")

    # Check & set up the MCMC parameters
    if (is.null(mcmcParam)) { # If none are supplied
        mcmcParam <- vector("list",4)
        mcmcParam[[1]] <- as.integer(3)
        mcmcParam[[2]] <- as.integer(16000)
        mcmcParam[[3]] <- as.integer(8000)
        mcmcParam[[4]] <- as.integer(8)
        names(mcmcParam) <- c("nChains","nIter","nBurnin","nThin")
    }
    else {
        if (!is.list(mcmcParam)) stop("The MCMC parameters must be supplied as a list\n")
        if (is.null(mcmcParam$nChains)) {
            mcmcParam$nChains <- as.integer(3) # The default will be set as 3 if no value is specified
            message(sprintf("No value supplied for $nChains: set to the default of %i.\n", mcmcParam$nChains))
        }
        else {
            if (!is.numeric(mcmcParam$nChains)) stop("The number of chains must be supplied as an number\n")
            mcmcParam$nChains <- as.integer(mcmcParam$nChains)
        }
        if (is.null(mcmcParam$nIter)) {
            mcmcParam$nIter <- as.integer(16000) # The default will be set as 16000 if no value is specified
            message("No value supplied for $nIter: set to the default of 16000.\n")
        }
        else {
            if (!is.numeric(mcmcParam$nIter)) stop("The number of iterations must be supplied as an number\n")
            mcmcParam$nIter <- as.integer(mcmcParam$nIter)
        }
        if (is.null(mcmcParam$nBurnin) || !is.numeric(mcmcParam$nBurnin) ) {
            mcmcParam$nBurnin <- as.integer(mcmcParam$nIter/2) # The default will be set as nIter/2 if no value is specified
            message(sprintf("Invalid/Missing value for nBurnin: set to the default %i.\n",mcmcParam$nBurnin))
        }
        if (is.null(mcmcParam$nThin) ||!is.numeric(mcmcParam$nThin)) {
            mcmcParam$nThin <- as.integer(8) # Set the default to 8
            message(sprintf("Invalid/Missing value for nThin: set to the default %i.\n",mcmcParam$nThin))
        }
        mcmcParam$nBurnin <- as.integer(mcmcParam$nBurnin)
        mcmcParam$nThin <- as.integer(mcmcParam$nThin)
    }
    # Finally, ensure that the mcmcParam arguments are in the correct order:
    ord <- c(which(names(mcmcParam)=="nChains"), which(names(mcmcParam)=="nIter"), which(names(mcmcParam)=="nBurnin"), which(names(mcmcParam)=="nThin"))
    mcmcParam <- mcmcParam[ord]

    # Check the model is permissable
    if (!tolower(model) %in% c("full","noprobes")) stop("Invalid Model Specified: Only 'Full' or 'noProbes' are permissable.\n")

    # So if the checks are all OK, the process can proceed
    if (J==1) {
        sims <- .Call("runSingleExonUniformMCMC", PM, conditions, exons, lambda, delta, mcmcParam, PACKAGE="BMEA")
    }
    else {
        if (tolower(model)=="full") sims <- .Call("runUniformMCMC", PM, conditions, exons, lambda, delta, mcmcParam, PACKAGE="BMEA")
        if (tolower(model)=="noprobes") sims <- .Call("runUniformMCMCnoProbe", PM, conditions, exons, lambda, delta, mcmcParam, PACKAGE="BMEA")
    }

    return(sims)

}
