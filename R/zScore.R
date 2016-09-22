#' @title Calculate the Z-Score for a given transcript/unit
#' 
#' @description Based on the background priors for each probe, 
#' the Z-score is a measure of whether there is any detectable signal above background (DABG)
#' 
#' @details If no signal is present for a group of probes, any observed intensites will be samples 
#' from the distribution for background signal. 
#' Genes with z-score below a pre-defined threshold (e.g. 1.96) can be declared as not expressed, 
#' and likewise for any given exon.
#' 
#' If cacluating the value on a cumulative basis, individual scores will first be calculated for each exon, 
#' and then on an increasing basis.
#' 
#' @param PM the raw observed intensities for the given unit. 
#' Must be a matrix with dimensions nrows=nProbes (K) & ncol= nChips (I).
#' @param lambda the vector of means for the log-normal distribution of the background signal.
#' @param delta the vector of standard deviations for the log-normal distribution of the background signal.
#' @param exons the exon structure of the gene as factors. 
#' If supplied, the function will return a value for each exon.
#' 
#' @return
#' A list with the components:
#' \itemize{
#' \item{$gene}{ the Z-Score for the entire gene}
#' \item{$exons}{ the Z-Score for each exon}
#' }
#' 
#' If the \code{exons} vector is not supplied, the \code{$exons} component will return \code{NULL}.
#' 
#' @examples
#' ## Assuming the PM matrix, and the vectors of means & sds
#' ## have already been obtained
#' PM <- matrix(exp(rnorm(120, 4, 2)), nrow=12, ncol=10) # A matrix of BG only pseudo-PM values
#' lambda <- matrix(4, nrow=12, ncol=10) # E[BG] at each probe
#' delta <- matrix(2, nrow=12, ncol=10) # SD[BG] at each probe
#' exons <- as.factor(rep(paste("Exon",1:3, sep=""), each=4)) # The exon structure
#' z <- zScore(PM, lambda, delta, exons) # The zScore
#' z
#' 
#' @export
zScore <- function (PM, lambda, delta, exons = NULL) {

    ###########
    ## Checks #
    ###########
    if (!is.numeric(PM) || !is.numeric(lambda) || !is.numeric(delta)) stop("All data must be numeric\n")
    if (!is.matrix(PM)) stop("Incorrect Format: PM must be a matrix\n")
    if (min(PM) < 0) stop("The PM matrix cannot contain negative numbers\n")
    if (length(which(is.na(PM))) != 0) stop("The PM matrix cannot contain NA values\n")
    if (missing(lambda)) stop("Missing Data: Means for the priors(lambda) must be provided\n")
    if (missing(delta)) stop("Missing Data: Standard Deviations for the priors (delta) must be provided\n")
    if (length(lambda) != length(delta))  stop("The dimensions of the background priors must be identical\n")
    ###################
    ## End Of Checks ##
    ###################

    ## Set the dimensions
    I <- ncol(PM)
    K <- nrow(PM)
    if (length(lambda) == K) {
        # Expand the bg vectors if necesary
        lambda <- as.vector(rep(lambda, times = I))
        delta <- as.vector(rep(delta, times = I))
    }
    if (length(lambda) != I * K) stop("Incompatible priors: The vector of means does not match the supplied data\n")

    ## Vectorise everything
    lambda <- as.vector(lambda)
    delta <- as.vector(delta)
    PM <- as.vector(PM)

    ## The gene-level value
    Z.gene <- mean((log(PM) - lambda)/delta)/(1/sqrt(I * K))
    Z.exon <- c()

    ## If exon-level is required
    if (!is.null(exons)) {
        ## Checks
        if (length(exons) != K) stop("Incompatible Data: The exons vector does not match the supplied intensities\n")
        if (!is.factor(exons)) stop("Error: The exons must be supplied as factors\n")
        ## Set the dimensions
        J <- length(levels(exons))
        exons <- rep(exons, times = I) # Expand to be a vector which corresponds in size to a vectorised PM matrix
        ind <- lapply(lapply(levels(exons), "==", exons), which) # Get a list of which probe is in which exon
        ## Calculate the scores for each exon separately
        Z.exon <- unlist(
                         lapply(
                                ind,
                                FUN = function(ind, PM, lambda, delta) mean((log(PM[ind]) - lambda[ind])/delta[ind])/(1/sqrt(length(ind))),
                                PM = PM, lambda = lambda, delta = delta
                                )
                         )
        names(Z.exon) <- levels(exons)

    }
    return(list(gene = Z.gene, 
                exons = Z.exon))
}
