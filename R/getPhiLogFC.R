#' @title Sample the log fold-change for each exon from a BMEA.MCMC dataset
#' 
#' @description 
#' Samples the posterior distribution of exon-specific log fold-change for each specified contrast from the given dataset.
#' 
#' @details 
#' Obtains the exon-specific posterior distribution for log fold-change for each supplied contrast 
#' by direct sampling from the MCMC output. 
#' A summary for each contrast with mean, sd & the quantiles are returned in the $summary argument 
#' as the default.
#' 
#' All data is retrieved using C and this function is a simple wrapper to the C function.
#' 
#' @param data an object of class BMEA.MCMC. 
#' Will contain a set of simulations for every condition-specific expression level.
#' @param contr.matrix a contrast matrix with rows representing each condition & columns representing the specific contrasts. 
#' The column names are recycled as rownames in the summary output & for the simulations if kept.
#' @param exonNames a character vector containing the exon (or group) names as specified on the cdf.
#' @param keepSims a Boolean variable that determines whether the sampled values are returned in the output 
#' or just the summary statistics. Defaults to FALSE.
#' 
#' @return A list with the following components:
#' \itemize{
#' \item{summary}{ a list with a separate component for each contrast. 
#' Each contrast-level component is a matrix with the summarised phi log fold-change for each exon within the transcript.}
#' \item{sims}{ a list with a separate component for each contrast. Each contrast-level component is a matrix of the sampled posterior disitributions for each contrast. Rows represent a kept iteration from the MCMC process, and columns represent each exon. Returns NULL if keepSims is set to FALSE.}
#' }
#' 
#' @useDynLib BMEA getPhiLogFC_C
#' 
#' @export
getPhiLogFC <- function(data, contr.matrix, exonNames=NULL, keepSims=FALSE) {
    # A wrapper to the C function

    # Check the input is of the correct form
    if (class(data)!="BMEA.MCMC") return(cat("The data must be of class 'BMEA.MCMC'\n"))

    # Set up the phi matrix
    phi.ind <- which(substr(dimnames(data)[[1]], start=1, stop=3)=="phi")
    phi <- t(matrix(data[phi.ind,,], nrow=length(phi.ind)))
    colnames(phi) <- dimnames(data)[[1]][phi.ind]
    H <- nrow(contr.matrix) # The number of conditions

    # Check the contrast matrix to the best of our ability
    if ((length(phi.ind)/H - as.integer(length(phi.ind)/H))!=0) return(cat("The specified contrasts do not match the number of parameters for phi.\n"))
    if (!is.null(exonNames)) {
        if (!is.character(exonNames)) return(cat("Exon names must be supplied as a character vector.\n"))
        J <- length(exonNames)
        if (length(phi.ind)!=H*J) return(cat("The dimensions of the contrast matrix, exonNames & data do not match.\n"))
    }

    # Get the data as required
    out <- .Call("getPhiLogFC_C", phi, contr.matrix, exonNames, keepSims, PACKAGE="BMEA")

    out

}
