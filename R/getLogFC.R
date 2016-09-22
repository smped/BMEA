#' @title Sample the log fold-change from a BMEA.MCMC dataset
#' 
#' @description 
#' Samples the posterior distribution of log fold-change for each specified contrast from the given dataset.
#' 
#' @details 
#' Obtains the posterior distribution for log fold-change for each contrast by direct sampling from 
#' the MCMC output. 
#' A summary for each contrast with mean, sd & the quantiles are returned in the $summary argument
#' as the default.
#' 
#' All data is retrieved using C and this function is a simple wrapper to the C function.
#' 
#' @param data an object of class BMEA.MCMC. 
#' Will contain a set of simulations for every condition-specific expression level.
#' @param contr.matrix a contrast matrix with rows representing each condition & columns 
#' representing the specific contrasts. 
#' The column names are recycled as rownames in the summary output & for thesimulations if kept.
#' @param keepSims logical variable. 
#' Determines whether the sampled values are returned in the output or just the summary statistics. 
#' Defaults to FALSE.
#' 
#' @return
#' A list with the following components:
#' \itemize{
#' \item{summary}{a summary matrix with mean, standard deviation and the quantiles for each contrast. 
#' Each row represents a different contrast.}
#' \item{sims}{a matrix of the sampled posterior disitributions for each contrast. 
#' Rows represent a kept iteration from the MCMC process, and columns represent each contrast. 
#' Returns NULL if keepSims is set to FALSE}
#' }
#' 
#' @export
getLogFC <- function(data, contr.matrix, keepSims=FALSE) {
    # A wrapper to the C function

    if (class(data)!="BMEA.MCMC") return(cat("The data must be of class 'BMEA.MCMC'\n"))

    mu.ind <- which(substr(dimnames(data)[[1]], start=1, stop=2)=="mu")
    mu <- t(matrix(data[mu.ind,,], nrow=length(mu.ind)))
    colnames(mu) <- dimnames(data)[[1]][mu.ind]

    if (nrow(contr.matrix)!=length(mu.ind)) return(cat("The specified contrasts do not match the number of parameters for mu.\n"))

    out <- .Call("getLogFC", mu, contr.matrix, keepSims, PACKAGE="BMEA")

    out

}
