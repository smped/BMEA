#' @title Summarise the output of the BMEA MCMC process
#' 
#' @description Takes the kept simulations from a BMEA MCMC process and provides the summary statistics.
#' 
#' @details 
#' Provides summary statistics, quantiles & convergence statistics for a BMEA.MCMC object.
#' 
#' To approximate normality, sigma_mu, sigma_p & sigma_S are log transformed when calculating rHat. 
#' If \code{transform = TRUE}, all S simulations are log transformed & phi simulations are logit transformed.
#' If the mixture prior for phi has been used, this should be set to \code{FALSE} as values of zero 
#' or one will be possible & will not be suitable for transformation.
#' 
#' @param data an object of class BMEA.MCMC
#' @param transform logical. If \code{transform = TRUE}, the S & phi simulations will be transformed to the 
#' log & logit scales respectively when calculating the convergence statistics.
#' 
#' @return A matrix with columns containing the mean, sd, quantiles (0.025, 0.25, 0.50, 0.75, 0.975), rHat & 
#' nEff for each BMEA model parameter.
#' 
#' @references 
#' Brooks, S.P. and Gelman, A. (1998) 
#' \emph{General Methods for Monitoring Convergence of Iterative Simulations} 
#' Journal of Computational and Graphical Statistics, Vol. 7, No. 4 (Dec), pp. 434-455  
#' 
#' @useDynLib BMEA summariseMCMC
#' 
#' @export 
summariseChains <- function(data, transform=TRUE) {
    # data must be the output from a BMEA run, with class BMEA.MCMC
    # transform indicates whether to logit transform phi & log transform S when calculating convergence statistics
    # returns a matrix with the summary statistics

    if (class(data)!="BMEA.MCMC") stop("The object must be of class 'BMEA.MCMC'\n")

    out <- .Call("summariseMCMC",data, transform, PACKAGE="BMEA") # TRUE indicates logit transformation for phi, FALSE will disable this feature

    return(out)

}
