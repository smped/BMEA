#' @title Define the bins, mean & sd for log background signal using the MAT model.
#' 
#' @description Define the boundaries for the background signal, 
#' based on the fitted MAT values for each probe, and the subsequent means & standard deviations.
#' 
#' @param bgParam the output from the function \code{fitBackgroundParameters}. 
#' Must contain the slot \code{$coef}.
#' @param nBins specifies the number of bins to form based on the fitted MAT values. 
#' Defaults to 20 bins.
#' @param fitIndividual logical.Determines whether bin means & sds are estimated for each array or for the entire dataset.
#' 
#' @details Given the fitted MAT values from the background set of probes in the \code{bgParam} argument, 
#' the boundaries for forming the background prior distributions are calculated using this function. 
#' Bins are formed using approximately equal numbers of probes.
#' 
#' @return 
#'   Returns a list with the components: 
#'   \itemize{
#'   \item{$boundaries}{the boundaries for the bins. The lower boundary is implied as zero & 
#'   the upper boundary is implies as \code{Inf}.}
#'   \item{$lambda}{the observed mean of the background probes classified into each bin}
#'   \item{$delta}{ the observed standard deviation of the background probes as classified into each bin}
#'   }
#'   
#'   All means and standard deviations are on the log-scale and refer to the normal distributions for log(Background).
#' 
#' @references 
#'  Kapur, K., Xing, Y., Ouyang, Z., Wong, WH. (2007) \emph{Exon arrays provide accurate assessments of gene expression} 
#'  Genome Biol. 8(5):R82 
#'   
#'  Johnson, W.E., Li, W., Meyer, C.A., Gottardo, R., Carroll, J.S., Brown, M., Liu, X.S. (2006) \emph{Model-based analysis of tiling-arrays for ChIP-chip.}
#'  Proc Natl Acad Sci USA 103:12457-12462
#' 
#' @seealso \code{\link{fitBackgroundParameters}}
#' 
#' @export
defineMatBins <- function(bgParam, nBins=20, fitIndividual=TRUE) {

  # The checks still need to be checked for this function

    # Check that bgParam is the correct output structure
    if(is.null(bgParam$fitted)) stop("Incorrect background: The object bgParam must contain fitted values\n")
    if(is.null(bgParam$observed)) stop("Incorrect background: The object bgParam must contain observed logIntensities\n")

    nChips <- ncol(bgParam$observed)
    quant <- seq(0,1,length.out=nBins+1) # The quantiles being sought

    # Set the output matrices & bin boundaries
    lambda <- delta <- matrix(nrow=nBins,ncol=nChips)
    boundaries <- quantile(bgParam$fitted,probs=quant) # Get the boundary points
    boundaries[c(1, nBins+1)] <- c(0, 16*log(2)) # Set as the min & max of all possible values

    if (fitIndividual) { # If fitting each array separately

        for (i in 1:nBins) {
            ind <- intersect(which(bgParam$fitted > boundaries[i]), 
                             which(bgParam$fitted <= boundaries[i+1])) # Find the probes for each bin
            lambda[i,] <- colMeans(bgParam$observed[ind,])
            delta[i,] <- matrixStats::colSds(bgParam$observed[ind,])
        }

        colnames(lambda) <- colnames(delta) <- colnames(bgParam$observed)
    }
    else { # If fitting the entire dataset

         for (i in 1:nBins) {
             ind <- intersect(which(bgParam$fitted > boundaries[i]), 
                              which(bgParam$fitted <= boundaries[i+1])) # Find the probes for each bin
             lambda[i,] <-mean(as.vector(bgParam$observed[ind,]))
             delta[i,] <- sd(as.vector(bgParam$observed[ind,]))
         }

    }

    list(boundaries=boundaries[-c(1,nBins+1)],
         lambda=lambda, 
         delta=delta)

}
