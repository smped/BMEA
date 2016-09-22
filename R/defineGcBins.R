#' @title Estimate the mean and standard deviation for the log background signal
#' 
#' @description Based on the supplied set of background probes, 
#' estimate the background mean and standard deviation for a given GC count.
#' 
#' @param celSet the set of arrays being processed, as an \code{AffymetrixCelSet}.
#' @param cdf an alternate cdf. 
#' This will override the one specified in \code{celSet}. 
#' If using a custom Affymetrix cdf, 
#' the original Affymetrix cdf should be specified here as the custom one may not contain the set of background probes.
#' @param bgProbes the tab delimited file containing the sequence data for each background probe.
#' @param nBins the number of bins to form. Defaults to 20
#' @param fitIndividual logical variable. 
#' If \code{fitIndividual=TRUE}, the parameters will be estimated on an array-specific basis, 
#' otherwise the parameters will be estimated across the entire celSet.
#' 
#' @details 
#' This function finds the log scale mean (lambda) and standard deviation (delta) of the specified \code{bgProbes} for each \code{gc_count}.
#' 
#' @return 
#'   A list with the following components:
#'   \itemize{ 
#'   \item{lambda}{the mean of the log-transformed \code{bgProbes} for a given \code{gc_count}}
#'   \item{delta}{the standard deviation of the log-transformed \code{bgProbes} for a given \code{gc_count}}
#'   }
#'   For easy downstream usage, the names/rownames for each component are specified as "GC4" through to "GC25."
#'   
#' @import aroma.affymetrix
#' 
#' @export
defineGcBins <- function(celSet, cdf=NULL, bgProbes="r2.genomic.bgp", nBins=20, fitIndividual=TRUE) {

    # Check that the celSet is actually a celSet
    if(class(celSet)[1]!="AffymetrixCelSet") stop("Incorrect data format: The input must be of Class: AffymetrixCelSet\n")
    parentName <- getParentName(celSet) # Get the name of the dataset/celSet
    chipType <- getName(getCdf(celSet))
    wd <- getwd() # The working directory

    # Check the cdf & set the correct one
    origCdf <- getCdf(celSet)
    if (is.null(cdf)) { # If no cdf is specified, extract the info
        cdf <- origCdf
    }
    else { # Otherwise set the the one that has been specified manually
        if (length(class(cdf))!=11 && class(cdf)[1]!="AffymetrixCdfFile") stop("The specified cdf file is not actually a cdf!\n")
        setCdf(celSet,cdf)
    }

    # Check that the bgProbes file exists & then load it
    bgFile <- paste(wd,"annotationData","chipTypes",chipType,paste(chipType,bgProbes,sep="."),sep="/")
    if(!file.exists(bgFile)) stop("Cannot find the specified bgProbes file:\n",bgFile,"\n")
    bgp <- read.table(bgFile,skip=7,header=TRUE, stringsAsFactors=FALSE) # Load the background probes

    # Check the format of the loaded bgp object. It should have 10 columns with one named "probe_sequence"
    if(ncol(bgp)!=10) stop("Invalid bgProbes File: Please check the file format. It should contain 10 columns\n")
    if(is.na(match("probe_sequence",colnames(bgp)))) stop("Invalid bgProbes File: File must have a column named probe_sequence\n")
    if(is.na(match("x",colnames(bgp)))) stop("Invalid bgProbes File: File must have a column named 'x'\n")
    if(is.na(match("y",colnames(bgp)))) stop("Invalid bgProbes File: File must have a column named 'y'\n")
    if(is.na(match("gc_count",colnames(bgp)))) stop("Invalid bgProbes File: File must have a column named gc_count\n")
    nProbes <- nrow(bgp) # The number of probes in the bgp file
    nr <- getDimension(cdf)[1] # The number of rows on the cdf
    indices <- xy2indices(bgp$x, bgp$y, nr=nr) # Get the cell indices

    # Get the other required information
    nChips <- nbrOfArrays(celSet) # Get the number of arrays
    fileNames <- getFullNames(celSet) # The arrayNames

    # Get the intensity data for the set of probes. The first col of the bgProbes must contain the cell index
    bgIntensities <- extractMatrix(celSet, indices)
    gc_count <- as.integer(bgp$gc_count) # This will range from 4 to 23
    levels <- unique(gc_count)
    if (!is.null(nBins)){ # If manually setting the number of bins
        if (nBins < length(levels)) { # Only change if reducing the number of bins
            levels <- as.integer(seq(median(levels)-nBins/2, by=1, length.out=nBins)) # This will give the vector of integers
        }
    }
    levels <- levels[order(levels)]

    if (fitIndividual) {

        lambda <- delta <- matrix(NA, nrow=length(levels), ncol=nChips)
        # Do the first bin separately as it needs a '<=' rather than ==
        ind <-  which(gc_count<=levels[1])
        lambda[1,] <- colMeans(log(bgIntensities[ind,]))
        delta[1,] <- colSds(log(bgIntensities[ind,]))
        # The next lot of bins can be done in a loop
        for (i in 2:(length(levels)-1)) {
            ind <- which(gc_count==levels[i])
            lambda[i,] <- colMeans(log(bgIntensities[ind,]))
            delta[i,] <- colSds(log(bgIntensities[ind,]))
        }
        # The final bin needs a >=
        ind <-  which(gc_count>=levels[length(levels)])
        lambda[length(levels),] <- colMeans(log(bgIntensities[ind,]))
        delta[length(levels),] <- colSds(log(bgIntensities[ind,]))

        rownames(lambda) <- rownames(delta) <- paste("GC",levels,sep="")
        colnames(lambda) <- colnames(delta) <- fileNames

    }
    else {

        lambda <- delta <- vector("numeric", length(levels))
        ## Do the first bin separately as it needs a '<=' rather than ==
        ind <-  which(gc_count<=levels[1])
        lambda[1] <- mean(as.vector(log(bgIntensities[ind,])))
        delta[1] <- sd(as.vector(log(bgIntensities[ind,])))
        for (i in 2:(length(levels)-1)) {
            ind <- which(gc_count==levels[i])
            lambda[i] <- mean(as.vector(log(bgIntensities[ind,])))
            delta[i] <- sd(as.vector(log(bgIntensities[ind,])))
        }
        ind <-  which(gc_count>=levels[length(levels)])
        lambda[length(levels)] <- mean(as.vector(log(bgIntensities[ind,])))
        delta[length(levels)] <- sd(as.vector(log(bgIntensities[ind,])))

        names(lambda) <- names(delta) <- paste("GC",levels,sep="")
    }

    # Reset the cdf to the original one as this will affect the celSet outside this function
    setCdf(celSet, origCdf)

    list(lambda=lambda, delta=delta)

}
