#' @title Estimate the parameters for the Background Signal.
#' 
#' @description Fits the background signal model for a set of exon arrays
#' 
#' @details 
#' This function fits the model parameters for background signal estimates, using the set of probes specified in the 'bgProbes' argument. 
#' Currently only the modified MAT model has been implemented.
#' 
#' This function is dependent on the file directory structure as used by aroma.affymetrix. 
#' The .CEL files must be placed in the \code{.../probeData/celSetName/chipType} directory, 
#' and supplied to the function as an \code{AffymetrixCelSet}. 
#' 
#' Quantile normalisation & optical background correction prior to processing with this function 
#' is optional but recommended.
#' 
#' The set of background probes to be used for fitting the model parameters must be located in the \code{.../annotationData/chipTypes/chipType} directory. 
#' The filename must begin with the chipType & end with the suffix as supplied in the bgProbes argument, e.g. "HuEx-1_0-st-v2.r2.genomic.bgp". 
#' The file format must be the same as the .bgp files as supplied by affymetrix, i.e. must have 10 columns, of which one is "probe_sequence".
#'  
#' @param celSet the set of .CEL files being analysed as an \code{AffymetrixCelSet}
#' @param cdf the cdf required for accessing the background probe intensities. 
#'  Some custom cdf files do not contain any background probes.
#' @param bgProbes the tab delimited file containing the sequence data for each background probe.
#' @param method specifies the model for background signal. 
#'  Currently only "MAT" is implemented.
#'  
#' @return 
#'  Returns a list with the following components:
#' \itemize{
#'  \item{$coef}{ the fitted model coefficients}
#'  \item{$fitted}{ the fitted values for each background probe}
#'  \item{$observed}{ the observed values for each background probe}
#'  \item{$method}{ a character string denoting which method was used for model fitting.}
#'  \item{$chipType}{ a character string denoting which chipType (& cdf) was used for fitting the model.}
#' }
#' 
#' @references 
#' H. Bengtsson, K. Simpson, J. Bullard, and K. Hansen, (2008) 
#' \emph{aroma.affymetrix: A generic framework in R for analyzing small to very large Affymetrix data sets in bounded memory}, 
#' Tech Report #745, Department of Statistics, University of California, Berkeley. 
#' 
#' Kapur, K., Xing, Y., Ouyang, Z., Wong, WH. (2007) 
#' \emph{Exon arrays  provide accurate assessments of gene expression} 
#' Genome Biol. 8(5):R82
#' 
#' Johnson, W.E., Li, W., Meyer, C.A., Gottardo, R., Carroll, J.S.,Brown, M., Liu, X.S. (2006) 
#' \emph{Model-based analysis of tiling-arrays for ChIP-chip.}
#' Proc Natl Acad Sci USA 103:12457-12462
#' 
#' @import aroma.affymetrix
#' @import affy
#' 
#' @export
fitBackgroundParameters <- function(celSet, cdf=NULL, bgProbes="r2.genomic.bgp", method="MAT") {

  # The checks in this process need to be checked...

      # Check that the celSet is actually a celSet
    if(class(celSet)[1]!="AffymetrixCelSet") stop("Incorrect data format: The input must be of Class: AffymetrixCelSet\n")

    # Get the cdf, or set as required
    if (is.null(cdf)) {
        cdf <- getCdf(celSet)
    }
    else {
        if (length(class(cdf))!=11 && class(cdf)[1]!="AffymetrixCdfFile") stop("The specified cdf file is not actually a cdf!\n")
        setCdf(celSet,cdf)
    }

    parentName <- getParentName(celSet) # Get the name of the dataset/celSet
    chipType <- getChipType(celSet)
    wd <- getwd() # The working directory

    # Check that a valid method has been selected
    if (method!="MAT") stop("Invalid Model Name: The bg model can currently only be MAT\n")

    # Check that the bgProbes file exists & then load it
    bgFile <- paste(wd, "annotationData", "chipTypes", chipType, paste(chipType, bgProbes, sep="."), sep="/")
    if(!file.exists(bgFile)) stop("Cannot find the specified bgProbes file:\n", bgFile, "\n")
    bgp <- read.table(bgFile, skip=7, header=TRUE, stringsAsFactors=FALSE) # Load the background probes

    # Check the format of the loaded bgp object. It should have 10 columns with one named "probe_sequence"
    if(ncol(bgp)!=10) stop("Invalid bgProbes File: Please check the file format. It should contain 10 columns\n")
    if(is.na(match("probe_sequence",colnames(bgp)))) stop("Invalid bgProbes File: File must have a column titled 'probe_sequence'\n")
    if(is.na(match("x",colnames(bgp)))) stop("Invalid bgProbes File: File must have a column named 'x'\n")
    if(is.na(match("y",colnames(bgp)))) stop("Invalid bgProbes File: File must have a column named 'y'\n")
    if(is.na(match("gc_count",colnames(bgp)))) stop("Invalid bgProbes File: File must have a column named gc_count\n")
    nProbes <- nrow(bgp) # The number of probes in the bgp file
    nr <- getDimension(cdf)[1] # The number of rows on the cdf
    indices <- affy::xy2indices(bgp$x, bgp$y, nr=nr) # Get the cell indices

    # Get the other required information
    nChips <- nbrOfArrays(celSet) # Get the number of arrays
    bgIntensities <- extractMatrix(celSet, indices)  # Get the intensity data for the set of probes

    if (method=="MAT") {

        # Define the model matrix as a vector for each probe sequence. This will be nProbes*80
        X <- setMatMatrix(as.character(bgp[,"probe_sequence"]))
        colnames(X) <- c("alpha", paste("beta.A", 1:25, sep=""), paste("beta.C", 1:25, sep=""), paste("beta.G", 1:25, sep=""), "gamma.A", "gamma.C", "gamma.G", "gamma.T")
        rownames(X) <- as.character(bgp[,1])

        # Fit for the entire dataset
        fullX <- matrix(rep(t(X), times=nChips), nrow=nChips*nProbes, byrow="T") # Expand X for the complete dataset
        matCoef <- solve(t(fullX)%*%fullX)%*%t(fullX)%*%log(as.vector(bgIntensities)) # Find the coefficients
        matCoef <- matrix(matCoef, ncol=1)
        colnames(matCoef) <- parentName
        fitted <- X%*%matCoef

        rownames(matCoef) <- colnames(X)
        rownames(fitted) <- rownames(X)

    }

    # Change the output to be a list of coefficients, fittedValues & observed bgIntensities
    list(coef=matCoef, 
         fitted=fitted, 
         observed=log(bgIntensities), 
         method=method, 
         chipType=chipType)

}
