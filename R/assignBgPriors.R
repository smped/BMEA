#' @title Assign the Background Signal Priors for each PM probe
#' 
#' @description The function assigns a log-normal prior distribution for background signal to each PM probe on the specified cdf. 
#' The prior is estimated based on sequence composition using either MAT or the GC count
#' 
#' @param celSet the set of CEL files to be analysed
#' @param seqFile the file with the sequence data for all PM probes. 
#' Must contain the xy co-ordinates of each probe & sequence data for each probe.
#' @param bgBins the output of either \code{defineMatBins} or  \code{defineGcBins}.
#' @param bgParam the model parameters as returned by the function \code{fitBackgroundParameters}. 
#' Not required if using \code{method="GC"}.
#' @param method only "MAT" & "GC" can be supplied. Defaults to "MAT".
#' @param path the root directory where files should be written. 
#' Defaults to the working directory
#' @param batchSize the number of probes to be written to file in batches. 
#' Calculated as a function of memory size.
#' @param overWrite \code{logical}. 
#' Determines whether to over-write any existing .CEL files with new information. 
#' If \code{overWrite = FALSE}, checks for pre-existing celFiles will be performed & 
#' if a complete set matching the supplied \code{celSet} is found, this will be used.
#' 
#' @details
#' This process will create two sets of .CEL files containing the mean (lambda) & standard deviation (delta) 
#' for the background priors under the BMEA model for each PM probe contained on the supplied cdf. 
#' The priors are based on the assumption that log(B) ~ N(lambda, delta).
#' 
#' The .CEL files will be written to the subdirectories of the path or working directory  \code{backgroundPriors/parentName,method,lambda/chipType} &
#' \code{backgroundPriors/parentName,method,delta/chipType}.
#'  
#' @return 
#'  Returns a list with components \code{$lambda} & \code{$delta}.
#'  Each component is an \code{AffymetrixCelSet} containing the respective means & sds.
#'  
#' @references 
#'  H. Bengtsson, K. Simpson, J. Bullard, and K. Hansen, (2008) \emph{aroma.affymetrix: A generic framework in R for analyzing small to very large Affymetrix data sets in bounded memory}, 
#'  Tech Report #745, Department of Statistics, University of California, Berkeley. 
#'   
#'  Kapur, K., Xing, Y., Ouyang, Z., Wong, WH. (2007) \emph{Exon arrays provide accurate assessments of gene expression} 
#'  Genome  Biol. 8(5):R82
#'   
#'  Johnson, W.E., Li, W., Meyer, C.A., Gottardo, R., Carroll, J.S., Brown, M., Liu, X.S. (2006) \emph{Model-based analysis of tiling-arrays for ChIP-chip.} 
#'  Proc Natl Acad Sci USA 103:12457-12462 
#' 
#' @export
assignBgPriors <- function(celSet, seqFile, bgBins=NULL, bgParam=NULL, method="MAT", path, 
                           batchSize=NULL, overWrite=FALSE) {

  # verbose should be incorporated into this process as well. Currently the default will print out everything
  
    # Check the method:
    method <- toupper(method) # Set it to upper case
    if (!(method=="MAT" || method=="GC")) stop("Background methods can only be MAT or GC.\n")

    # Check the celSet:
    if(class(celSet)[1]!="AffymetrixCelSet") stop("Incorrect data format: The input must be of Class: AffymetrixCelSet\n")
    # Get the celSet info
    celPath <- file.path(getwd(), getPath(celSet)) # Where the celFiles are
    celNames <- getNames.GenericDataFileSet(celSet) # The names of the files in the celSet
    nChips <- nbrOfArrays(celSet)
    parentName <- getParentName(celSet)
    chipType <- getChipType(celSet)

    # Set up the final paths & create them if needed
    if (missing(path)) {
        lambdaPath <- file.path("backgroundPriors", paste(parentName, method, "lambda", sep=","), chipType) # Where the lambda files should be
        deltaPath <- file.path("backgroundPriors", paste(parentName, method, "delta", sep=","), chipType) # Where the delta files should be
    }
    else {
        lambdaPath <- file.path(path, "backgroundPriors", paste(parentName, method, "lambda", sep=","), chipType) # Where the lambda files should be
        deltaPath <- file.path(path, "backgroundPriors", paste(parentName, method, "delta", sep=","), chipType) # Where the delta files should be
    }
    if (!file.exists(lambdaPath)) dir.create(lambdaPath, recursive=TRUE) # Create the directory if it doesn't exist
    if (!file.exists(deltaPath)) dir.create(deltaPath, recursive=TRUE) # Create the directory if it doesn't exist

    # Set up some control variables
    lambdaDone <- FALSE
    deltaDone <- FALSE

    if (!overWrite) { # If not overwriting, check to see if the celSets already exist. If so, use them for the output

        # Lambda
        message("Checking for pre-existing lambda CEL files...")
        lambdaNames <- list.files(lambdaPath)
        if (length(lambdaNames)>=nChips){ # If there are enough files in the directory
            # If the celNames match every file in the celSet, they already exist & shouldn't be overwritten
            if (length(which(toupper(paste(celNames, ".CEL", sep="")) %in% toupper(lambdaNames)))==nChips) {
                lambdaDone <- TRUE
                message("files found!\n")
                lambda <- AffymetrixCelSet$byPath(lambdaPath, checkChipType=FALSE)
            }
        }
        if (!lambdaDone) message("none found.\n")
        flush.console()

        #Delta
        message("Checking for pre-existing delta CEL files...")
        deltaNames <- list.files(deltaPath)
        if (length(deltaNames)>=nChips){ # If there are enough files in the directory
            # If the celNames match every file in the celSet, they already exist & shouldn't be overwritten
            if (length(which(toupper(paste(celNames,".CEL", sep="")) %in% toupper(deltaNames)))==nChips) {
                deltaDone <- TRUE
                message("files found!\n")
                delta <- AffymetrixCelSet$byPath(deltaPath, checkChipType=FALSE)
            }
        }
        if (!deltaDone) message("none found.\n")
        flush.console()

        if (lambdaDone && deltaDone) { # If both sets of CEL files exist

            message("Both background prior celSets exist & were not overwritten.\n")
            # Return a list with the two celSets
            out <- list(lambda=lambda, delta=delta)
            class(out) <- "AffymetrixCelSetList"
            return(out)

        }
    }

    #####################################################################
    # If pre-existing files are not used the function will keep running #
    #####################################################################

    # Check the bgBins & bgParam objects
    if (is.null(bgBins)) stop("The bgBins object must be supplied.\n")
    if (length(which(!c("lambda","delta") %in% names(bgBins)))!=0) stop("The bgBins object must contain the components lambda & delta.\n")
    if (length(bgBins$lambda) != length(bgBins$delta)) stop("The number of bins must be the same for lambda & delta.\n")
    if (method=="GC") {
        if (is.matrix(bgBins$lambda) && nrow(bgBins$lambda) != 20) stop("The bgBins object should contain 20 bins if using the GC method.\n")
        if (is.vector(bgBins$lambda) && length(bgBins$lambda) != 20) stop("The bgBins object should contain 20 bins if using the GC method.\n")
        nBins <- 20
    }
    else {
        # The bgParam argument also needs to be checked if GC bins are not being used
        if (is.null(bgParam)) stop("The bgParam object must be supplied if not using GC bins for the priors.\n")
        if (is.null(bgParam$coef)) stop("The bgParam object must contain model co-efficients.\n")
        if (toupper(bgParam$method)!=method) stop(sprintf("Error: The model used for fitting the coefficients (%c) is not the same the requested method (%c).\n", bgParam$method, method))
        if (ncol(bgParam$coef)!=1) stop("The model coefficients must be a matrix with a single column.\n")
        if (is.null(bgBins$boundaries)) stop("Bin boundaries must be specified in the bgBins object.\n")
        nBins <- length(bgBins$boundaries) + 1
        if (nBins != (nrow(bgBins$lambda))) stop("Bin boundaries & means do not match.\n")
        if (nBins != (nrow(bgBins$delta))) stop("Bin boundaries & standard deviations do not match.\n")
    }
    if (is.null(bgBins$lambda)) stop("Bin specific means must be supplied in the bgBins object.\n")
    if (is.null(bgBins$delta)) stop("Bin specific standard deviations must be supplied in the bgBins object.\n")
    if (!is.vector(bgBins$lambda) && ncol(bgBins$lambda)!=nChips) stop("The model coefficients do not match the number of arrays in the celSet.\n")

    # Check that the bgBins object is compatible with the celSet
    if (sum(dim(bgBins$delta)- dim(bgBins$lambda))!=0) stop("The dimensions of lambda & delta do not match.\n")
    # If the dimensions are not either 1 or nChips
    if (!is.vector(bgBins$lambda) && ncol(bgBins$lambda)!=nChips) stop("The dimensions of lambda & delta are not suitable for the number of chips.\n")
    # Modify the bgBins objects to ensure that ncol=nChips. If this is already the case, no change will occur
    if (is.vector(bgBins$lambda)) binNames <- names(bgBins$lambda)
    if (is.matrix(bgBins$lambda)) binNames <- rownames(bgBins$lambda)
    if (ncol(bgBins$lambda)==1) {
        bgBins$lambda <- matrix(bgBins$lambda, nrow=nBins, ncol=nChips)
        bgBins$delta <- matrix(bgBins$delta, nrow=nBins, ncol=nChips)
        rownames(bgBins$lambda) <- rownames(bgBins$delta) <- binNames
    }

    # Create the CEL files that hold the priors.
    message("Background Priors: Creating CEL files...")
    flush.console()
    for (i in 1:nChips) {

        # Read the header from the files in the celSet
        hdr <- readCelHeader(file.path(celPath, paste(celNames[i], "CEL", sep=".")))
        # Create the current file path for lambda & delta
        lambdaFile <- file.path(lambdaPath, paste(celNames[i], "CEL", sep="."))
        deltaFile <- file.path(deltaPath, paste(celNames[i], "CEL", sep="."))

        if (!overWrite) { # If not overwriting any existing files
            # Only write the file if it doesn't exist already
            if (!file.exists(lambdaFile) && !lambdaDone) {
                hdr$filename <- lambdaFile
                createCel(lambdaFile, hdr)
            }
            if (!file.exists(deltaFile) && !deltaDone) {
                hdr$filename <- deltaFile
                createCel(deltaFile, hdr)
            }
        }
        else { # Otherwise, just create new ones
            if (!lambdaDone) {
                hdr$filename <- lambdaFile
                createCel(lambdaFile, hdr, overwrite=TRUE)
            }
            if (!deltaDone) {
                hdr$filename <- deltaFile
                createCel(deltaFile, hdr, overwrite=TRUE)
            }
        }

    }
    message("done\n")

    # seqFile is the name of the sequence file. This must be in the .../annotationData/chipTypes/chipType folder
    cdf <- getCdf(celSet) # Get the cdf
    seqFile <- file.path("annotationData", "chipTypes", getName(cdf), basename(seqFile)) # include the directory info
    if(!file.exists(seqFile)) stop("The specified sequence file cannot be found.\n")

    # Read in the sequence data
    message("Reading probe sequence data...")
    flush.console()
    # Check the colnames of seqFile
    seqFileCols <- tolower(scan(seqFile, nlines=1, sep="\t", what="character"))
    xCol <- intersect(grep("x", seqFileCols), grep("probe", seqFileCols)) # Find the column with the X co-ordinate
    if (length(xCol)!=1) stop(sprintf("Unable to locate the probe X co-ordinates in %s\n",seqFile))
    yCol <- intersect(grep("y", seqFileCols), grep("probe", seqFileCols)) # The column with the Y co-ordinate
    if (length(yCol)!=1) stop(sprintf("Unable to locate the probe Y co-ordinates in %s\n",seqFile))
    seqCol <- grep("sequence", seqFileCols) # The column with the sequence data
    if (length(seqCol)!=1) stop(sprintf("Unable to locate the sequence data in %s\n",seqFile))
    seqData <- read.table(seqFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    cdfCols <- nbrOfColumns(cdf)
    seqCells <- affy::xy2indices(as.integer(seqData[,xCol]), as.integer(seqData[,yCol]), nc=cdfCols)
    seqData <- data.frame(seqCells, seqData[,seqCol], stringsAsFactors=FALSE)
    colnames(seqData) <- c("Probe ID","probe sequence")
    message("done\n")
    # seqData will only have the columns "Probe ID" & "probe sequence"

    # Check that every cell on the cdf is represented by some sequence data
    cellInd <- unlist(getCellIndices(cdf)) # The cell indices on the cdf
    if (length(which(!(cellInd %in% seqData[,"Probe ID"])))!=0) stop("The probe sequence file doesn't match the cdf.\n")
    # As the process is set by cycling through the seqData object, ensure that this exactly matches the cdf
    seqData <- seqData[match(cellInd,seqData[,"Probe ID"]),]

    # Set up the index vector as the length of probes being considered.
    nProbes <- length(cellInd)
    ind <- 1:nProbes

    # Set the batchSize depending on the available memory:
    if (is.null(batchSize)) {
        if (memory.limit()<=2048) {
            batchSize <- as.integer(2^19/nChips)
        }
        else {
            if (memory.limit()<=4096) {
                batchSize <- as.integer(2^20/nChips)
            }
            else   batchSize <- as.integer(2^21/nChips)
        }
    }

    # Set the batches
    batch <- as.integer(0)
    nBatches <- ceiling(nProbes/batchSize)

    # Now fit the sequences in batches & write to the .CEL files
    # If using the MAT model
    if (method=="MAT") {

        while (length(ind)>0) {

            batch <- batch+1
            subset <- 1:min(batchSize, length(ind))
            subSeq <- seqData[ind[subset], "probe sequence"] # The next batch of sequences
            modMatrixMat <- setMatMatrix(subSeq) # Get the model matrix for the sequences
            message(sprintf("Fitting Values for batch %i of %i...", batch, nBatches))
            flush.console()
            fittedMat <- modMatrixMat%*%bgParam$coef # Get the fitted values for the current batch of probes
            # Check the dimensions of the fittedMat object!!!
            prbIds <- as.numeric(seqData[ind[subset], "Probe ID"]) # The probe IDs

            message("done\n")
            message("Updating CEL files...")
            flush.console()
            tempBins <- findInterval(fittedMat[,1], bgBins$boundaries) + 1 # The 1 needs to be added as the lower boundary is -Inf & not listed
            for (i in 1:nChips) {

                if (!lambdaDone) {
                    lambda <- bgBins$lambda[tempBins, i]
                    lambdaFile <- file.path(lambdaPath, paste(celNames[i], "CEL", sep="."))
                    updateCel(lambdaFile, indices=prbIds, intensities=exp(lambda)) # Must be on the exp scale such that 0<lambda<16log(2)
                }
                if (!deltaDone) {
                    delta <- bgBins$delta[tempBins, i]
                    deltaFile <- file.path(deltaPath, paste(celNames[i], "CEL", sep="."))
                    updateCel(deltaFile, indices=prbIds, intensities=exp(delta))
                }
            }
            message("done\n")

            fittedMat <- c() # Clear from the memory for efficiency
            prbIds <- c() # Clear from the memory for efficiency
            ind <- ind[-subset] # Remove the subset of probes from the ind vector

        }
    }
    flush.console()

    # If using GC bins
    if (method=="GC") {

        # Get the min & max GC counts from the bg probes
        minGc <- min(as.numeric(substr(rownames(bgBins$lambda), start=3, stop=nchar(rownames(bgBins$lambda)))))
        maxGc <- max(as.numeric(substr(rownames(bgBins$lambda), start=3, stop=nchar(rownames(bgBins$lambda)))))

        while (length(ind)>0) {

            batch <- batch+1
            subset <- 1:min(batchSize, length(ind))
            subSeq <- seqData[ind[subset], "probe sequence"]
            prbIds <- as.numeric(seqData[ind[subset], "Probe ID"]) # The probe IDs
            message(sprintf("Getting GC counts for batch %i of %i.\n", batch, nBatches))
            tempBins <- gcCounts(subSeq) # Get the counts for the current batch of sequences
            tempBins[which(tempBins<minGc)] <- minGc # The lowest GC count in the bg probes
            tempBins[which(tempBins>maxGc)] <- maxGc # The highest count in the bg probes
            tempBins <- paste("GC", tempBins, sep="")
            message("Updating CEL files...")
            if (!lambdaDone) {
                lambda <- bgBins$lambda[tempBins, ]
                for (i in 1:nChips) {
                    lambdaFile <- file.path(lambdaPath, paste(celNames[i], "CEL", sep="."))
                    updateCel(lambdaFile, indices=prbIds, intensities=exp(lambda[,i]))
                }
            }
            flush.console()

            if (!deltaDone) {
                delta <- bgBins$delta[tempBins, ]
                for (i in 1:nChips) {
                    deltaFile <- file.path(deltaPath, paste(celNames[i], "CEL", sep="."))
                    updateCel(deltaFile, indices=prbIds, intensities=exp(delta[,i]))
                }
            }
            message("done\n")
            flush.console()

            prbIds <- c() # Clear from the memory for efficiency
            ind <- ind[-subset] # Remove the subset of probes from the ind vector

        }
    }

    message(sprintf("Background Priors have been written to:\n %s \n %s \n", lambdaPath, deltaPath))

    # Form the celSets
    lambda <- AffymetrixCelSet$byPath(lambdaPath, checkChipType=FALSE)
    delta <- AffymetrixCelSet$byPath(deltaPath, checkChipType=FALSE)
    out <- list(lambda=lambda, delta=delta)
    class(out) <- "AffymetrixCelSetList" # Ensure that it has the correct class attribute

    # Return a list with the two celSets
    return(out)

}