#' @title Fit the BMEA model for a large number of units
#' 
#' @description Fit a large number of units in batches using the BMEA model
#' 
#' @details 
#' This is the high-level function which fits the specified units in batches, 
#' as determined by \code{batchSize}. 
#' Each batch of units (i.e. genes) is fitted, then the parameters written to disk. 
#' Generally writing to disk can be a bottleneck when it comes to processing time,
#' so a wise choice for \code{batchSize} can impact on the speed of the process. 
#' 
#' @param celSet the celSet which is being processed
#' @param bgCelSet a list with components \code{lambda} & \code{delta}. 
#' Each component must be a celSet containing the information for the background signal priors
#' @param units the units to be processed
#' @param batchSize the number of units to process before writing to disk. Defaults to 16 (2^4)
#' @param conditions a vector of factors
#' @param contMatrix the contrasts
#' @param ... used for passing arguments to \code{fitBmeaBatch}, \code{writeBmeaBatch} & \code{runMCMC.BMEA}
#' @param verbose controls the level of detail to be displayed whilst the process is running
#' 
#' @return \code{NULL}
#' 
#' @seealso 
#'  \code{\link{fitBmeaBatch}}, \code{\link{writeBmeaBatch}}, \code{\link{runMCMC.BMEA}}
#' 
#' @import aroma.affymetrix
#' 
#' @export
fitBmea <- function(celSet, bgCelSet, units=NULL, batchSize=2^4, conditions, contMatrix, ..., verbose=0) {

    # Check the celSet
    if (length(grep("AffymetrixCelSet",class(celSet)))==0) return(cat("The supplied celSet is not actually a celSet.\n"))
    cdf <- getCdf(celSet)
    nChips <- nbrOfArrays(celSet)

    # Check the bgCelSet
    if (length(which(c("lambda","delta") %in% attributes(bgCelSet)$names))!=2) return(cat("The supplied bgCelSet must have components $delta & $lambda.\n"))
    if (length(grep("AffymetrixCelSet",class(bgCelSet$lambda)))==0) return(cat("The supplied bgCelSet$lambda is not a celSet.\n"))
    if (length(grep("AffymetrixCelSet",class(bgCelSet$delta)))==0) return(cat("The supplied bgCelSet$delta is not a celSet.\n"))
    if (getChipType(bgCelSet$lambda)!=getChipType(bgCelSet$delta)) return(cat("The chipTypes in the bgCelSet do not match.\n"))
    if (nbrOfArrays(bgCelSet$lambda)!=nbrOfArrays(bgCelSet$delta)) return(cat("The nbrOfArrays in the bgCelSets do not match.\n"))

    # Compare the two celSets
    if (getChipType(celSet) != getChipType(bgCelSet$lambda)) return(cat("The bgCelSet is a different chipType to the main celSet.\n"))
    if (nChips != nbrOfArrays(bgCelSet$lambda)) return(cat("The nbrOfArrys in the bgCelSet doesn't match that of the main celSet.\n"))

    # Check the conditions match the celSets
    if (!is.factor(conditions)) return(cat("The conditions must be specified as a vector of factors.\n"))
    if (length(conditions) != nChips) return(cat("The number of conditions doesn't match the number of arrays.\n"))
    # Check that the contMatrix matches the conditions:
    if (nrow(contMatrix)!=length(levels(conditions))) return(cat("The contrast matrix doesn't match the conditions.\n"))
    nCont <- ncol(contMatrix)

    # Check the units are on the cdf
    if (is.null(units)) {
        units <- seq(1,nbrOfUnits(cdf),1) # If no units are supplied, do the entire set of units
        cat("Units were not specified. Fitting all units!\n")
    }
    if (!is.numeric(units)) return(cat("The units must be specified numerically.\n"))
    if (max(units) > nbrOfUnits(cdf)) return(cat(sprintf("Invalid unit(s). The maximum unit number possible is %i.\n",nbrOfUnits(cdf))))

    # Once everything checks out, start with the process
    nBatches <- as.integer(ceiling(length(units)/batchSize))
    batch <- as.integer(0) # Set the first batch
    ind <- units # Pass the units to another variable which is modified
    batchLog <- list(units=NULL,parameters=NULL,time=NULL)

    while (length(ind)>0) {

        batch <- batch+1 # Only used for displaying progress if verbose>=1
        subset <- 1:min(batchSize, length(ind)) # Set the units to be the next batch
        curUnits <- ind[subset]

        # Fit the batch & write to disk
        if (verbose>0) cat(sprintf("Fitting batch %i of %i...\n", batch, nBatches))
        batchData <- fitBmeaBatch(celSet=celSet, bgCelSet=bgCelSet, units=curUnits, conditions=conditions, contMatrix=contMatrix, ...)
        if (verbose>0) cat("Updating parameter CEL files...\n")
        # The writeBmeaBatch function should now also have a log file associated with it
        lastBatch <- writeBmeaBatch(batchData, ..., verbose=verbose)
        ind <- ind[-subset]

        batchLog$units <- c(batchLog$units,lastBatch$units)
        if (is.null(batchLog$parameters)) batchLog$parameters <- c(batchLog$parameters,lastBatch$parameters)
        batchLog$time <- c(batchLog$time,lastBatch$time)

    }

    return(batchLog)

}