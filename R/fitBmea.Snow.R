#' @title Fit the BMEA process for a number of genes utilising the snow package
#' 
#' @description Enables genes to be fit in parallel, utilising the architecture of the snow package
#' 
#' @param celSet the celSet which is being processed
#' @param bgCelSet a list with components \code{lambda} & \code{delta}. 
#' Each component must be a celSet containing the information for the background signal priors
#' @param cl A cluster as formed by \code{makeCluster} or similar from the package \code{snow}
#' @param units the units to be processed
#' @param batchSize the number of units to process before writing to disk. Defaults to 16 (2^4)
#' @param conditions a vector of factors
#' @param contMatrix the contrasts
#' @param ... used for passing arguments to \code{fitBmeaBatch}, \code{writeBmeaBatch} & \code{runMCMC.BMEA}
#' @param nodePaths Where each node in the cluster will write it's output
#' 
#' @details This process will split the supplied units into vectors of equal length, 
#' with one for each node in the supplied cluster. 
#' The BMEA process will be fit in series for the vector of units on each parallel node. 
#' Depending on the computing architecture being used & the size of the dataset, 
#' this process can still take a number of hours, to a number of days.
#' 
#' Due to the architecture of the \code{snow} package, all arguments will be first sent to
#' .GlobalEnv with the prefix 'bmeaTemp.'. 
#' Any objects in .GlobalEnv with this prefix will be deleted as part of the process, 
#' so users should avoid using this string as part of the name of any objects.
#'  
#' On each node, the genes will be fit in series & the results saved to disk periodically as a batch,
#' with the results saved as CEL files. 
#' The number of genes fit in each batch is determined by the \code{batchSize} argument. 
#' Each time the results for a batch of units are written to disk, 
#' a log file with the name \code{parentName.log.txt} will be updated to ensure minimum data is lost in the case of any interuptions. 
#' These log files are recorded in a node-specific manner in the directory "nodePath/bmeaData". 
#' The results for each parameter in the \code{paramToSave} argument will be written to CEL files 
#' on each node in the directories \code{nodePath/bmeaData/modelData}, 
#' whilst the contrasts specified in the contrast matrix will be written to the directories \code{nodePath/bmeaData/contrastData}.
#' 
#' If the \code{mcmcParam} argument is not specified as part of the '...' the defaults of \code{nChains=3},
#' \code{nIter=16000}, \code{nBurnin=8000} & \code{nThin=8} will be used. 
#' The default values for \code{zGene} (4.265) & \code{zExon} (1.645) will also be used if not specified.
#' 
#' If the \code{paramToSave} argument is not specified, the default of the chip effects ("c"), 
#' condition-specific effects ("mu") and the exon-level terms ("phi") will be saved. 
#' If the \code{paramToWrite} argument is not specified, it will default to the same as the \code{paramToSave} argument.
#' 
#' After the process has been run, the saved data for the sets of parallel units can be merged via the function \code{\link{mergeNodes}}.
#' 
#' @return 
#' Returns a list object with the elements:
#' \itemize{
#' \item{$units}{ a list with each sub-element containing the units fitted on each node. 
#' Each sub element is named with the provided \code{nodePath}}
#' \item{$times}{ a list with the times each batch was written to disk for each node.}
#' \item{$parameters}{ the model parameters that have been saved on each node}
#' }
#' 
#' @seealso 
#'  \code{\link{mergeNodes}}, \code{\link{AffymetrixCelSet}}, 
#'  \code{\link{AffymetrixCelSetList}}, \code{\link{makeCluster}}
#' 
#' 
#' @import snow
#' @import aroma.affymetrix
#' 
#' @export
fitBmea.Snow <- function(celSet, bgCelSet, cl, units=NULL, batchSize=NULL, conditions, contMatrix, ..., nodePaths=NULL) {

    ###########################
    # First check the cluster #
    ###########################
    if(length(grep("cluster",class(cl)))==0) stop("The supplied cluster must actually be a cluster.\n")
    nNodes <- length(cl)

    # Check &/or set the nodePaths
    if (is.null(nodePaths)) {
        nodePaths <- paste("node",1:nNodes,sep="")
    }
    else {
        if (length(nodePaths)<nNodes) { # If some names haven't been specified, create them
            for (i in (length(nodePaths)+1):nNodes) nodePaths[i] <- paste("node",i,sep="")
        }
    }
    nodePaths <- nodePaths[1:nNodes] # Get rid of any extra nodePaths

    ###############################
    # Check the celSet & bgCelSet #
    ###############################
    if (length(grep("AffymetrixCelSet",class(celSet)))==0) stop("Error: Supplied celSet is not a celSet.\n")
    nChips <- nbrOfArrays(celSet)
    chipType <- getChipType(celSet)
    cdf <- getCdf(celSet)
    if (class(bgCelSet)!="AffymetrixCelSetList") stop("Error: Supplied bgClSet is not an AffymetrixCelSetList.\n")
    if (length(match(names(bgCelSet), c("lambda","delta")))!=2) stop("Error: bgCelSet must only have elements lambda & delta.\n")
    if (length(grep("AffymetrixCelSet",class(bgCelSet[[1]])))==0 || length(grep("AffymetrixCelSet",class(bgCelSet[[2]])))==0) stop("Error: Both elements of the bgCelSet object must be AffymetrixCelSets.\n")
    if (getChecksum(getCdf(bgCelSet[[1]])) != getChecksum(getCdf(bgCelSet[[2]]))) stop("Error: Both elements of the bgCelSet object must have identical cdf files.\n")
    if (nbrOfArrays(bgCelSet[[1]])!=nChips || nbrOfArrays(bgCelSet[[2]])!=nChips) stop("Error: The bgCelSet has a different number of arrays to the parent celSet.\n")
    if (getChipType(bgCelSet[[1]])!=chipType || getChipType(bgCelSet[[2]])!=chipType) stop("Error: The chipType must be the same for the bgCelSet & the parent celSet.\n")

    ###########################
    # Check & setup the units #
    ###########################
    maxUnits <- nbrOfUnits(cdf)
    if (is.null(units)) {
        units <- seq(1, maxUnits, by=1)
        nUnits <- length(units)
    }
    else {
        if (!is.numeric(units)) stop("Error: The units must be supplied numerically.\n")
        if (max(units)>maxUnits) stop(sprintf("Error: Invalid units. The maximum available unit is %i\n.", maxUnits))
        nUnits <- length(units)
        if (nUnits<nNodes) nNodes <- nUnits # If there are less units than nodes, just use the minimum required
    }
    nodeUnits <- vector("list", nNodes)
    for (i in 1:nNodes){ # This will place every unit into each node in intervals of nNodes
        tempUnits <- seq(i,nUnits, by=nNodes)
        nodeUnits[[i]] <- units[tempUnits]
    }

    #####################################
    # Check the conditions & contMatrix #
    #####################################
    if (length(conditions)!=nChips) stop("Error: The conditions vector does not match the celSet.\n")
    if (!is.factor(conditions)) stop("Error: The conditions vector must be supplied as factors.\n")
    if (length(levels(conditions))!=length(rownames(contMatrix))) stop("Error: The conditions vector & contrast matrix are not compatible.\n")
    if (length(which(is.na(match(levels(conditions),rownames(contMatrix)))))!=0) stop("Error: The conditions do not match the rownames of the contrast matrix.\n")

    ##############################
    # Check the dot-dot-dot args #
    ##############################
    dotArgs <- list(...)
    if (length(dotArgs)==0) {
        # If none have been supplied
        dotArgs <- c()
    }
    # Check & create the possible named arguments here.

    # paramToSave
    if(is.null(dotArgs$paramToSave)) paramToSave <- c("c","mu","phi")
    else paramToSave <- dotArgs$paramToSave
    # Now ensure it can only be the model parameters:
    allParam <- c("S","sigma_S","c","mu","sigma_mu","p","sigma_p","phi")
    ind <- which(!is.na(match(allParam, paramToSave)))
    paramToSave <- allParam[ind]

    # keepSims
    if(is.null(dotArgs$keepSims)) keepSims <- FALSE
    else keepSims <- dotArgs$keepSims
    if (!is.logical(keepSims)) stop("Error: keepSims must be either TRUE/FALSE.\n")

    # zGene & zExon
    if(is.null(dotArgs$zGene)) zGene <- 4.265
    else zGene <- dotArgs$zGene
    if(is.null(dotArgs$zExon)) zExon <- 1.645
    else zExon <- dotArgs$zExon

    # paramToWrite should default to the same as paramToSave
    if(is.null(dotArgs$paramToWrite)) paramToWrite <- paramToSave
    else paramToWrite <- dotArgs$paramToWrite
    # Make sure it doesn't try to write any parameters not save
    paramToWrite <- paramToWrite[which(!is.na(match(paramToWrite, paramToSave)))]

    # mcmcParam
    if (is.null(dotArgs$mcmcParam)) {
        # Set to defaults if nor included
        mcmcParam <- list(nChains=3, nIter=16000, nBurnin=8000, nThin=8)
    }
    else {
        mcmcParam <- dotArgs$mcmcParam
        if (length(which(is.na(match(c("nChains","nIter","nBurnin","nThin"), names(mcmcParam)))))!=0) stop("Error: Invalid mcmcParam object.\n")
    }

    ####################
    # Set up the nodes #
    ####################

    # The key problem here is that objects need to be exported from the Global environment
    # A solution would be to send the key objects to the .GlobalEnv with a prefix & then remove them once they
    # have been exported. They could then be renamed within the cluster as well

    if (is.null(batchSize)) batchSize <- 2^4

    # Use the prefix "bmeaTemp" to store the objects in the Global Environment
    bmeaTemp.celSet <<- celSet
    bmeaTemp.bgCelSet <<- bgCelSet
    bmeaTemp.conditions <<- conditions
    bmeaTemp.contMatrix <<- contMatrix
    bmeaTemp.batchSize <<- batchSize
    bmeaTemp.nodeUnits <<- nodeUnits
    bmeaTemp.nodePaths <<- nodePaths
    bmeaTemp.mcmcParam <<- mcmcParam
    bmeaTemp.paramToSave <<- paramToSave
    bmeaTemp.keepSims <<- keepSims
    bmeaTemp.zGene <<- zGene
    bmeaTemp.zExon <<- zExon
    bmeaTemp.paramToWrite <<- paramToWrite

    # Export the objects to workspaces in each node
    snow::clusterEvalQ(cl, library(BMEA))
    snow::clusterExport(cl, list("bmeaTemp.celSet", "bmeaTemp.bgCelSet", "bmeaTemp.conditions", "bmeaTemp.contMatrix", "bmeaTemp.batchSize", "bmeaTemp.nodeUnits", "bmeaTemp.nodePaths", "bmeaTemp.mcmcParam", "bmeaTemp.paramToSave", "bmeaTemp.keepSims", "bmeaTemp.zGene", "bmeaTemp.zExon", "bmeaTemp.paramToWrite"))

    # Rename the objects on each node
    snow::clusterEvalQ(cl, celSet <- bmeaTemp.celSet)
    snow::clusterEvalQ(cl, bgCelSet <- bmeaTemp.bgCelSet)
    snow::clusterEvalQ(cl, conditions <- bmeaTemp.conditions)
    snow::clusterEvalQ(cl, contMatrix <- bmeaTemp.contMatrix)
    snow::clusterEvalQ(cl, batchSize <- bmeaTemp.batchSize)
    snow::clusterEvalQ(cl, mcmcParam <- bmeaTemp.mcmcParam)
    snow::clusterEvalQ(cl, paramToSave <- bmeaTemp.paramToSave)
    snow::clusterEvalQ(cl, keepSims <- bmeaTemp.keepSims)
    snow::clusterEvalQ(cl, zGene <- bmeaTemp.zGene)
    snow::clusterEvalQ(cl, zExon <- bmeaTemp.zExon)
    snow::clusterEvalQ(cl, paramToWrite <- bmeaTemp.paramToWrite)
    # Do the node specific objects:
    for (i in 1:nNodes){
        i <<- i
        snow::clusterExport(cl, "i")
        snow::clusterEvalQ(cl[i], units <- bmeaTemp.nodeUnits[[i]]) # Send the units
        snow::clusterEvalQ(cl[i], nodePath <- bmeaTemp.nodePaths[[i]]) # Send the nodePaths
    }

    # Remove the superfluous objects from each node
    snow::clusterEvalQ(cl[1:nNodes], rm(list=ls(pattern="bmeaTemp")))
    # Remove the tempBmea objects from the global environment
    rm(list=ls(pattern="bmeaTemp", envir=.GlobalEnv), envir = .GlobalEnv)

    # Split every cluster level variable using clusterSplit & then run clusterApply.
    out <- snow::clusterEvalQ(cl, fitBmea(celSet=celSet, bgCelSet=bgCelSet, units=units, batchSize=batchSize, conditions=conditions, contMatrix=contMatrix, mcmcParam=mcmcParam, paramToSave=paramToSave, keepSims=keepSims, zGene=zGene, zExon=zExon, paramToWrite=paramToWrite, nodePath=nodePath))

    # Format the output:
    units <- vector("list", nNodes)
    names(units) <- nodePaths
    times <- units
    for (i in 1:nNodes) {
        units[[i]] <- out[[i]]$units
        times[[i]] <- out[[i]]$time
    }
    parameters <- out[[1]]$parameters
    
    list(units=units, 
         times=times, 
         parameters=parameters)
}