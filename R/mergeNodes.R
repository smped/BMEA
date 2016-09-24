#' @title Merge the CEL files written on multiple nodes into a single set
#' 
#' @description Takes the node specific sets of CEL files containing BMEA contrast & 
#' model parameter data and merges them into a single set of CEL files.
#' 
#' @details During analysis using parallel nodes, 
#' each node will write to a batch of CEL files in its own separate directory. 
#' This function collects the units that have been run on each node and copies the saved results 
#' from the CEL files on each node to a main set of CEL files in the root directory.
#' 
#' CEL files with the contrast data (i.e. logFC & phiLogFC) will be written to the directories 
#' \code{bmeaData/contrastData/parentName,logFC/chipType} & \code{bmeaData/contrastData/parentName,phiLogFC/chipType} respectively. 
#' If no output directory is set via \code{outDir} the current working directory will be used.
#' 
#' CEL files with model parameters will be written to the directories \code{bmeaData/modelData/parentName,parameterName/chipType}.
#' 
#' @param celSet the parent celSet being analysed
#' @param nodePaths where to find all the node specific CEL files
#' @param outDir the directory to merge the node-level information to.
#' @param ... not used
#' @param paramToWrite the model parameters to collect from each node & write to the merged CEL files
#' @param nodeUnits a \code{\link{list}} with a vector of units for each node that correspond to those processed on each node. Overrides any information contained in the log files. 
#' The units will be extracted from the log files if this is not supplied
#' @param visible if \code{visible=TRUE} a \code{\link{list}} will be returned, 
#' otherwise the function will invisibly return \code{TRUE}. Defaults to FALSE
#' @param verbose determines the amount of information to display. Defaults to \code{TRUE}
#' 
#' @return  
#' If \code{visible=FALSE} will invisibly return \code{TRUE}, otherwise a \code{\link{list}} will be 
#' returned with the following elements:
#' \itemize{
#' \item{$written}{ \code{TRUE}}
#' \item{$parameters}{ contains the model parameters that have been written to the set of CEL files.}
#' \item{$nodeUnits}{ the units processed on each node}
#' \item{$nodePaths}{ the paths the node specific CEL files were sourced from}
#' \item{$logFiles}{ the contents of the log files from each node}
#' }
#' 
#' @seealso \code{\link{AffymetrixCelSet}}
#' 
#' @export
mergeNodes <- function(celSet, nodePaths, outDir=NULL, ..., paramToWrite=NULL, nodeUnits=NULL, visible=FALSE, verbose=TRUE) {

  
  # This function is designed to write the separate node outputs into a single set of CEL files
  
  # The required objects will be:
  # 1 - the original celSet
  # 2 - the nodePaths
  
  # The optional parameters will be:
  # 1 - A list of units from each node
  # 2 - A vector of parameters to merge
  # 3 - A vector of contrasts to merge (or a contrast matrix)
  # 4 - verbose
  
  
    ##########################
    # First check the celSet #
    ##########################
    if (length(grep("AffymetrixCelSet",class(celSet)))==0) stop("The supplied celSet is not  a celSet.\n")

    # Get the required info from the celSet
    parentName <- getParentName(celSet) # The basic name of the dataset
    nChips <- nbrOfArrays(celSet)
    # The cdf assigned to the celSet & it's related info
    cdf <- getCdf(celSet)
    chipType <- getChipType(cdf) # The chipType
    nUnits <- nbrOfUnits(cdf) # The number of units on the cdf
    cdfPath <- getPath(cdf) # The directory where the cdf is
    cdfType <- getName(cdf) # The original chipType without any further tags
    if (any(grepl(paste(chipType,"monocell.[Cc][Dd][Ff]",sep=","), 
                  list.files(cdfPath)))) { # Of course this should exist as it has been used to write the files on each node
        monoCdf <- AffymetrixCdfFile$byChipType(paste(chipType,",monocell",sep=""))
    }
    else {
        monoCdf <- createMonocellCdf(cdf)
    }

    #######################
    # Check the nodePaths #
    #######################
    # Check it is a list
    if (!is.list(nodePaths)) stop("The node paths must be specified as a list.\n")
    nNodes <- length(nodePaths)

    if (verbose) message("Checking directory structures..."); flush.console()
    if (length(which(!file.exists(file.path(nodePaths, "bmeaData"))))!=0) stop("The nodePaths are not setup as specified.\n")
    # Check the contrast directories for logFC & phiLogFC
    if (length(which(!file.exists(file.path(nodePaths, "bmeaData", "contrastData", paste(parentName,"logFC",sep=",")))))!=0) stop("The contrastData directories are not setup as required for logFC.\n")
    if (length(which(!file.exists(file.path(nodePaths, "bmeaData", "contrastData", paste(parentName,"phiLogFC",sep=",")))))!=0) stop("The contrastData directories are not setup as required for phiLogFC.\n")
    # Check the directories for the modelData parameters to write
    if (length(which(!file.exists(file.path(nodePaths, "bmeaData", "modelData"))))!=0) stop("The modelData directories are not setup as required.\n")
    
    #############################################
    # Check the parameters to write & the logFC #
    #############################################
    if (is.null(paramToWrite)) {
        # If not specified, find what is there in node1 & write all model parameters
        paramToWrite <- substr(dir(file.path(nodePaths[[1]], "bmeaData","modelData")),start=nchar(parentName)+2, stop=nchar(dir(file.path(nodePaths[[1]], "bmeaData","modelData"))))
        # This won't need to be checked as it's been extracted from what's there
    }
    else {
        # Check that the parameters requested have been written on the first node. If any can't be found, remove them from the list
        if (length(which(!file.exists(file.path(nodePaths[[1]], "bmeaData","modelData", paste(parentName, paramToWrite, sep=",")))))!=0) message("Some requested parameters cannot be found & will be ignored.\n")
        paramToWrite <- paramToWrite[which(file.exists(file.path(nodePaths[[1]], "bmeaData","modelData", paste(parentName, paramToWrite, sep=","))))]
    }
    nParam <- length(paramToWrite)
    if (nParam==0) message("NB: No parameter files will be written as no valid parameters have been requested.\n"); flush.console()
    modFiles <- vector("list",nParam)
    names(modFiles) <- paramToWrite
    for (i in 1:nParam) {
        modFiles[[i]] <- dir(file.path(nodePaths[[1]], "bmeaData","modelData", paste(parentName, paramToWrite, sep=","), cdfType)[[i]]) # Get the names of the files in the directory from node1 as these will be replicated across all nodes
    }

    # Compare the contrast directories & make sure they have the same structure
    # The basic structure of all the logFC &phiLogFC directories should be the same:
    lfcFiles <- dir(file.path(nodePaths[[1]], "bmeaData","contrastData",paste(parentName,"logFC",sep=","),cdfType))
    # Now compare all the directories from the other nodes
    for (i in 1:nNodes) {
        # Check the modelData directories
        for (j in 1:nParam) {
            if(length(which(!modFiles[[j]] %in% dir(file.path(nodePaths[[i]], "bmeaData","modelData",paste(parentName,paramToWrite[j],sep=","),cdfType))))!=0) stop("The modelData directories are not identical.\n")
        }
        # Check the logFC directories
        if (length(which(!lfcFiles %in% dir(file.path(nodePaths[[i]], "bmeaData","contrastData",paste(parentName,"logFC",sep=","),cdfType))))!=0) stop("The logFC directories are not identical.\n")
        # Check the phiLogFC folders
        if (length(which(!lfcFiles %in% dir(file.path(nodePaths[[i]], "bmeaData","contrastData",paste(parentName,"phiLogFC",sep=","),cdfType))))!=0) stop("The phiLogFC directories are not identical.\n")
    }
    if (verbose) message("...done!\n"); flush.console()
    # Now get the contrast names, although this may not be neccessary in reality
    contNames <- unique(matrix(unlist(strsplit(lfcFiles,",")),ncol=2,byrow=TRUE)[,1])

    ##########################################
    # Check the directories to be written to #
    ##########################################

    # First check that the directories exist & create them if required
    if (verbose) message("Checking write directories..."); flush.console()
    if (is.null(outDir))  outDir <- getwd()
    # Set the root directory for contrasts & modelData
    contRoot <- file.path(outDir, "bmeaData","contrastData")
    modelRoot <- file.path(outDir, "bmeaData","modelData")

    # Now set the entire set of paths & create any neccessary
    writePaths <- file.path(contRoot,paste(parentName,c("logFC","phiLogFC"),sep=","),cdfType)
    writePaths <- c(writePaths, file.path(modelRoot, paste(parentName, paramToWrite, sep=","), cdfType))
    for (i in 1:length(writePaths)) {
        if (!file.exists(writePaths[i])) dir.create(writePaths[i], recursive=TRUE) # Create the directory if it doesn't exist
    }
    if (verbose) message("done!\n"); flush.console()

    # If the files don't exist, then create them. It's probably best to check them one at a time in a loop
    if (verbose) message("Checking for pre-existing contrast Cel files...\n"); flush.console()
    for (i in 1:length(lfcFiles)) {

        # Create the file if it doesn't exist. All logFC files will use the monocell cdf
        curFile <- file.path(contRoot, paste(parentName,"logFC",sep=","), cdfType, lfcFiles[i])
        if (!file.exists(curFile)) {
            message(sprintf("\tCreating %s ...", curFile)); flush.console()
            tempHdr <- createCelHeader(lfcFiles[i], monoCdf)
            createCel(curFile, tempHdr, overwrite=TRUE)
            message(sprintf("done!\n", curFile))
        }
        else {
            chk <- unlist(readCelHeader(curFile)) %in% unlist(readCelHeader(file.path(nodePaths[[1]],"bmeaData","contrastData",paste(parentName,"logFC",sep=","),cdfType, lfcFiles[i])))
            # The header sections 1, 9 & 10 willnot match as they contain directory & timestamp information. Check the rest of the data quite stringently.
            # If the two files check out, then life is good we can move on.
            if (length(which(!chk[-c(1,9,10)]))!=0) {
                # Spit out a message & create a new file
                message(sprintf("The file %s will be overwritten as it did not match the source files.\n",curFile))
                message(sprintf("\tCreating %s ...", curFile)); flush.console()
                tempHdr <- createCelHeader(lfcFiles[i], monoCdf)
                createCel(curFile, tempHdr, overwrite=TRUE)
                message(sprintf("done!\n", curFile))
            }
        }

        # Repeat for phiLogFC
        curFile <- file.path(contRoot, paste(parentName,"phiLogFC",sep=","), cdfType, lfcFiles[i])
        if (!file.exists(curFile)) {
            message(sprintf("\tCreating %s ...", curFile)); flush.console()
            tempHdr <- createCelHeader(lfcFiles[i], monoCdf)
            createCel(curFile, tempHdr, overwrite=TRUE)
            message(sprintf("done!\n", curFile))
        }
        else {
            chk <- unlist(readCelHeader(curFile)) %in% unlist(readCelHeader(file.path(nodePaths[[1]],"bmeaData","contrastData",paste(parentName,"phiLogFC",sep=","),cdfType, lfcFiles[i])))
            # The header sections 1, 9 & 10 willnot match as they contain directory & timestamp information. Check the rest of the data quite stringently.
            # If the two files check out, then life is good we can move on.
            if (length(which(!chk[-c(1,9,10)]))!=0) {
                # Spit out a message & create a new file
                message(sprintf("The file %s will be overwritten as it did not match the source files.\n",curFile))
                message(sprintf("\tCreating %s ...", curFile)); flush.console()
                tempHdr <- createCelHeader(lfcFiles[i], monoCdf)
                createCel(curFile, tempHdr, overwrite=TRUE)
                message(sprintf("done!\n", curFile))
            }
        }

    }
    if (verbose) message("Checking for pre-existing contrast Cel files...done!\n"); flush.console()

    # Now check the modelData files
    if (verbose) message("Checking for pre-existing parameter Cel files...\n"); flush.console()
    for (j in 1:nParam) {
        for (i in 1:length(modFiles[[j]])) {
            curParam <- names(modFiles[j])
            curFile <- file.path(modelRoot, paste(parentName, curParam, sep=","), cdfType, modFiles[[j]][i])
            if (!file.exists(curFile)) {
                message(sprintf("\tCreating %s ...", curFile)); flush.console()
                if (tolower(curParam)=="s" || tolower(curParam)=="p") {
                    tempHdr <- createCelHeader(modFiles[[j]][i], cdf)
                    createCel(curFile, tempHdr, overwrite=TRUE)
                }
                else {
                    tempHdr <- createCelHeader(modFiles[[j]][i], monoCdf)
                    createCel(curFile, tempHdr, overwrite=TRUE)
                }
                message(sprintf("done!\n", curFile))

            }
            else {
                chk <- unlist(readCelHeader(curFile)) %in% unlist(readCelHeader(file.path(nodePaths[[1]],"bmeaData","modelData",paste(parentName, curParam,sep=","),cdfType, modFiles[[j]][i])))
                # The header sections 1, 9 & 10 willnot match as they contain directory & timestamp information. Check the rest of the data quite stringently.
                # If the two files check out, then life is good we can move on.
                if (length(which(!chk[-c(1,9,10)]))!=0) {
                    # Spit out a message & create a new file
                    message(sprintf("The file %s will be overwritten as it did not match the source files.\n",curFile))
                    message(sprintf("\tCreating %s ...", curFile)); flush.console()
                    if (tolower(curParam)=="s" || tolower(curParam)=="p") {
                        tempHdr <- createCelHeader(modFiles[[j]][i], cdf)
                        createCel(curFile, tempHdr, overwrite=TRUE)
                    }
                    else{
                        tempHdr <- createCelHeader(modFiles[[j]][i], monoCdf)
                        createCel(curFile, tempHdr, overwrite=TRUE)
                    }
                    message(sprintf("done!\n", curFile))
                }
            }

        }
    }
    if (verbose) message("Checking for pre-existing parameter Cel files...done!\n"); flush.console()

    #######################################################################
    # Get the nodeUnits from the logFiles, unless they have been supplied #
    #######################################################################
    if (is.null(nodeUnits)) { # Create the object if not supplied
        if (verbose) message("Reading nodeUnits from log files...")
        nodeUnits <- vector("list", nNodes)
        names(nodeUnits) <- unlist(nodePaths)
        logFiles <- nodeUnits # Create a copy for saving the complete set of information
        for (i in 1:nNodes) {
            if (!file.exists(file.path(nodePaths[[i]], "bmeaData", paste(parentName, "log","txt",sep=".")))) stop(sprintf("Missing log file: %s\n", file.path(nodePaths[[i]], "bmeaData", paste(parentName, "log","txt",sep="."))))
            # The units will be in the first column
            logFiles[[i]] <- read.table(file.path(nodePaths[[i]], "bmeaData", paste(parentName, "log","txt",sep=".")), header=TRUE, stringsAsFactors=FALSE)
            nodeUnits[[i]] <- logFiles[[i]][,1]
        }
        if (verbose) message("done!\n")
    }
    else { # If it is supplied, check it matches the nodes
        if (length(nodeUnits)!=nNodes) stop("The nodeUnits do not match the number of nodes.\n")
        logFiles <- vector("list", nNodes)
        names(logFiles) <- unlist(nodePaths)
        for (i in 1:nNodes) {
            if (!file.exists(file.path(nodePaths[[i]], "bmeaData", paste(parentName, "log","txt",sep=".")))) stop(sprintf("Missing log file: %s\n", file.path(nodePaths[[i]], "bmeaData", paste(parentName, "log","txt",sep="."))))
            # The units will be in the first column
            logFiles[[i]] <- read.table(file.path(nodePaths[[i]], "bmeaData", paste(parentName, "log","txt",sep=".")), header=TRUE, stringsAsFactors=FALSE)
        }
    }

    # Check that no units are duplicated
    if (length(which(duplicated(unlist(nodeUnits))))!=0) stop("Some units appear to be present on more than one node.\n")

    ##########################################################################
    # Now extract the data from each node & write it to the master CEL files #
    ##########################################################################
    for (i in 1:nNodes) {

        if (verbose) message(sprintf("Extracting data from %s:\n",nodePaths[[i]]))
        # NB: For the contrastData files, the monocell cdf will always be used:
        curUgcMap <- getUnitGroupCellMap(monoCdf, units=nodeUnits[[i]])

        # Starting with logFC copy them one celFile at a time. Copy all cells, even though for logFC, we really only need the first one
        if (verbose) message("\tUpdating logFC..."); flush.console()
        for (j in 1:length(lfcFiles)){
            outFile <- file.path(contRoot,paste(parentName, "logFC",sep=","), cdfType, lfcFiles[j])
            readFile <- file.path(nodePaths[[i]], "bmeaData", "contrastData", paste(parentName, "logFC", sep=","), cdfType, lfcFiles[j])
            tempData <- readCelIntensities(readFile, indices=curUgcMap$cell)
            updateCel(outFile, curUgcMap$cell, as.vector(tempData))
        }
        if (verbose) message("done!\n")

        # Now do phiLogFC
        if (verbose) message("\tUpdating phiLogFC..."); flush.console()
        for (j in 1:length(lfcFiles)){
            outFile <- file.path(contRoot,paste(parentName, "phiLogFC",sep=","), cdfType, lfcFiles[j])
            readFile <- file.path(nodePaths[[i]], "bmeaData", "contrastData", paste(parentName, "phiLogFC", sep=","), cdfType, lfcFiles[j])
            tempData <- readCelIntensities(readFile, indices=curUgcMap$cell)
            updateCel(outFile, curUgcMap$cell, as.vector(tempData))
        }
        if (verbose) message("done!\n"); flush.console()

        # And now do the parameters.
        # Most will use the moncell Cdf (i.e. c, mu, sigmaMu, sigmaP, sigmaS, phi), however some (i.e. S & p) will use the full cdf
        # Update in this order to ensure that the units don't have to be reeextracted unless absolutely necessary (i.e. if S & p are being written)
        # Move S & p to the end of the list of modelFiles
        if ("s" %in% tolower(names(modFiles))) {
            sInd <- match("s", tolower(names(modFiles)))
            modFiles <- modFiles[c(seq(1:length(modFiles))[-sInd],sInd)]
        }
        if ("p" %in% tolower(names(modFiles))) {
            pInd <- match("p", tolower(names(modFiles)))
            modFiles <- modFiles[c(seq(1:length(modFiles))[-pInd],pInd)]
        }

        # Now update the model parameter files:
        for (j in 1:nParam) {
            curParam <- names(modFiles)[j]
            if (verbose) message(sprintf("\tUpdating %s...",curParam)); flush.console()
            # If we need to use the main cdf
            if (tolower(curParam)=="s" || tolower(curParam)=="p") curUgcMap <- getUnitGroupCellMap(cdf, units=nodeUnits[[i]])
            for (k in 1:length(modFiles[[j]])) {
                outFile <- file.path(modelRoot,paste(parentName, curParam ,sep=","), cdfType, modFiles[[j]][k])
                readFile <- file.path(nodePaths[[i]], "bmeaData", "modelData", paste(parentName, curParam ,sep=","), cdfType, modFiles[[j]][k])
                tempData <- readCelIntensities(readFile, indices=curUgcMap$cell)
                updateCel(outFile, curUgcMap$cell, as.vector(tempData))
            }
            if (verbose) message("done!\n")
        }
        flush.console()

    }

    # Tidy up the logFile(s)
    logFile <- c()
    for (i in 1:nNodes) {
        logFile <- rbind(logFile, cbind(i,logFiles[[i]]))
    }
    colnames(logFile) <- c("node",colnames(logFiles[[1]]))

    # Write the merged logFile to disk
    write.table(logFile, file.path(outDir, "bmeaData",paste(parentName,"log","txt",sep=".")), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

    # Return the output with some useful information if visible=TRUE, otherwise just invisibly return 'TRUE'
    if (visible) {
        return(list(written=TRUE, parameters=paramToWrite, nodeUnits=nodeUnits, nodePaths=nodePaths, logFile=logFile))
    }
    else {
        invisible(TRUE)
    }

}
