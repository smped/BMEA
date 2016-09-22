#' @title Write the parameters to disk for a batch of genes
#' 
#' @description Writes the quantiles & summary statistics to disk for a batch of genes
#' 
#' @details 
#' This function will write the requested model parameters to disk, 
#' along with any contrasts simulated, for a batch of genes. 
#' The \code{batchData} input must be the output from \code{fitBmeaBatch}.
#' 
#' If \code{nodePath=NULL}, model parameters will be written to the subdirectory \code{bmeaData/modelData/parentName,parameterName/chipName},
#' otherwise the output path will be written to \code{nodePath/bmeaData/modelData/parentName,parameterName/chipName}.
#' 
#' If contrasts have been supplied, the contrast data will be written to the subdirectories 
#' \code{bmeaData/contrastData/parentName,logFC/chipName} & \code{bmeaData/contrastData/parentName,phiLogFC/chipName}. 
#' If a \code{nodePath} has been supplied, this path will similarly appear at the beginning of the specified paths.
#' 
#' If \code{writeLog=TRUE}, a log file will be written to the directory \code{bmeaData} with the name \code{parentName.log.txt}. 
#' The units, parameters and time & date informastion are stored in this file.
#' 
#' @param batchData contains the results for a batch of genes after analysis with BMEA. 
#' Must be the output from \code{\link{fitBmeaBatch}}}
#' @param paramToWrite the model parameters to write to disk. 
#' Any non-model parameters will be ignored
#' @param ... required for passing arguments from higher level functions
#' @param nodePath the path to be used if being analysed on a specific node.
#' @param writeLog specifies whether a log file will be written for each batch of genes
#' @param verbose determines how much information will be displayed
#' 
#' @return 
#' Returns a list with components:
#' \itemize{
#' \item{$units}{ the units that have just been written to disk}
#' \item{$parameters}{ the parameters for shich the units have been written}
#' \item{$time}{ the system time that the process was completed at, in POSIXct format}
#' }
#' 
#' @seealso \code{\link{fitBmeaBatch}}
#' 
#' @export
writeBmeaBatch <- function(batchData, paramToWrite=c("c","mu","phi"), ..., nodePath=NULL, writeLog=TRUE, verbose=0) {

    # The batchData must of class BMEA.Batch
    if (!isBMEA.Batch(batchData)) return(cat("The supplied data is not a BMEA.Batch object.\n"))
    if (is.null(batchData$units$units)) return(cat("No units have been recorded in the supplied BMEA.Batch object.\n"))
    units <- as.integer(batchData$units$units)
    conditions <- batchData$conditions

    if (verbose>0) cat(sprintf("Updating units %i through %i...\n",units[1], units[length(units)]))

    # Check the paramToWrite matches the paramToSave slot from the batch & get rid of any that don't match
    batchData$paramToSave <- tolower(batchData$paramToSave) # And change the other to lower case as well
    nullParam <- which(!tolower(paramToWrite) %in% batchData$paramToSave) # Find any requested parameters which have not been saved
    if (length(nullParam)!=0) {
        if (length(nullParam)==1) cat(sprintf("The parameter %s cannot be written to disk.\n",paramToWrite[nullParam]))
        if (length(nullParam)>1) {
            for (i in 1:length(nullParam)) {
                cat(sprintf("The parameter %s cannot be written to disk.\n",paramToWrite[nullParam[i]]))
            }
        }
        paramToWrite <- paramToWrite[-nullParam]
    }
    paramToWrite <- tolower(paramToWrite)

    # Check for the cdf & monocell cdf
    celSet <- batchData$celSet
    if (class(celSet)[1]!="AffymetrixCelSet") return(cat("The supplied celSet is not actually a celSet.\n"))
    cdf <- getCdf(celSet)
    chipName <- getName(cdf) # NB: This will be the root name, e.g. "HuEx-1_0-st-v2" of the chipType
    chipType <- getChipType(cdf) # This will be based on the cdf, e.g. "HuEx-1_0-st-v2,U-Ensembl49,G-Affy"
    cdfPath <- getPath(cdf) # This is where all the annotation files will be
    celPath <- getPath(celSet)
    celNames <- getNames(celSet)
    parentName <- getParentName(celSet)
    nChips <- nbrOfArrays(celSet)

    # Setup &/or check the monoCdf
    if (!file.exists(file.path(cdfPath,paste(chipType,",monocell.CDF",sep="")))) {
        # Create the monocell cdf
        monoCdf <- createMonoCell(cdf)
    }
    else { # Otherwise check the monocell cdf
        monoCdf <- AffymetrixCdfFile$byChipType(paste(chipType,",monocell",sep=""))
        if (!isMonocellCdf(monoCdf)) return(cat("The existing monocall CDF is not a monocell CDF.\n"))
    }

    # Set up where to store the model parameters:
    if (is.null(nodePath)) {
        modelPath <- file.path("bmeaData", "modelData")
    }
    else {
        modelPath <- file.path(nodePath, "bmeaData", "modelData")
    }
    if (!file.exists(modelPath)) dir.create(modelPath, recursive=TRUE)

    # Set up where to store the contrasts:
    if (is.null(nodePath)) {
        contrastPath <- file.path("bmeaData", "contrastData")
    }
    else {
        contrastPath <- file.path(nodePath, "bmeaData", "contrastData")
    }
    if (!file.exists(contrastPath)) dir.create(contrastPath, recursive=TRUE)

    # Create the necessary CEL files:
    if (is.null(paramToWrite)) {
        cat("Some parameters must be written to disk. Writing defaults of c, mu & phi.\n")
        paramToWrite <- c("c","mu","phi")
    }

    # Each parameter needs to be checked. Use a set of boolean variables:
    writeS <- ifelse(length(which(paramToWrite=="s"))!=0, TRUE, FALSE)
    writeSigmaS <- ifelse(length(which(paramToWrite=="sigma_s"))!=0, TRUE, FALSE)
    writeC <- ifelse(length(which(paramToWrite=="c"))!=0, TRUE, FALSE)
    writeMu <- ifelse(length(which(paramToWrite=="mu"))!=0, TRUE, FALSE)
    writeSigmaMu <- ifelse(length(which(paramToWrite=="sigma_mu"))!=0, TRUE, FALSE)
    writeP <- ifelse(length(which(paramToWrite=="p"))!=0, TRUE, FALSE)
    writeSigmaP <- ifelse(length(which(paramToWrite=="sigma_p"))!=0, TRUE, FALSE)
    writePhi <- ifelse(length(which(paramToWrite=="phi"))!=0, TRUE, FALSE)

    # Get the ugcMaps for the standard CEL files & the monocell CEL files:
    ugcMap <- getUnitGroupCellMap(cdf, units=units, retNames=TRUE)
    monoUgcMap <- getUnitGroupCellMap(monoCdf, units=units, retNames=TRUE)

    ####################
    # MODEL PARAMETERS #
    ####################

    # Name the statistics which are saved for each model parameter
    savedStats <- c("mean","sd","2.5%","25%","50%","75%","97.5%","rHat")
    if (verbose>0) cat("Updating model parameter files...\n")

    # S_hijk (Signal)
    # There will be an entry for every probe in every unit, so any .CEL files will need to be standard .CEL files
    # The data is already on the linear scale, so won't need to be transformed
    # There will need to be nbrOfArrays*8 CEL files (Mean, StDev, 2.5%, 25%, 50%, 75%, 97.5%, rHat)
    # Store in /bmeaData/modelData/parentName/chipType
    # Name each CEL file celName,S,mean.CEL, celName,S,sd.CEL, celName,S,025.CEL ... celName,S,975.CEL, celName,S,rHat.CEL
    # Generally not useful for posterior inference, so don't write!!!
    # S is also VERY cumbersome to save is it requires the larger (standard) CEL files & every unit must be written!
    if (writeS) {

        # Check & create the directory
        tags <- paste(parentName,"S",sep=",")
        sPath <- file.path(modelPath,tags,chipName)
        if (!file.exists(sPath)) dir.create(sPath, recursive=TRUE)

        # Check for the files, as the directory has been checked above:
        sNames <- paste(rep(celNames,each=8), savedStats, sep=",")
        sNames <- file.path(sPath, paste(sNames,"CEL",sep="."))
        # Reformat the sNames into a matrix (nChips * 8)
        sNames <- matrix(sNames, nrow=nChips, ncol=8, byrow=TRUE)
        if (length(which(!file.exists(sNames)))!=0) { # If one or more of the files don't exist
            for (i in 1:nChips) { # Step through all the chips
                for (j in 1:8) { # There are 8 statistics to keep
                    if (!file.exists(sNames[i,j])) {
                        hdr <- createCelHeader(filename=sNames[i,j], cdf=cdf, tags=tags)
                        createCel(sNames[i,j],hdr)
                    }
                }
            }
        }

        # Get the S values from the BMEA.Batch object as an array
        # Each original chip can be accessed by the 2nd dimension, i.e. cel 1 = S[,1,]
        # Each of the 8 stats can be accessed by the 3rd dimesion, i.e. rHat = S[,,8]
        S <- c()
        nS <- table(ugcMap$unit) # The number of units
        for (i in 1:length(units)) {
            rows <- 1:(nS[i]*nChips) # S will always be the first rows if it has been saved!
            S <- abind(S, array(batchData$summaries[[i]][rows,1:8],dim=c(nS[i],nChips,8)),along=1)
        }
        dimnames(S) <- list(ugcMap$cell,celNames,savedStats)

        for (i in 1:nChips) {
            for (j in 1:8){
                # The S parameter doesn't need to be exponentiated
                updateCel(sNames[i,j], indices=ugcMap$cell, intensities=S[,i,j])
            }
        }
    }

    # sigma_S (Signal Sd)
    # There will only be one value per unit, so CEL files can be monocell CEL files
    # The data can be exponentiated, although it won't be a huge issue.
    # There will be 8 moncell CEL files (Mean, StDev, 2.5%, 25%, 50%, 75%, 97.5%, rHat)
    # Store in /bmeaData/modelData/parentName,tags,sigma_S/chipType
    # Name the CEL files sigma_S,Mean.CEL
    # Generally not useful for inference, so don't write
    if (writeSigmaS) {

        # Check the conditions vector:
        if (length(conditions)!=nChips) return(cat("The conditions vector doesn't match the celSet.\n"))

        # Check & create the directory
        tags <- paste(parentName,"sigma_S",sep=",")
        sigma_SPath <- file.path(modelPath,tags,chipName)
        if (!file.exists(sigma_SPath)) dir.create(sigma_SPath, recursive=TRUE)

        # Format the sigma_SNames into a matrix (nSigmaS * 8)
        sigma_SNames <- paste(rep("sigma_S",8), savedStats,sep=",")
        sigma_SNames <- file.path(sigma_SPath, paste(sigma_SNames,"CEL",sep="."))

        # If one or more of the files don't exist, create any missing ones.
        if (length(which(!file.exists(sigma_SNames)))!=0) {
            for (j in 1:8) { # There are 8 statistics to keep
                if (!file.exists(sigma_SNames[j])) { # If the file doesn't exist
                    # Update any chip-specific header entries
                    hdr <- createCelHeader(filename=sigma_SNames[j], cdf=monoCdf, tags=tags)
                    createCel(sigma_SNames[j],hdr, overwrite=TRUE)
                }
            }
        }

        # Get the sigma_S values from the BMEA.Batch object as an matrix
        # Each of the 8 stats can be accessed by the column, i.e. rHat = c[,8]
        sigma_S <- c()
        for (i in 1:length(units)) {
            rows <- grep("sigma_S",rownames(batchData$summaries[[i]]))
            sigma_S <- rbind(sigma_S, batchData$summaries[[i]][rows,1:8])
        }
        dimnames(sigma_S) <- list(batchData$units$unitNames,savedStats)

        indices <- monoUgcMap$cell[match(as.character(batchData$units$unitNames),monoUgcMap$unit)] # This will update the first group for each unit
        for (j in 1:8){
            # The sigma_S values should be exponentiated as the parameter is on the log-scale
            updateCel(sigma_SNames[j], indices=indices, intensities=exp(sigma_S[,j]))
        }
    }

    # c (Chip Expression Levels)
    # There will only be one value per unit, so CEL files can be monocell CEL files
    # The data can be exponentiated. There may be an issue with non-expressed genes
    # There will be nbrOfArrays*8 CEL files (Mean, StDev, 2.5%, 25%, 50%, 75%, 97.5%, rHat)
    # Store in /bmeaData/modelData/parentName,tags,c/chipType
    # Name the cel files celName,c,Mean.CEL etc.
    # Must be written!!!
    if (writeC) {

        # Check & create the directory
        tags <- paste(parentName,"c",sep=",")
        cPath <- file.path(modelPath,tags,chipName)
        if (!file.exists(cPath)) dir.create(cPath, recursive=TRUE)

        # Format the cNames into a matrix (nChips * 8)
        cNames <- paste(rep(celNames, each=8), savedStats,sep=",")
        cNames <- file.path(cPath, paste(cNames,"CEL",sep="."))
        cNames <- matrix(cNames, nrow=nChips, ncol=8, byrow=TRUE)

        # If one or more of the files don't exist, create any missing ones.
        if (length(which(!file.exists(cNames)))!=0) {
            for (i in 1:nChips) { # Step through all the chips
                for (j in 1:8) { # There are 8 statistics to keep
                    if (!file.exists(cNames[i,j])) { # If the file doesn't exist
                        hdr <- createCelHeader(filename=cNames[i,j], cdf=monoCdf, tags=tags)
                        # Create the monocell .CEL file
                        createCel(cNames[i,j],hdr, overwrite=TRUE)
                    }
                }
            }
        }

        # Get the c values from the BMEA.Batch object as an array
        # Each original chip can be accessed by the 2nd dimension, i.e. cel 1 = c[,1,]
        # Each of the 8 stats can be accessed by the 3rd dimesion, i.e. rHat = c[,,8]
        c <- c()
        for (i in 1:length(units)) {
            rows <- grep("c",rownames(batchData$summaries[[i]]))
            c <- abind(c, array(batchData$summaries[[i]][rows,1:8],dim=c(1,nChips,8)),along=1)
        }
        dimnames(c) <- list(batchData$units$unitNames,celNames,savedStats)

        indices <- monoUgcMap$cell[match(as.character(batchData$units$unitNames),monoUgcMap$unit)] # This will update the first group for each unit
        for (i in 1:nChips) {
            for (j in 1:8){
                # Exponentiate the values as this parameter is on the log scale
                updateCel(cNames[i,j], indices=indices, intensities=exp(c[,i,j]))
            }
        }
    }

    # mu (Condition-specific Expression Levels)
    # There will be one value per unit so use monocell CEL files
    # The data can be exponentiated
    # There will be nConditions*8 CEL files
    # Store in /bmeaData/modelData/parentName,tags,mu/chipType
    # Name the cel files condition,mu,Mean.CEL etc.
    # Should be written!!!
    if (writeMu) {

        # Check the conditions vector:
        if (length(conditions)!=nChips) return(cat("The conditions vector doesn't match the celSet.\n"))
        nMu <- length(levels(conditions)) # The number of conditions

        # Check & create the directory
        tags <- paste(parentName,"mu",sep=",")
        muPath <- file.path(modelPath,tags,chipName)
        if (!file.exists(muPath)) dir.create(muPath, recursive=TRUE)

        # Format the muNames into a matrix (nMu * 8)
        muNames <- paste(rep(levels(conditions), each=8), savedStats,sep=",")
        muNames <- file.path(muPath, paste(muNames,"CEL",sep="."))
        muNames <- matrix(muNames, nrow=nMu, ncol=8, byrow=TRUE)

        # If one or more of the files don't exist, create any missing ones.
        if (length(which(!file.exists(muNames)))!=0) {
            for (i in 1:nMu) { # Step through all the chips
                for (j in 1:8) { # There are 8 statistics to keep
                    if (!file.exists(muNames[i,j])) { # If the file doesn't exist
                        hdr <- createCelHeader(muNames[i, j], cdf=monoCdf, tags=tags)
                        # Create the monocell .CEL file
                        createCel(muNames[i,j],hdr, overwrite=TRUE)
                    }
                }
            }
        }

        # Get the mu values from the BMEA.Batch object as an array
        # Each original chip can be accessed by the 2nd dimension, i.e. cel 1 = c[,1,]
        # Each of the 8 stats can be accessed by the 3rd dimesion, i.e. rHat = c[,,8]
        mu <- c()
        for (i in 1:length(units)) {
            rows <- grep("mu",rownames(batchData$summaries[[i]]))[1:nMu] # If sigma_mu is saved, this needs to be removed here
            mu <- abind(mu, array(batchData$summaries[[i]][rows,1:8],dim=c(1,nMu,8)),along=1)
        }
        dimnames(mu) <- list(batchData$units$unitNames,levels(conditions),savedStats)

        indices <- monoUgcMap$cell[match(as.character(batchData$units$unitNames),monoUgcMap$unit)] # This will update the first group for each unit
        for (i in 1:nMu) {
            for (j in 1:8){
                # Exponentiate the values
                updateCel(muNames[i,j], indices=indices, intensities=exp(mu[,i,j]))
            }
        }

    }

    # sigma_mu (Sd around condition-specific mean)
    # Same as sigma_S (except using mu)
    # Not particularly useful
    if (writeSigmaMu) {

        # Check the conditions vector:
        if (length(conditions)!=nChips) return(cat("The conditions vector doesn't match the celSet.\n"))

        # Check & create the directory
        tags <- paste(parentName, "sigma_mu",sep=",")
        sigma_muPath <- file.path(modelPath, tags, chipName)
        if (!file.exists(sigma_muPath)) dir.create(sigma_muPath, recursive=TRUE)

        # Format the sigma_muNames into a matrix (nSigmaS * 8)
        sigma_muNames <- paste(rep("sigma_mu",8), savedStats,sep=",")
        sigma_muNames <- file.path(sigma_muPath, paste(sigma_muNames,"CEL",sep="."))

        # If one or more of the files don't exist, create any missing ones.
        if (length(which(!file.exists(sigma_muNames)))!=0) {
            for (j in 1:8) { # There are 8 statistics to keep
                if (!file.exists(sigma_muNames[j])) { # If the file doesn't exist
                    hdr <- createCelHeader(sigma_muNames[j], cdf=monoCdf, tags=tags)
                    # Create the monocell .CEL file
                    createCel(sigma_muNames[j],hdr, overwrite=TRUE)
                }
            }
        }

        # Get the sigma_mu values from the BMEA.Batch object as an matrix
        # Each of the 8 stats can be accessed by the column, i.e. rHat = c[,8]
        sigma_mu <- c()
        for (i in 1:length(units)) {
            rows <- grep("sigma_mu",rownames(batchData$summaries[[i]]))
            sigma_mu <- rbind(sigma_mu, batchData$summaries[[i]][rows,1:8])
        }
        dimnames(sigma_mu) <- list(batchData$units$unitNames,savedStats)

        indices <- monoUgcMap$cell[match(as.character(batchData$units$unitNames),monoUgcMap$unit)] # This will update the first group for each unit
        for (j in 1:8){
            # Exponentiate the values
            updateCel(sigma_muNames[j], indices=indices, intensities=exp(sigma_mu[,j]))
        }
    }

    # p (Probe Effects)
    # There will be one value per probe so standard CEL files are required
    # The data MUST be exponentiated! (Negative values)
    # There will be 8 CEL files
    # Store in /bmeaData/modelData/parentName,tags,p/chipType
    # Name the celFiles p,Mean.CEL etc...
    # May be useful...
    if (writeP) {

        # Check & create the directory
        tags <- paste(parentName,"p",sep=",")
        pPath <- file.path(modelPath,tags,chipName)
        if (!file.exists(pPath)) dir.create(pPath, recursive=TRUE)

        # Check for the files, as the directory has been checked above:
        pNames <- paste(rep("probeEffects",8), savedStats,sep=",")
        pNames <- file.path(pPath, paste(pNames,"CEL",sep="."))
        if (length(which(!file.exists(pNames)))!=0) { # If one or more of the files don't exist
            for (j in 1:8) { # There are 8 statistics to keep
                if (!file.exists(pNames[j])) {
                    hdr <- createCelHeader(pNames[j], cdf=cdf, tags=tags)
                    createCel(pNames[j],hdr)
                }
            }
        }

        # Get the S values from the BMEA.Batch object as an array
        # Each original chip can be accessed by the 2nd dimension, i.e. cel 1 = S[,1,]
        # Each of the 8 stats can be accessed by the 3rd dimesion, i.e. rHat = S[,,8]
        p <- c()
        nP <- table(ugcMap$unit) # The number of units
        for (i in 1:length(units)) {
            rows <- setdiff(grep("p",rownames(batchData$summaries[[i]])),c(grep("sigma_p",rownames(batchData$summaries[[i]])),grep("phi",rownames(batchData$summaries[[i]]))))
            p <- rbind(p, batchData$summaries[[i]][rows,1:8])
        }
        dimnames(p) <- list(ugcMap$cell,savedStats)

        for (j in 1:8){
            # These MUST be exponentiated as there may be negative values!!!
            updateCel(pNames[j], indices=ugcMap$cell, intensities=exp(p[,j]))
        }
    }

    # sigma_p (Sd for p)
    # Same as sigma_S (except using p)
    # Not particularly useful
    if (writeSigmaP) {

        # Check the conditions vector:
        if (length(conditions)!=nChips) return(cat("The conditions vector doesn't match the celSet.\n"))

        # Check & create the directory
        tags <- paste(parentName,"sigma_p",sep=",")
        sigma_pPath <- file.path(modelPath, tags, chipName)
        if (!file.exists(sigma_pPath)) dir.create(sigma_pPath, recursive=TRUE)

        # Format the sigma_pNames into a matrix (nSigmaS * 8)
        sigma_pNames <- paste(rep("sigma_p",8), savedStats,sep=",")
        sigma_pNames <- file.path(sigma_pPath, paste(sigma_pNames,"CEL",sep="."))

        # If one or more of the files don't exist, create any missing ones.
        if (length(which(!file.exists(sigma_pNames)))!=0) {
            for (j in 1:8) { # There are 8 statistics to keep
                if (!file.exists(sigma_pNames[j])) { # If the file doesn't exist
                    hdr <- createCelHeader(sigma_pNames[j], cdf=monoCdf, tags=tags)
                    # Create the monocell .CEL file
                    createCel(sigma_pNames[j],hdr, overwrite=TRUE)
                }
            }
        }

        # Get the sigma_p values from the BMEA.Batch object as an matrix
        # Each of the 8 stats can be accessed by the column, i.e. rHat = c[,8]
        sigma_p <- c()
        for (i in 1:length(units)) {
            rows <- grep("sigma_p",rownames(batchData$summaries[[i]]))
            sigma_p <- rbind(sigma_p, batchData$summaries[[i]][rows,1:8])
        }
        dimnames(sigma_p) <- list(batchData$units$unitNames,savedStats)

        indices <- monoUgcMap$cell[match(as.character(batchData$units$unitNames),monoUgcMap$unit)] # This will update the first group for each unit
        for (j in 1:8){
            # Exponentiate again
            updateCel(sigma_pNames[j], indices=indices, intensities=exp(sigma_p[,j]))
        }
    }

    # phi (proportion of trancripts with each exon)
    # There will be one value per group, so monocell CEL files are required
    # The data can be exponentiated
    # There will be nConditions*8 monocell CEL files
    # Store in /bmeaData/modelData/parentName,tags,phi/chipType
    # Name the celFiles condition,phi,mean.CEL etc
    # Should be written!!!
    if (writePhi) {

        # Check the conditions vector:
        if (length(conditions)!=nChips) return(cat("The conditions vector doesn't match the celSet.\n"))
        nH <- length(levels(conditions)) # The number of conditions

        # Check & create the directory
        tags <- paste(parentName, "phi",sep=",")
        phiPath <- file.path(modelPath, tags, chipName)
        if (!file.exists(phiPath)) dir.create(phiPath, recursive=TRUE)

        # Format the phiNames into a matrix (nH * 8)
        phiNames <- paste(rep(levels(conditions), each=8), savedStats,sep=",")
        phiNames <- file.path(phiPath, paste(phiNames,"CEL",sep="."))
        phiNames <- matrix(phiNames, nrow=nH, ncol=8, byrow=TRUE)

        # If one or more of the files don't exist, create any missing ones.
        if (length(which(!file.exists(phiNames)))!=0) {
            for (i in 1:nH) { # Step through all the chips
                for (j in 1:8) { # There are 8 statistics to keep
                    if (!file.exists(phiNames[i,j])) { # If the file doesn't exist
                        hdr <- createCelHeader(phiNames[i,j], cdf=monoCdf, tags=tags)
                        # Create the monocell .CEL file
                        createCel(phiNames[i,j],hdr, overwrite=TRUE)
                    }
                }
            }
        }

        # Get the phi values from the BMEA.Batch object as an array
        # Each original chip can be accessed by the 2nd dimension, i.e. cel 1 = c[,1,]
        # Each of the 8 stats can be accessed by the 3rd dimesion, i.e. rHat = c[,,8]
        phi <- c()
        for (i in 1:length(units)) {
            rows <- grep("phi",rownames(batchData$summaries[[i]]))
            nJ <- length(rows)/nH
            phi <- abind(phi, array(batchData$summaries[[i]][rows,1:8],dim=c(nJ,nH,8)),along=1)
        }
        dimnames(phi) <- list(monoUgcMap$group,levels(conditions),savedStats)

        indices <- monoUgcMap$cell # This will update each entire unit
        for (i in 1:nH) {
            for (j in 1:8){
                # These terms do not need to be exponentiated as the phi term is not on the log scale
                updateCel(phiNames[i,j], indices=indices, intensities=phi[,i,j])
            }
        }
    }

    if (verbose>0) cat("Updating model parameter files...done\n")

    #############
    # CONTRASTS #
    #############

    # There will need to be a directory for each contrast
    # Store in /bmeaData/contrastData/parentName,tags,contrastName/chipType

    # There may be issues with white-space in the contrasts names

    # The phiLogFC & logFC can be saved in the same directories as monocell CEL files
    # There will be 9 files for each (Mean, Sd, 2.5, 25, 50, 75, 97.5, maxP, B)
    # Data must be exponentiated

    if (is.null(names(batchData$logFC))) return(cat("The contrasts must be named\n"))
    contNames <- trim(names(batchData$logFC))
    nCont <- length(contNames)
    # These change when writing the contrasts
    # Clearly there will be issues for the B statistics...
    savedStats <- c("mean","sd","2.5%","25%","50%","75%","97.5%","maxP","B")
    if (verbose>0) cat("Updating contrast parameter files...\n")

    #########
    # logFC #
    #########
    tags <- paste(parentName, "logFC",sep=",")

    # Check & create the directory
    logFCPath <- file.path(contrastPath, tags, chipName)
    if (!file.exists(logFCPath)) dir.create(logFCPath, recursive=TRUE)

    # Format the logFCNames into a matrix (nH * 8)
    logFCNames <- paste(rep(contNames, each=9), savedStats,sep=",")
    logFCNames <- file.path(logFCPath, paste(logFCNames,"CEL",sep="."))
    logFCNames <- matrix(logFCNames, nrow=nCont, ncol=9, byrow=TRUE)

    # If one or more of the files don't exist, create any missing ones.
    if (length(which(!file.exists(logFCNames)))!=0) {
        for (i in 1:nCont) {
            for (j in 1:9) { # There are 9 statistics to keep
                if (!file.exists(logFCNames[i,j])) { # If the file doesn't exist
                    hdr <- createCelHeader(logFCNames[i,j], cdf=monoCdf, tags=tags)
                    # Create the monocell .CEL file
                    createCel(logFCNames[i,j],hdr, overwrite=TRUE)
                }
            }
        }
    }

    # This will update the first group for each unit
    indices <- monoUgcMap$cell[match(as.character(batchData$units$unitNames),monoUgcMap$unit)]
    for (i in 1:nCont) {

        # Get the intensity data. Inf can be written to the files, so everything is OK
        logFC <- batchData$logFC[[i]]
        for (j in 1:9) {
            updateCel(logFCNames[i,j], indices=indices, intensities=exp(logFC[,j]))
        }

    }

    ############
    # phiLogFC #
    ############
    tags <- paste(parentName, "phiLogFC",sep=",")

    # Check & create the directory
    phiLogFCPath <- file.path(contrastPath, tags, chipName)
    if (!file.exists(phiLogFCPath)) dir.create(phiLogFCPath, recursive=TRUE)

    # Format the phiLogFCNames into a matrix (nH * 8)
    phiLogFCNames <- paste(rep(contNames, each=9), savedStats,sep=",")
    phiLogFCNames <- file.path(phiLogFCPath, paste(phiLogFCNames,"CEL",sep="."))
    phiLogFCNames <- matrix(phiLogFCNames, nrow=nCont, ncol=9, byrow=TRUE)

    # If one or more of the files don't exist, create any missing ones.
    if (length(which(!file.exists(phiLogFCNames)))!=0) {
        for (i in 1:nCont) {
            for (j in 1:9) { # There are 9 statistics to keep
                if (!file.exists(phiLogFCNames[i,j])) { # If the file doesn't exist
                    hdr <- createCelHeader(phiLogFCNames[i,j], cdf=monoCdf, tags=tags)
                    # Create the monocell .CEL file
                    createCel(phiLogFCNames[i,j],hdr, overwrite=TRUE)
                }
            }
        }
    }

    # This will update the entire group for each unit
    indices <- monoUgcMap$cell
    for (i in 1:nCont) {

        # Get the intensity data. Inf can be written to the files, so everything is OK
        phiLogFC <- batchData$phiLogFC[[i]]
        for (j in 1:9) {
            updateCel(phiLogFCNames[i,j], indices=indices, intensities=exp(phiLogFC[,j]))
        }

    }

    if (verbose>0) cat("Updating contrast parameter files...done\n")
    cat(sprintf("Updated all CEL files for units %i through %i\n",as.integer(units[1]), as.integer(units[length(units)])))

    ############
    # Log File #
    ############

    if (writeLog) {

        # Check to see if the file exists & create it if necessary
        logPath <- ifelse(is.null(nodePath), file.path("bmeaData"), file.path(nodePath,"bmeaData"))
        logFile <- file.path(logPath,paste(parentName,"log.txt",sep="."))
        if (!file.exists(logFile)) write.table(data.frame("units","parameters","date","time"),logFile,sep="\t",row.names=FALSE,col.names=FALSE)

        # Get the list of units & the time the batch was written
        logUnits <- data.frame(units=units,parameters=paste(paramToWrite,collapse=","),date=format(Sys.time(), "%Y-%m-%d"),time=format(Sys.time(), "%H:%M:%S"))

        # Append the latest batch of units to the end of those already there
        write.table(logUnits, logFile, append=TRUE, col.names=FALSE, row.names=FALSE)

    }

    #####################################################################
    # Return the units & the model parameters that were written to disk #
    #####################################################################

    return(list(units=units,parameters=paramToWrite, time=Sys.time()))

}
