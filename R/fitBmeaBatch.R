#' @title Fit the BMEA model for a single batch of units
#' 
#' @description Fits the BMEA model sequentially for more than one unit
#' 
#' @details This is the function used to fit the BMEA model to a batch of units (or genes). 
#' Each unit is tested to see if it contains multiple exons, and is expressed detectably 
#' above background before analysis. 
#' For single exon genes, all exon-level terms are omitted from the model, 
#' as the PLM model used for conventional 3' Arrays holds for these genes & can be used
#' with minimal computational effort.
#' 
#' Units that are not fitted are also removed from the output vector of units.
#' 
#' Restricting the parameters to be saved, via the \code{paramToSave} argument can significantly 
#' save the memory requirements for large batches of genes. 
#' This will default to the parameters "c", "mu" & "phi". 
#' The signal parameter "S" is the most demanding on memory resources &
#' is generally not advised to be saved unless it is of specific interest.
#' 
#' @param celSet an \code{AffymetrixCelSet} with the data to be fit
#' @param bgCelSet a list with components \code{$lambda} & \code{$delta}. 
#' Each of these must be an \code{AffymetrixCelSet} containing the means & standard deviations
#' for the background signal priors
#' @param units the units (i.e. genes) to be fit
#' @param conditions a vector of factors specifying which cell-type/condition each array in the \code{celSet} belongs to
#' @param contMatrix a contrast matrix for the summarised output
#' @param ... used for passing further arguments such as \code{mcmcParam} to \code{runMCMC.BMEA}
#' @param paramToSave the model parameters to be saved for downstream processing. 
#' The parameters "c", "mu" & "phi" will always be saved.
#' @param keepSims logical variable. 
#' If \code{TRUE} all sims from the process & contrasts will be kept
#' @param zGene the zScore below which a gene is classified as not detectable above background
#' @param zExon the zScore below which an exon is classified as not detectable above background
#' 
#' @return 
#'  An object of class("BMEA.Batch"), which is a list with the following components:
#'  
#'  \itemize{
#'   \item{$celSet}{ the \code{celSet} being analysed, as supplied to the function}
#'   \item{$summaries}{ a \code{list} with a component for each unit. 
#'   Each component contains the summary statistics for the unit, including the convergence statistics "rHat" & "nEff."}
#'   \item{$logFC}{ a \code{list} with a component for each contrast supplied in \code{contMatrix}. 
#'   Each row contains the summary statistics for a single unit, for that contrast}
#'   \item{$phiLogFC}{ a \code{list} with a component for each contrast supplied in \code{contMatrix}. 
#'   Each row represents an exon (group).}
#'   \item{$conditions}{ the cell-types (or conditions) as factors, as supplied to the function}
#'   \item{$units}{ a \code{data.frame} with the units fit & the  corresponding unitNames.}
#'   \item{$paramToSave}{ the parameters requested to be saved.}
#'   \item{$sims}{ a \code{list} with a component for each unit. 
#'   If \code{keepSims=FALSE}, will return \code{NULL} for each component.}
#'   }
#'   
#' @seealso 
#' \code{\link{fitBmeaSingle}}, \code{\link{writeBmeaBatch}}
#' 
#' @import aroma.affymetrix
#' 
#' @export
fitBmeaBatch <- function(celSet, bgCelSet, units, conditions, contMatrix, ..., paramToSave=c("c","mu","phi"), keepSims=FALSE, zGene=4.265, zExon=1.645) {

    # celSet:-      The celSet with the expression data
    # bgCelSet:-    The celSet with the background priors
    # units:-       Which units to analyse
    # conditions:-  Specifies which array belongs to which cellType/treatment
    # contMatrix:-  The contrast matrix for obtaining logFC & phiLogFC
    # ...:-         Used for passing any information to runMCMC.BMEA such as 'mcmcParam'
    # paramToSave:- Specifies which model terms to save in the output.
    # keepSims:-    Should sims be returned in the output, or just the summary information
    # zGene:-       The threshold below which the expression level is declared not detectable above BG for the entire gene across all arrays
    # zExon:-       The threshold below which the exon is declared not detectable above BG for the entire exon across all arrays

    # Check the celSet is a celSet
    if (class(celSet)[1]!="AffymetrixCelSet") stop("The supplied celSet is not a celSet.\n")
    nChips <- nbrOfArrays(celSet)
    chipType <- getChipType(celSet)
    cdf <- getCdf(celSet)

    # Check the bgCelSet is a list with components $lambda & $delta, each of which is a celSet of the same type & size
    if (length(which(c("lambda","delta") %in% attributes(bgCelSet)$names))!=2) stop("The supplied bgCelSet must have components $delta & $lambda.\n")
    if (class(bgCelSet$lambda)[1]!="AffymetrixCelSet") stop("The supplied bgCelSet$lambda is not a celSet.\n")
    if (class(bgCelSet$delta)[1]!="AffymetrixCelSet") stop("The supplied bgCelSet$delta is not a celSet.\n")
    if (getChipType(bgCelSet$lambda)!=getChipType(bgCelSet$delta)) stop("The chipTypes in the bgCelSet do not match.\n")
    if (nbrOfArrays(bgCelSet$lambda)!=nbrOfArrays(bgCelSet$delta)) stop("The nbrOfArrays in the bgCelSets do not match.\n")

    # Compare the two celSets
    if (chipType != getChipType(bgCelSet$lambda)) stop("The bgCelSet is a different chipType to the main celSet.\n")
    if (nChips != nbrOfArrays(bgCelSet$lambda)) stop("The nbrOfArrys in the bgCelSet doesn't match that of the main celSet.\n")

    # Check the conditions match the celSets
    if (!is.factor(conditions)) stop("The conditions must be specified as a vector of factors.\n")
    if (length(conditions) != nChips) stop("The number of conditions doesn't match the number of arrays.\n")
    # Check that the contMatrix matches the conditions:
    if (nrow(contMatrix)!=length(levels(conditions))) stop("The contrast matrix doesn't match the conditions.\n")
    nCont <- ncol(contMatrix)
    nCond <- nrow(contMatrix)

    # Get the ugcMap for the batch of units
    if (!is.numeric(units)) stop("The units must be specified numerically.\n")
    if (max(units) > nbrOfUnits(cdf)) stop(sprintf("Invalid unit(s). The maximum unit number possible is %i.\n",nbrOfUnits(cdf)))
    units <- unique(units) # Ensure there are no duplicates
    nUnits <- length(units)
    unitNames <- getUnitNames(cdf,units)
    ugcMap <- getUnitGroupCellMap(cdf, units=units, force=TRUE, retNames=TRUE) # Using force=TRUE ensures any old cahced values are not used

    # Check the parameters to save
    if (is.null(paramToSave)) {
        paramToSave <- c("S","sigma_S","c","mu","sigma_mu","p","sigma_p","phi")
    }
    # Change to lower case
    paramToSave <- tolower(paramToSave)
    writeS <- ifelse("s" %in% paramToSave, TRUE, FALSE)
    writeSigmaS <- ifelse("sigma_s" %in% paramToSave, TRUE, FALSE)
    # writeC <- ifelse("c" %in% paramToSave, TRUE, FALSE) # This must always be saved
    # writeMu <- ifelse("mu" %in% paramToSave, TRUE, FALSE) # This must always be saved
    writeSigmaMu <- ifelse("sigma_mu" %in% paramToSave, TRUE, FALSE)
    writeP <- ifelse("p" %in% paramToSave, TRUE, FALSE)
    writeSigmaP <- ifelse("sigma_p" %in% paramToSave, TRUE, FALSE)
    # writePhi <- ifelse("phi" %in% paramToSave, TRUE, FALSE) # This must always be saved too

    # Extract the necessary data from the celSets in one batch:
    PM <- extractMatrix(celSet, cells=ugcMap$cell, field="intensities")
    lambda <- log(extractMatrix(bgCelSet$lambda, cells=ugcMap$cell, field="intensities"))
    delta <- log(extractMatrix(bgCelSet$delta, cells=ugcMap$cell, field="intensities"))

    # Set up the output. Now the logFC will need to be setup for each contrast. The rownames should be the unitID
    # The phiLogFC should be similar except the rownames will be the unitID & exonID joined with sep="."
    out <- vector("list",8)
    names(out) <- c("celSet","summaries","logFC","phiLogFC","conditions","units","paramToSave","sims")
    out$celSet <- celSet
    out$summaries <- vector("list", nUnits) # A new summary for each unit
    out$logFC <- vector("list",nCont)
    names(out$logFC) <- colnames(contMatrix)
    out$phiLogFC <- vector("list",nCont)
    names(out$phiLogFC) <- colnames(contMatrix)
    out$conditions <- conditions
    out$units <- data.frame(units, unitNames)
    out$paramToSave <- paramToSave
    out$sims <- vector("list",nUnits)

    exonRowNames <- c()

    for (i in 1:nUnits) {

        curUnit <- unitNames[i]
        curRows <- which(ugcMap$unit==curUnit)
        exons <- factor(ugcMap$group[curRows], levels=unique(as.character(ugcMap$group[curRows]))) # Set the exons vector
        nProbes <- length(exons) # The number of probes in the current unit
        nExons <- length(levels(exons))
        exonRowNames <- c(exonRowNames, paste(curUnit,levels(exons),sep=".")) # The cumulative exonIDs

        # Check for detectable expression. Also skip any genes with only one probe
        if (nProbes!=1) {
            z <- zScore(PM[curRows,], lambda[curRows,], delta[curRows,], exons)
        }
        else {
            z <- list(gene=-Inf, exon=-Inf)
        }

        if (z$gene <= zGene) {
            # Print a message
            if (z$gene==-Inf){
                cat(sprintf("%s only contains one probe and has been omitted from the analysis.\n", curUnit))
            }
            else {
                cat(sprintf("%s has a zScore of %.2f which is <=%.2f and is likely not expressed in any samples.\n", curUnit,z$gene, zGene))
            }

            # Now set up the output to be the same format:
            data <- NULL
            summary <- matrix(0, nrow=prod(length(curRows)*nChips) + nChips + nCond + nExons*nCond + length(curRows) + 3, ncol=9)
            colnames(summary) <- c("Mean","StDev","2.5%","25%","50%","75%","97.5%","rHat","nEff")
            summary[,c("rHat","nEff")] <- 1
            rownames(summary) <- c(paste(paste("S[",rep(1:nChips,each=length(curRows)),sep=""),",",1:length(curRows),"]",sep=""), "sigma_S", paste("c[",1:nChips,"]",sep=""), paste("mu[",1:nCond,"]",sep=""), "sigma_mu", paste("p[",1:length(curRows),"]",sep=""), "sigma_p", paste(rep(paste("phi[",1:nCond,",",sep=""),each=nExons),1:nExons,"]",sep=""))

            logFC <- vector("list",2)
            names(logFC) <- c("summary","sims")
            logFC$summary <- matrix(0,nrow=nCont,ncol=9)
            colnames(logFC$summary) <- c("Mean","Sd","2.5%","25%","50%","75%","97.5%","maxP","B")
            rownames(logFC$summary) <- colnames(contMatrix)

            phiLogFC <- vector("list",2)
            names(phiLogFC) <- c("summary","sims")
            phiLogFC$summary <- vector("list",nCont)
            names(phiLogFC$summary) <- colnames(contMatrix)
            for (j in 1:nCont) {
                phiLogFC$summary[[j]] <- matrix(0,nrow=nExons,ncol=9)
                colnames(phiLogFC$summary[[j]]) <- c("Mean","Sd","2.5%","25%","50%","75%","97.5%","maxP","B")
                rownames(phiLogFC$summary[[j]]) <- levels(exons)
            }

                # Now add the latest run to the output:
            out$summaries[[i]] <- summary
            for (j in 1:nCont) {
                out$logFC[[j]] <- rbind(out$logFC[[j]],logFC$summary[j,]) # Add the lastest logFC
                out$phiLogFC[[j]] <- rbind(out$phiLogFC[[j]],phiLogFC$summary[[j]]) # Add the latest phiLogFC
            }
            if (keepSims) out$sims[[i]] <- data

        }

        else { # If the gene is expressed in some detectable quantity

            # Setup the process to use the zExons value
            filter <- which(z$exon < zExon) # The undetectable exons
            if (length(filter)!=0){
                rm <- exons %in% levels(exons)[filter] # The data exons & the PM, lambda & delta matrices to remove
                exonInd <- seq(1,nExons)[-filter] # The exon numbers that are fitted
                curRows <- curRows[!rm]
                exons <- exons[!rm,drop=TRUE]
            }

            # Fit the model
            data <- runMCMC.BMEA(PM[curRows,], conditions, exons, lambda[curRows,], delta[curRows,], ...)

            # Summarise the output
            summary <- summariseChains(data, transform=TRUE)

            # Get the contrasts
            logFC <- getLogFC(data, contMatrix, keepSims=keepSims)
            phiLogFC <- getPhiLogFC(data, contMatrix, exonNames=levels(exons), keepSims=keepSims)

            # If any exons have bene removed
            if (length(filter)!=0) {

                newSummary <- c()

                # If exons have been filtered out, setup the output correctly using the 'write' variables
                if (writeS) {
                    tempS <- matrix(0, nrow=nProbes*nChips, ncol=9) # Create the object for the S values
                    tempS[,8:9] <- 1
#                    sRows <- as.vector(matrix( seq(1, nProbes*nChips), nrow=nProbes, byrow=FALSE)[curRows,]) # Set the rows which should have data
                    sRows <- as.vector(matrix( seq(1, nProbes*nChips), nrow=nProbes, byrow=FALSE)[curRows-min(curRows)+1,]) # Set the rows which should have data
                    tempS[sRows,] <- summary[seq(1, length(curRows)*nChips),] # Extract the data form the original summary object
                    newSummary <- rbind(newSummary, tempS)
                    rownames(newSummary) <- paste(rep(paste("S[", seq(1, nChips),",", sep=""), each=nProbes), seq(1, nProbes), "]", sep="")
                }

                # Sigma_S
                if (writeSigmaS) {
                    newSummary <- rbind(newSummary, summary["sigma_S",])
                    rownames(newSummary)[nProbes*nChips + 1] <- "sigma_S"
                }

                # c & mu are always written so add them automatically
                newSummary <- rbind(newSummary, summary[grep("c", rownames(summary)),]) # The chip effects
                newSummary <- rbind(newSummary, summary[grep("mu", rownames(summary))[1:nCond],]) # mu

                # Sigma_Mu
                if (writeSigmaMu) {
                    newSummary <- rbind(newSummary, summary["sigma_mu",])
                    rownames(newSummary)[nrow(newSummary)] <- "sigma_mu"
                }

                # Probe effects
                if (writeSigmaP) {
                    tempP <- matrix(0, nrow=nProbes, ncol=9)
                    tempP[,8:9] <- 1
                    pRows <- setdiff(grep("p", rownames(summary)), c(grep("sigma_p",rownames(summary)),grep("phi",rownames(summary))))
                    tempP[curRows-min(curRows)+1,] <- summary[pRows,]
                    rownames(tempP) <- paste("p[", seq(1, nProbes),"]", sep="")
                    newSummary <- rbind(newSummary, tempP)
                }

                # Sigma_P
                if (writeSigmaP) {
                    newSummary <- rbind(newSummary, summary["sigma_p",])
                    rownames(newSummary)[nrow(newSummary)] <- "sigma_p"
                }

                # Phi is always written
                tempPhi <- matrix(0, nrow=nExons*nCond, ncol=9)
                tempPhi[,8:9] <- 1
                rownames(tempPhi) <- paste(rep(paste("phi[",1:nCond,",",sep=""), each=nExons), 1:nExons, "]", sep="")
                tempPhi[as.vector(matrix(1:(nCond*nExons), nrow=nExons, ncol=nCond)[exonInd,]),] <- summary[grep("phi", rownames(summary)),]
                newSummary <- rbind(newSummary, tempPhi)

                # And now assign everything to the output object
                out$summaries[[i]] <- newSummary

                # Setup the phiLogFC object for omitted exons & assign the logFC
                newPhiLogFC <- vector("list", nCont)
                for (j in 1:nCont) {

                    # Add in the logFC to the output
                    out$logFC[[j]] <- rbind(out$logFC[[j]],logFC$summary[j,]) # Add the lastest logFC
                    # The phiLogFC
                    tempPhi <- matrix(0, nrow=nExons, ncol=9)
                    tempPhi[exonInd,] <- phiLogFC$summary[[j]]
                    out$phiLogFC[[j]] <- rbind(out$phiLogFC[[j]], tempPhi)
                }

                ############################################
                # The rownames still need to be checked!!! #
                ############################################

                if (keepSims) out$sims[[i]] <- data

            }
            else {# If no exons have been filtered out it's very straightforward

                # Strip out any data that is not required. This will significantly save on memory usage for large batches:
                # NB: The parameters c, mu & phi are always saved!
                if (!writeS) summary <- summary[-(1:(length(curRows)*nChips)),]
                if (!writeSigmaS) summary <- summary[-grep("sigma_S",rownames(summary)),]
                if (!writeSigmaMu) summary <- summary[-grep("sigma_mu",rownames(summary)),]
                if (!writeP) summary <- summary[-setdiff(grep("p",rownames(summary)),c(grep("sigma_p",rownames(summary)),grep("phi",rownames(summary)))),]
                if (!writeSigmaP) summary <- summary[-grep("sigma_p",rownames(summary)),]

                # Now add the latest run to the output:
                out$summaries[[i]] <- summary
                for (j in 1:nCont) {
                    ######################################################################################
                    # NB: Getting the rownames sorted for any genes that are not fit, is crucial here!!! #
                    # It's not checked yet either, although those units are omitted from the output      #
                    ######################################################################################
                    out$logFC[[j]] <- rbind(out$logFC[[j]],logFC$summary[j,]) # Add the lastest logFC
                    out$phiLogFC[[j]] <- rbind(out$phiLogFC[[j]],phiLogFC$summary[[j]]) # Add the latest phiLogFC
                }
                if (keepSims) out$sims[[i]] <- data
            }
        }
    }

    # Once all units have been fit, tidy the output
    names(out$summaries) <- unitNames
    names(out$sims) <- unitNames
    for (j in 1:nCont) {
        rownames(out$logFC[[j]]) <- unitNames
        rownames(out$phiLogFC[[j]]) <- exonRowNames
    }

    class(out) <- "BMEA.Batch"
    return(out)

}