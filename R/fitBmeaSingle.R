#' @title Fit & summarise the BMEA model for a single gene(unit)
#' 
#' @description This function fits the BMEA model & provides the summarised output for a single gene.
#' 
#' @details 
#' This function runs all the necessary checks & fits the BMEA model for a single unit (i.e. gene).
#' The \code{bgCelSet} must be written to disk beforehand, during the preparatory stages.
#' 
#' The value of \code{zGene} can be set to \code{NULL} if the filtering out of low-expressed genes 
#' is not desired. 
#' This value can be changed to any value to effectively restrict the range of expression values in 
#' the dataset to those with a high signal-to-noise ratio. 
#' The value of \code{zExon} can also be set to \code{NULL} to bypass the filtering of exons which 
#' are not confidently detectable above background.
#' 
#' @param celSet an \code{AffymetrixCelSet} with the data to be fit
#' @param bgCelSet a list with components \code{$lambda} & \code{$delta}. 
#' Each of these must be an \code{AffymetrixCelSet} containing the means & standard deviations for the background signal priors
#' @param unit the unit (i.e. gene) to be fit
#' @param conditions a vector of factors specifying which cell-type/condition each array in the \code{celSet} belongs to
#' @param contMatrix a contrast matrix for the summarised output
#' @param ... used for passing variables such as \code{mcmcParam} to \code{runMCMC.BMEA}
#' @param keepSims logical variable. If \code{TRUE} all sims from the process & contrasts will be kept
#' @param zGene the zScore below which the gene is determined to be not detectable above background
#' @param zExon the zScore below which the exon is determined to be not detectable above background
#' 
#' @return 
#' Returns a list with the following components:
#'  
#' \itemize{
#' \item{$summary}{ a matrix with the summary for all the model parameters, with convergence statistics \code{rHat} & \code{nEff}}
#' \item{$logFC}{ a matrix with a separate row for each contrast. Includes mean, sd & key quantiles, with the values maxP & B. maxP refers to the maximum of the proportion of samples which were >=0 or <=0. B is calculated using the weak inequalities log(p>=0) - log(p<=0) as non-zero values for phiLogFC are possible under a (yet to be implemented) mixture model for phi.}
#' \item{$phiLogFC}{ a list with a component for each contrast. Each component is a matrix as above for logFC, but with each row representing an exon.}
#' \item{$sims}{ if \code{keepSims=TRUE}, this will return a list with components \code{$data} for the simulation output, \code{$logFC} for the sampled values of log fold-change and \code{$phiLogFC} for the exon-level sampled values for the change in phi. Otherwise, returns \code{NULL}.}
#' }
#' 
#' @import aroma.affymetrix
#'
#' @export 
fitBmeaSingle <- function(celSet, bgCelSet, unit, conditions, contMatrix, ..., keepSims=FALSE, zGene=4.265, zExon=1.645) {

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
    if (nbrOfArrays(bgCelSet$lambda)!=nbrOfArrays(bgCelSet$delta)) stop("The nbrOfArrays in the gCelSets do not match.\n")

    # Compare the two celSets
    if (chipType != getChipType(bgCelSet$lambda)) stop("The bgCelSet is a different chipType to the main celSet.\n")
    if (nChips != nbrOfArrays(bgCelSet$lambda)) stop("The nbrOfArrys in the bgCelSet doesn't match that of the main celSet.\n")

    # Check the conditions match the celSets
    if (!is.factor(conditions)) stop("The conditions must be specified as a vector of factors.\n")
    if (length(conditions) != nChips) stop("The number of conditions doesn't match the number of arrays.\n")

    # Check that the contMatrix matches the conditions:
    nCond <- nrow(contMatrix)
    if (nCond!=length(levels(conditions))) stop("The contrast matrix doesn't match the conditions.\n")

    # Get the ugcMap & check the requested unit
    if (length(unit)!=1) {
        warning("This function can only process one unit. Only the first argument will be used.\n")
        unit <- unit[1]
    }
    if (!is.numeric(unit)) stop("The unit must be specified numerically.\n")
    if (unit > nbrOfUnits(cdf)) stop(sprintf("Invalid Unit. The maximum unit number is %i.\n", nbrOfUnits(cdf)))
    ugcMap <- getUnitGroupCellMap(cdf,units=unit,retNames=TRUE)
    exons <- factor(ugcMap$group, levels=unique(as.character(ugcMap$group))) # Set the exons vector & keep the names in the correct order
    unitName <- levels(ugcMap$unit)[1] # Get the unitID

    # Extract the necessary data from the celSets:
    PM <- extractMatrix(celSet, cells=ugcMap$cell, field="intensities")
    lambda <- log(extractMatrix(bgCelSet$lambda, cells=ugcMap$cell, field="intensities"))
    delta <- log(extractMatrix(bgCelSet$delta, cells=ugcMap$cell, field="intensities"))

    # Get the key remaining indices
    nExons <- length(levels(exons)) # nExons
    nProbes <- nrow(PM) # nProbes
    nCont <- ncol(contMatrix)

    # Setup the zScore thrsholds
    if (is.null(zGene)) zGene <- -Inf
    if (is.null(zExon)) zExon <- -Inf

    # Calculate the zScore & exit the function if the gene is not expressed
    z <- zScore(PM, lambda, delta, exons)
    if (z$gene<zGene) { # If the gene is not detectable

        message(sprintf("The gene %s has a zScore of %.2f which is <=%.2f and is likely not expressed in any of the cell types.\n", unitName, z$gene, zGene))

        # Now set up the output to be the same format:
        summary <- matrix(0, nrow=nChips*nProbes + nChips + nCond + nExons*nCond + nProbes + 3, ncol=9)
        colnames(summary) <- c("Mean","StDev","2.5%","25%","50%","75%","97.5%","rHat","nEff")
        summary[,c("rHat","nEff")] <- 1
        rownames(summary) <- c(paste(paste("S[",rep(1:nChips,each=nProbes),sep=""),",",1:nProbes,"]",sep=""), "sigma_S", paste("c[",1:nChips,"]",sep=""), paste("mu[",1:nCond,"]",sep=""), "sigma_mu", paste("p[",1:nProbes,"]",sep=""), "sigma_p", paste(rep(paste("phi[",1:nCond,",",sep=""),each=nExons),1:nExons,"]",sep=""))

        logFC <- matrix(0,nrow=nCont,ncol=9)
        colnames(logFC) <- c("Mean","Sd","2.5%","25%","50%","75%","97.5%","maxP","B")
        rownames(logFC) <- colnames(contMatrix)

        phiLogFC <- vector("list",nCont)
        names(phiLogFC) <- colnames(contMatrix)
        for (i in 1:nCont) {
            phiLogFC[[i]] <- matrix(0,nrow=nExons,ncol=9)
            colnames(phiLogFC[[i]]) <- c("Mean","Sd","2.5%","25%","50%","75%","97.5%","maxP","B")
            rownames(phiLogFC[[i]]) <- levels(exons)
        }

        # Format the output & exit the function
        out <- list(summary=summary, logFC=logFC, phiLogFC=phiLogFC, sims=NULL)
        return(out)

    }

    ###################################
    # Now check the exon level zScore #
    ###################################

    #Test the exon level zScore as well & reduce the object if required
    filter <- which(z$exon < zExon) # The undetectable exons
    if (length(filter)!=0){
        rm <- exons %in% levels(exons)[filter] # The data exons & the PM, lambda & delta matrices to remove
        exonInd <- seq(1,nExons)[-filter] # The exon numbers that are fitted

        # Fit the model
        data <- runMCMC.BMEA(PM[!rm,], conditions, exons[!rm, drop=TRUE], lambda[!rm,], delta[!rm,], ...)

        # Get the summary object
        summary <- summariseChains(data, transform=TRUE)

        # Create a separate phi object for the complete phi dataset
        phiRows <- grep("phi", rownames(summary)) # Find the summary rows where phi is. NB: They will always be the last few rows...
        phiRows <- matrix(phiRows, ncol=nCond, byrow=FALSE)
        tempPhi <- matrix(0, nrow=nCond*nExons, ncol=9) # The default value is zero
        tempPhi[,8:9] <- 1 # The default rHat & nEff values are 1
        for (h in 1:nCond) {
            tempPhi[(h-1)*nExons + exonInd,] <- summary[phiRows[,h],]
        }

        # For total accuracy, the Signal & probe effects rows will also need to be adjusted.
        pRows <- setdiff(grep("p", rownames(summary)), c(grep("sigma", rownames(summary)), phiRows))
        if (length(pRows)==0) model <- "noProbes"
        else model <- "Full"
        sRows <- as.vector(matrix(seq(1, nChips*nProbes), ncol=nChips, byrow=FALSE)[!rm,])

        # Modify the summary object to incorporate all exons
        if (model!="noProbes"){ # If the full model was fitted
            newSummary <- matrix(0, nrow=nChips*nProbes + nChips + nCond + nExons*nCond + nProbes + 3, ncol=9)
            rownames(newSummary) <- c(paste(paste("S[",rep(1:nChips,each=nProbes),sep=""),",",1:nProbes,"]",sep=""), "sigma_S", paste("c[",1:nChips,"]",sep=""), paste("mu[",1:nCond,"]",sep=""), "sigma_mu", paste("p[",1:nProbes,"]",sep=""), "sigma_p", paste(rep(paste("phi[",1:nCond,",",sep=""),each=nExons),1:nExons,"]",sep=""))
            newSummary[,8:9] <- 1
            newSummary[c("sigma_S","sigma_mu","sigma_p"),] <- summary[c("sigma_S","sigma_mu","sigma_p"),] # The variances
            newSummary[nChips*nProbes + nChips + nCond + 2 + seq(1, nProbes)[!rm], ] <- summary[pRows,] # The probe effects
            newSummary[nChips*nProbes + nChips + nCond + nProbes + 3 + seq(1, nCond*nExons), ] <- tempPhi # The phi values
        }
        else{ # If no probe effects were included
            newSummary <- matrix(0, nrow=nChips*nProbes + nChips + nCond + nExons*nCond + 2, ncol=9)
            rownames(newSummary) <- c(paste(paste("S[",rep(1:nChips,each=nProbes),sep=""),",",1:nProbes,"]",sep=""), "sigma_S", paste("c[",1:nChips,"]",sep=""), paste("mu[",1:nCond,"]",sep=""), "sigma_mu", paste(rep(paste("phi[",1:nCond,",",sep=""),each=nExons),1:nExons,"]",sep=""))
            newSummary[,8:9] <- 1
            newSummary[c("sigma_S","sigma_mu"),] <- summary[c("sigma_S","sigma_mu"),] # The variances
            newSummary[nChips*nProbes + nChips + nCond + 2 + seq(1, nCond*nExons), ] <- tempPhi # The phi values
        }
        colnames(newSummary) <- c("Mean","StDev","2.5%","25%","50%","75%","97.5%","rHat","nEff")
        newSummary[sRows,] <- summary[seq(1, length(sRows)),] # The signal estimates
        newSummary[nChips*nProbes + 1 + seq(1, nChips + nCond), ] <- summary[length(sRows) + 1 +seq(1, nChips + nCond),] # c & mu
        summary <- newSummary

        # Get the contrasts
        logFC <- getLogFC(data, contMatrix, keepSims=keepSims)
        phiLogFC.fit <- getPhiLogFC(data, contMatrix, exonNames=levels(exons[!rm, drop=TRUE]), keepSims=keepSims)

        # Edit the phiLogFC object so it includes all exons
        phiLogFC <- vector("list", 2)
        names(phiLogFC) <- c("summary","sims")
        phiLogFC$summary <- vector("list", nCont)
        names(phiLogFC$summary) <- names(phiLogFC.fit$summary)
        if (keepSims) phiLogFC$sims <- phiLogFC.fit$sims
        for (i in 1:nCont) {
            phiLogFC$summary[[i]] <- matrix(0, nrow=nExons, ncol=9)
            phiLogFC$summary[[i]][exonInd,] <- phiLogFC.fit$summary[[i]]
            colnames(phiLogFC$summary[[i]]) <- colnames(phiLogFC.fit$summary[[i]])
            rownames(phiLogFC$summary[[i]]) <- levels(exons)
        }

    }
    else { # If all exons are detectable...

        # Fit the model
        data <- runMCMC.BMEA(PM, conditions, exons, lambda, delta, ...)

        # Summarise the output
        summary <- summariseChains(data, transform=TRUE)

        # Get the contrasts
        logFC <- getLogFC(data, contMatrix, keepSims=keepSims)
        phiLogFC <- getPhiLogFC(data, contMatrix, exonNames=levels(exons), keepSims=keepSims)

    }

    # Setup the output
    if (keepSims) sims <- list(data=data, logFC=logFC$sims, phiLogFC=phiLogFC$sims)
    else sims <- NULL

    out <- list(summary=summary, logFC=logFC$summary, phiLogFC=phiLogFC$summary, sims=sims)
    return(out)

}