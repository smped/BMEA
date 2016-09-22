#' @title Draws a boxplot of the phi distribution for all exons/groups.
#' 
#' @description Draws a boxplot of the posterior phi distributions for all exons within a transcript.
#' Can be either a contrast (phiLogFC) or a parameter (phi).
#' 
#' @param data can be either an \code{array}, \code{matrix} or \code{\link{AffymetrixCelSetList}}
#' @param unit specifies which unit to plot if \code{data} is an \code{\link{AffymetrixCelSetList}}.
#' Ignored if \code{data} is an \code{array} or \code{matrix}.
#' @param which specifies which condition/cell-type or contrast to plot. 
#' Defaults to the first condition or contrast
#' @param maxB the maximum possible B-statistic from the entire dataset. 
#' Should be \code{log(2*(nChains*(nIter-nBurnin))/nThin+1)} from the \code{mcmcParam} object in the general workspace. 
#' Only used if plotting a contrast.
#' @param bThresh the maximum B-statistic to display in colour
#' @param nProbes the number of probes per group. 
#' Will over-ride any values found on an associated \code{AffymetrixCdfFile} for a given unit
#' @param grpNames the names of the group-level probesets. 
#' Will over-ride any found on an \code{AffymetrixCdfFile}
#' @param featName the feature being plotted. 
#' Should be a cell-type/condition if plotting raw phi values, or the contrast name if plotting a contrast. 
#' Not required if extracting data from an \code{\link{AffymetrixCelSetList}}
#' @param showGroups logical. Determines whether to diplay the groupNames on the x-axis.
#' @param unchCol the default fill colour for the boxplots
#' @param upCol the colour to use for contrast values greater than \code{bThresh}
#' @param downCol the colour to use for contrast values less than  -\code{bThresh}
#' @param ... passes any further specified graphical parameters
#' @param widths the width of the boxes. 
#' Over-rides any values calculated from an \code{AffymetrixCdfFile}
#' @param main the title of the plot
#' @param medLwd the width of the line to draw for the medians
#' @param lty the line-type for the whiskers
#' @param las see \code{\link{plot}}
#' @param cex.x \code{cex.axis} for the x-axis labels
#' @param cex.y \code{cex.axis} for the y-axis labels
#' @param h values at which a horizontal line will be drawn
#' @param abCol the colour of any horizontal-lines
#' @param hLab text characters to denote significantly up/down changed  groups
#' @param col.axis the colour for the main axes & any zero line
#' 
#' @details 
#' Plots the posterior distributions for the exon proportions (phi) as specified by the BMEA model. 
#' Either the distributions for a given cell type, or a given contrast will be plotted.
#' 
#' Whiskers will extend from the 75th to 97.5th quantile and from the 25th to the 2.5th quantile. 
#' The box will display the IQR.
#' 
#' In the case of a contrast, any significant groups will be indicated by default with an asterisk. 
#' These will only appear for groups with a B-statistic equal to the maximum possible value. 
#' Groups with B-statistics less than this value, but above the specified value for \code{bThresh} will be coloured using \code{upCol} or \code{downCol}.
#' 
#' @seealso \code{\link{AffymetrixCelSetList}}, \code{\link{AffymetrixCdfFile}}, \code{\link{boxplot}}
#' 
#' @import aroma.affymetrix
#' 
#' @export
boxplotPhi <- function(data, unit, which=1, maxB=log(6001), bThresh=log(2997.5/3.5), nProbes, 
                       grpNames, featName=NULL, showGroups=TRUE, unchCol="white", upCol="green", downCol="red",
                       ..., widths,  main=NULL, medLwd=2, lty=1, las=2, cex.x=0.9, cex.y=0.9, h=0, 
                       abCol="blue", hLab="*", col.axis=c("black","grey")){


    dataType <- class(data)[1]
    if(!dataType %in% c("array","matrix","AffymetrixCelSetList")){
        stop("Supplied object must be of type 'array', 'matrix' or 'AffymetrixCelSetList'\n")
    }

    if (dataType=="AffymetrixCelSetList") {
        # Extract the key data, i.e. The matrix of means, quantiles, B & nProbes, groupNames, unitName & featureName
        if (which>length(data)) stop(sprintf("Invalid choice of 'which'. The maximum permissable is %i\n", length(data)))
        monoCdf <- getCdf(data[[which]])
        # If plotting phi itself, no log transformation will be required.
        celNames <- getFullNames(data[[which]])
        if (!is.na(match("B",strsplit(celNames[length(celNames)], split=",")[[1]]))){ # If it is a contrast celSetList
            log <- TRUE
            plotType <- "contrast"
        }
        else{ # If it is the cellType specific phi Values
            log <- FALSE
            plotType <- "parameter"
        }

        # Check the unit
        if (missing(unit)) stop(sprintf("If extracting data from an %s, the unit must be specified", dataType))
        if (length(unit)>1) {
            unit <- unit[1]
            message("More than one unit specified. All apart from the the first will be ignored\n")
        }
        if (is.character(unit)) { # If a unit is specified by name
            unitName <- unit
            unit <- grep(unitName, getUnitNames(monoCdf))
            if (length(unit)==0) stop(sprintf("Invalid Unit: %s\n", unitName))
        }
        else{ # Otherwise just use the number
            if (unit>nbrOfUnits(monoCdf)) stop("Invalid unit\n")
            unitName <- getUnitNames(monoCdf, unit)
        }

        # Get the key info from the monocell cdf
        ugcMap <- getUnitGroupCellMap(monoCdf, unit, retNames=TRUE)
        if (missing(grpNames)) grpNames <- as.character(ugcMap$group)
        if (missing(nProbes)) nProbes <- rep(1, length(grpNames))
        # Now the full cdf
        fullCdfName <- gsub(",monocell", "", getFilename(monoCdf))
        fullCdfPath <- file.path("annotationData", "chipTypes", getName(monoCdf), fullCdfName)
        if (file.exists(fullCdfPath) || missing(nProbes)) nProbes <- readCdfNbrOfCellsPerUnitGroup(fullCdfPath, units=unit)[[1]]
        # Get the feature name & actual values
        featName <- names(data[which])
        phiVals <- extractMatrix(data[[which]], ugcMap$cell)
        # Sort out the column names.
        nCol <- length(unlist(strsplit(celNames[1], split=",")))
        colnames(phiVals) <- unlist(strsplit(celNames, split=","))[seq(nCol, by=nCol, length.out=length(celNames))]
        if (log) phiVals <- log(phiVals) # Log Transform if required (contrasts)
    }

    if (dataType=="matrix") phiVals <- data

    if (dataType=="array") {
        if (which>dim(data)[3]) stop(sprintf("Invalid choice of 'which'. The maximum permissable is %i\n", dim(data)[3]))
        phiVals <- data[,,which]
        dataType <- "matrix"
    }

    # The checks will be the same
    if (dataType=="matrix") {
        if (nrow(phiVals)>300) stop(sprintf("%i groups? The dataset appears to be too large for a single unit\n",nrow(phiVals)))
        if (missing(nProbes)) nProbes <- rep(2, nrow(phiVals))
        if (missing(grpNames)) grpNames <- rownames(phiVals)
        if (!missing(unit)) {
            if (is.character(unit)) unitName <- unit
            if (is.numeric(unit)) unitName <- paste("unit", unit[1])
        }
        else unitName <- c()

        if (length(which(c("mean", "2.5%", "25%", "50%", "75%", "97.5%", "b") %in% tolower(colnames(phiVals))))==7) plotType="contrast"
        if (length(which(c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "rhat" ) %in% tolower(colnames(phiVals))))==8) plotType="parameter"
    }

    # Set the colours
    fillCols <- rep(rgb(t(col2rgb(unchCol)), maxColorValue=255), nrow(phiVals)) # Sets the unchanged colour in hexadecimal

    # SAort out the widths
    if (missing(widths)){
        widths <- nProbes-1 # The widths
    }
    else {
        if (!is.numeric(widths)) stop("widths must be numeric\n")
        if (length(widths)!=nrow(phiVals)) widths <- rep(widths, nrow(phiVals))[1:nrow(phiVals)]
    }

    # The various x co-ordinates
    xl <- c(0, cumsum(widths+1))[-(nrow(phiVals)+1)]+1 # x-left
    xr <- xl+widths # x-right
    xm <- xl+widths/2 # The middle

    # Save the original margins
    orig.mar <- par()$mar

    # Set any non-expressed exons to have different coloured labels:
    noExp <- which(phiVals[,"50%"]==0)
    if (length(col.axis)==1) col.axis <- rep(col.axis, 2)

    # Now draw the plot
    plot.new()
    if (showGroups) par(mar=c(9.5,4,4,2)+0.1)
    plot.window(xlim=c(0,max(xr)), ylim=range(phiVals[,c("2.5%","97.5%")]))
    box()
    segments(xm-widths/3, phiVals[,"2.5%"], xm+widths/3, ...)
    segments(xm-widths/3, phiVals[,"97.5%"], xm+widths/3, ...)
    segments(xm, phiVals[,"2.5%"], xm, phiVals[,"25%"], lty=lty, ...)
    segments(xm, phiVals[,"75%"], xm, phiVals[,"97.5%"], lty=lty, ...)
    if (showGroups){
        labels <- grpNames
    }
    else {
        labels <- 1:nrow(phiVals)
        title(xlab="Group")
    }

    if (length(noExp)!=0) {
        axis(1, at=xm[-noExp], labels=labels[-noExp], las=las, cex.axis=cex.x, col.axis=col.axis[1], ...)
        axis(1, at=xm[noExp], labels=labels[noExp], las=las, cex.axis=cex.x, col.axis=col.axis[2], ...)
    }
    else{
        axis(1, at=xm, labels=labels, cex.axis=cex.x, col.axis=col.axis[1], ...)
    }
    axis(2, las=las, cex.axis=cex.y, ...)
    abline(h=h, col=abCol,...)

    if (plotType=="contrast") {

        # Find any others where the 95% CCI doeszn't overlap zero
        phiVals[which(abs(phiVals[,"B"])<bThresh),"B"] <- 0 # Set any values to zero that are below the significance threshold
        alpha <- 255*(phiVals[,"B"]/maxB)^2 # The transparency for the colours
        upRgb <- t(col2rgb(upCol))
        dnRgb <- t(col2rgb(downCol))
        posB <- which(phiVals[,"B"]>0)
        negB <- which(phiVals[,"B"]<0)
        if (length(posB)>0) fillCols[posB] <- rgb(upRgb[,1], upRgb[,2], upRgb[,3], alpha[posB], maxColorValue=255)
        if (length(negB)>0) fillCols[negB] <- rgb(dnRgb[,1], dnRgb[,2], dnRgb[,3], alpha[negB], maxColorValue=255)

        # Add the boxes...
        rect(xl, phiVals[,"25%"], xr, phiVals[,"75%"], col=fillCols, ...)
        # Add the medians
        segments(xl, phiVals[,"50%"], xr, lwd=medLwd)

        if (is.null(main)) {
            if (is.null(unitName)) {
                if (is.null(featName)) mn <- substitute(expression(paste("Change in log ", phi)))
                else mn <- substitute(expression(paste("Change in log ", phi, " for ", featName)), list(featName=featName))
            }
            else {
                if (is.null(featName)) mn <- substitute(expression(paste(unitName, ": Change in log ", phi)), list(unitName=unitName))
                else  mn <- substitute(expression(paste(unitName, ": Change in log ", phi, " for ", featName)), list(featName=featName, unitName=unitName))
            }
            title(main=eval(mn), ylab=expression(paste(Delta," log ",phi)))
        }
        else {
            title(main=main, ..., ylab=expression(paste(Delta," log ",phi)))
        }

        # Now add asterisks for the points at the extreme (maxB)
        exPts <- which(abs(phiVals[,"B"])>=maxB)
        if (length(exPts)>1) {
            text( xm[exPts], sign(phiVals[exPts,"B"])*rowMaxs(abs(phiVals[exPts,c("2.5%","97.5%")])), labels=hLab, pos=sign(phiVals[exPts,"B"])*1+2, offset=0.1)
        }
        if (length(exPts)==1) {
            text( xm[exPts], sign(phiVals[exPts,"B"])*max(abs(phiVals[exPts,c("2.5%","97.5%")])), labels=hLab, pos=sign(phiVals[exPts,"B"])*1+2, offset=0.1)
        }

    }

    # Otherwise plot the condition specific values
    if (plotType=="parameter") {

        # Add the boxes...
        rect(xl, phiVals[,"25%"], xr, phiVals[,"75%"], col=fillCols, ...)
        # Add the medians
        segments(xl, phiVals[,"50%"], xr, lwd=medLwd)

        # Add the title
        if (is.null(main)) {
            if (is.null(unitName)) {
                if (is.null(featName)) mn <- substitute(expression(phi))
                else mn <- substitute(expression(paste(phi, " for ", featName)), list(featName=featName))
            }
            else {
                if (is.null(featName)) mn <- substitute(expression(paste(unitName, ": ", phi)), list(unitName=unitName))
                else  mn <- substitute(expression(paste(unitName, ": ", phi, " for ", featName)), list(featName=featName, unitName=unitName))
            }
            title(main=eval(mn), ylab=expression(phi))
        }
        else {
            title(main=main, ..., ylab=expression(phi))
        }
    }

    par(mar=orig.mar) # Reset the margins

}
