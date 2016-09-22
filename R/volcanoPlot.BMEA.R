#' @title Draw a volcano plot from a BMEA array
#' 
#' @description Draws a volcano plot incorporating the possiblity of points with an infinite value on 
#' the y-axis.
#'
#' @details Will draw a volcano plot from the results. Only a single contrast can be plotted, 
#' and any requested points will be highlighted, up to the maximum specified by \code{hLight.n}
#' 
#' @param data the results from analysis with the BMEA model. 
#' Can be an \code{\link{AffymetrixCelSetList}}, \code{array} or \code{matrix}}
#' @param plotVal determines which summary statistic to plot. 
#' Defaults to mean, but can also accept median or any other saved summary statistic named in the array
#' @param cont determines which contrast to plot. Equates to the the third dimension of any BmeaArray provided
#' @param hLight specify any points to be highlighted by unit ID
#' @param hLightCol the colour to use when highlighting points
#' @param hLight.lab.cex the scaling factor for point labelling
#' @param hLight.pt.cex the scaling factor to use when highlighting points
#' @param hLight.n the number of points to highlight
#' @param ... any further parameters to be passed to \code{\link{plot}}
#' @param pch see \code{\link{plot}}
#' @param cex see \code{\link{plot}}
#' @param main see \code{\link{title}}
#' 
#' @seealso \code{\link{extractBmeaArray}}
#' 
#' @export
volcanoPlot.BMEA <- function(data, plotVal="mean", cont=1, hLight=NULL, hLightCol="blue", hLight.lab.cex=0.7, hLight.pt.cex=1, hLight.n=10, pch=16, cex=0.35, main, ...) {
  
  # The data should be a BmeaArray as extracted by extractBmeaArray
  # The plotVal should select either the mean or median
  # Cont can select which of the contrasts to plot by name or number
  # Highlighted points should be specified by unitName (i.e. the value in dimnames(data)[[1]])
  
  # Ensure only the first value for cont is used
  if (length(cont)>1) {
    cont <- cont[1]
    cat("Warning: Only the first requested contrast will be plotted.\n")
  }
  # Ensure a correct value has been requested for the x-axis
  if (!plotVal %in% c("mean", "median", "50%")) stop("Invalid 'plotVal'\n")
  if (plotVal=="median") plotVal <- "50%"
  
  # Check the class of the data
  dataClass <- class(data)
  if (length(which(c("AffymetrixCelSetList", "array", "matrix") %in% dataClass))!=1) stop("Invalid data format. Must be an AffymetrixCelSetList, array or matrix\n")
  
  # Set the default for the tags
  tags <- c()
  
  if (dataClass=="AffymetrixCelSetList") {
    if (cont>length(data)) stop(sprintf("There are only %i contrasts in the data object\n", length(data)))
    # Ensure 'cont' is passed through the function as a numeric argument
    contNames <- names(data)
    if (!is.numeric(cont)) {
      cont <- match(cont, contNames)
      nas <- which(is.na(cont))
      if (length(nas)!=0) {
        cont <- cont[-which(is.na(cont))]
        if (length(cont)==0) stop("Could not find the named contrast\n")
      }
    }
    # Get the celSet info
    celNames <- getFullNames(data[[cont]])
    splitNames <- unlist(strsplit(celNames, split=","))
    params <- splitNames[seq(length(splitNames)/length(celNames), length(splitNames), by=length(splitNames)/length(celNames))]
    contName <- names(data)[cont]
    tempCs <- extract(data[[cont]], which(params %in% c(plotVal, "B")))
    tags <- getTags(tempCs)
    # Extract the data
    ugcMap <- getUnitGroupCellMap(getCdf(tempCs), retNames=TRUE)
    plotData <- extractMatrix(tempCs, ugcMap$cell)
    if (length(grep("phi", tags))!=0) rownames(plotData) <- as.character(ugcMap$group)
    else rownames(plotData) <- as.character(ugcMap$unit)
    colnames(plotData) <- c(plotVal, "B")
    plotData <- log(plotData[which(plotData[,1]!=0),])
  }
  
  if (dataClass=="array") {
    
    contNames <- dimnames(data)[[3]]
    if (!is.numeric(cont)) {
      cont <- match(cont, contNames)
      nas <- which(is.na(cont))
      if (length(nas)!=0) {
        cont <- cont[-which(is.na(cont))]
        if (length(cont)==0) stop("Could not find the named contrast\n")
      }
    }
    if (cont > dim(data)[3]) stop(sprintf("There are only %i contrasts in the data object\n", dim(data)[3]))
    if (!is.null(dimnames(data)[[3]])) contName <- dimnames(data)[[3]][cont]
    if (length(c(plotVal, "B") %in% dimnames(data)[[2]])!=2) stop("Couldn't find the correct columns in the array/matrix\n")
    plotData <- data[,c(plotVal, "B"),cont]
    
  }
  
  if (dataClass=="matrix") {
    if (length(c(plotVal, "B") %in% colnames(data))!=2) stop("Couldn't find the correct columns in the array/matrix\n")
    plotData <- data[,c(plotVal, "B")]
    if (is.numeric(cont)) contName <- paste("Contrast", cont)
    else contName <- cont
    
  }
  
  # The data will now be an Nx2 matrix, so plot it
  plotData[,"B"] <- abs(plotData[,"B"]) # Ensure the B stat is positive
  plot(plotData, pch=pch, cex=cex,...)
  if (missing(main)) {
    if (is.null(tags)) main <- contName
    else main <- paste(contName,": ", tags[length(tags)], sep="")
  }
  title(main=main,...)
  
  # Check the highlight points.
  if (!is.null(hLight)) {
    int <- intersect(hLight, rownames(plotData))
    if (length(int)==0) hLight <- NULL # If there is no match set it to NULL
    else hLight <- int[1:min(hLight.n, length(int))] # Restrict it to maximum permissable (10 by default)
    hLightCol <- rep(hLightCol, length(int))
    posPts <- which(plotData[int,1]>0)
    negPts <- which(plotData[int,1]<0)
    points(plotData[int,], col=hLightCol, pch=pch, cex=cex*hLight.pt.cex)
    if (length(posPts)>0) text(plotData[int[posPts], ], labels=int[posPts], col=hLightCol[posPts], cex=hLight.lab.cex, adj=c(0.8,-0.5)) # Align to the left top
    if (length(negPts)>0) text(plotData[int[negPts], ], labels=int[negPts], col=hLightCol[negPts], cex=hLight.lab.cex, adj=c(0.2,-0.5)) # Align to the right top
  }
  
}