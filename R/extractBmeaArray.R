#' @title Extract the data from an AffymetrixCelSetList as an array or matrix
#' 
#' @description Extract the data saved to an AffymetrixCelSetList for a given set of units
#' 
#' @details During the BMEA process, contrast, model parameters & background priors are all saved as AffymetrixCel Files. 
#' This enables simple extraction of saved data for any combination of units & files.
#' 
#' @param celSetList an object of class \code{\link{AffymetrixCelSetList}}
#' @param units the units to be extracted. Can be either a numeric vector or
#'  a character string with the unitNames as on the appropriate cdf. 
#'  If not supplied, all units will be returned.
#' @param which specify which elements of the celSetList should have data extracted from them. 
#'  If not supplied, data from all elements of the celSetList will be returned
#' @param firstOnly if \code{firstOnly=TRUE}, only the first cell from each unit will be extracted. 
#'  If \code{firstOnly=FALSE} all cells from each unit will be extracted.
#' @param log should the data be log-transformed
#' 
#' @return Returns either a matrix or array, depending on which data is requested
#' 
#' @import abind
#' @import aroma.affymetrix
#' 
#' @export
extractBmeaArray <- function(celSetList, units=NULL, which=NULL, firstOnly=TRUE, log=TRUE) {

    # Check the celSetList
    if (class(celSetList)!="AffymetrixCelSetList") stop("The supplied celSetList is not an AffymetrixCelSetList\n")
    nChips <- unlist(lapply(celSetList, nbrOfArrays))
    if (length(unique(nChips))!=1) stop("The number of arrays is not equal across each list element\n")
    else nChips <- nChips[[1]]
    cdfList <- lapply(celSetList, getCdf)
    chkSums <- unlist(lapply(cdfList, getChecksum))
    if (length(unique(chkSums))!=1) stop("Different checksums: the elements of the celSetList appear to have different cdf files\n")
    cdf <- cdfList[[1]] # Set the cdf
    listNames <- names(celSetList)
    unitNames <- getUnitNames(cdf)

    # Check the units:
    if (is.null(units)) {
        idxs <- as.integer(seq(1, nbrOfUnits(cdf), by=1))
    }
    else {
        # If the units have been supplied numerically:
        if (is.numeric(units)) {
            idxs <- as.integer(units)
            if (max(idxs)>nbrOfUnits(cdf)) stop(sprintf("Invalid units: The maximum number is %i\n", nbrOfUnits(cdf)))
        }
        else { # Otherwise check the character strings:
            if (!is.character(units)) stop("Units can only be supplied as numeric or character vectors\n")
            idxs <- match(units, unitNames)
            if(length(which(is.na(idxs)))!=0) stop(sprintf("Invalid unit: Cannot find unit %s\n", as.character(units[which(is.na(idxs))[1]])))
        }
    }

    # And sort out which components of the list to extract:
    if (is.null(which)){
        which <- as.integer(seq(1, length(celSetList), by=1))
    }
    else {
        if (is.numeric(which)) {
            if (max(which)>length(celSetList)) stop(sprintf("Invalid Selection: The celSetList only has %i elements\n", length(celSetList)))
            else which <- as.integer(which)
        }
        else {
            if (!is.character(which)) stop("The 'which' component can only be supplied as a numeric or character vector\n")
            temp <- which
            which <- match(which, listNames) # Set it as an index vector
            if (length(which(is.na(which)))!=0) stop(sprintf("Invalid Selection: Cannot find list element %s\n", temp[which(is.na(which))[1]]))
        }
    }

    # Now the checks have been passed, extract the data:
    ugcMap <- getUnitGroupCellMap(cdf, units=idxs, retNames=TRUE)
    if (firstOnly) cells <- ugcMap$cell[match(unitNames[idxs],ugcMap$unit)]
    else cells <- ugcMap$cell

    # Use the names of the files as the column names:
    celNames <- lapply(celSetList, getFullNames)
    # Check to see if the listNames appear in the celNames:
    chk <- seq(1, nChips*length(celSetList)) %in% as.vector(sapply(names(celSetList), FUN=grep, unlist(celNames)))
    # If they appear in each celFile remove them, otherwise use the complete names for the column names
    if (length(which(!chk))==0) {
        colNames <- substr(celNames[[1]], start=nchar(names(celSetList)[[1]])+2, stop= nchar(celNames[[1]]))
    }
    else {
        colNames <- celNames[[1]]
    }

    # Set the rownames:
    if (firstOnly) rowNames <- unitNames[idxs]
    else rowNames <- paste(ugcMap$unit, ugcMap$group, sep=".")

    # Now extract the data
    if (length(which)>1) {
        out <- lapply(celSetList[which], extractMatrix, cells=cells)
        out <- abind(out, along=3)
        dimnames(out) <- list(rowNames,colNames,listNames[which])
    }
    else {
        out <- extractMatrix(celSetList[[which]], cells=cells)
        dimnames(out) <- list(rowNames, colNames)
    }
    if (log) out <- log(out) # Log transform if required

    out

}