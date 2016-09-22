#' @title Export a BED file(s) for a given set of units
#' 
#' @description Enables easy uploading of probeset structure & AS events to a genome browser
#' 
#' @details 
#' This function enables easy display of the unit & group (probeset) structure of a given gene, or set of genes. 
#' In the case of custom generated CDF files this can be very useful for displaying the structure.
#'   
#' If scores are not provided the displayed tracks will default to black,
#' however if scores are provided, any probesets with score above the threshold will be displayed in 
#' colours of varying intensity, using those supplied in \code{upCol} & \code{downCol}.
#' 
#' All files with the genomic co-ordinates must also be present in the same directory as the supplied cdf. 
#' The name of this file must end with "mapping.txt" and must also contain the full cdf name somewhere in it's name.
#' 
#' @param cdf an \code{AffymetrixCdfFile}
#' @param units the units to be exported
#' @param scoreCel optional. An AffymetrixCelFile containing the group-level, (exon-level) scores for the given units.
#' @param thresh The value above which the scores are indicated using \code{upCol} & \code{downCol}
#' @param singleBed logical. If \code{TRUE}, all units will be output as a single BED file, 
#' otherwise each unit will be written as an individual BED file
#' @param upCol The colour to display for probesets with a positive score
#' @param downCol The colour to display for probesets with a negative score
#' @param ucscTracks The additional tracks to display in the UCSC browser
#' @param ... not used
#'
#' 
#' @return 
#' Returns \code{TRUE} invisibly
#'   
#' Files will be written to the current working directory in the .BED format. 
#' If \code{singleBed=TRUE} a single file will be written including the current cdf name & any tags supplied via the scoreCel file. 
#' If \code{singleBed=FALSE} an individual .BED file will be written for each unit. 
#' The name will consist of any tags on the scoreCel file plus the unit ID as annotated on the cdf.
#'
#' 
#' @seealso 
#' \code{\link{AffymetrixCelFile}}, \code{\link{AffymetrixCdfFile}} 
#' 
#' @export 
units2bed <- function(cdf, units, scoreCel=NULL, thresh=7.3, singleBed=TRUE, upCol="green", downCol="red", ucscTracks=c("ensGene", "knownAlt", "refGene"), ...) {

    # Check the cdf
    if (length(grep("AffymetrixCdfFile",class(cdf)))==0) return(cat("The supplied cdf is incorrect. See help(\"AffymetrixCdfFile\")\n"))

    # Check the units
    if (missing(units)) return(cat("The units must be supplied\n"))
    if (!is.numeric(units)) return(cat("Units must be supplied numerically\n"))
    if (max(units) > nbrOfUnits(cdf)) cat("Units > %i supplied and will be omitted\n", nbrOfUnits(cdf))
    units <- units[which(units <= nbrOfUnits(cdf))]
    nUnits <- length(units)

    # Get the key Cdf info
    chipType <- getChipType(cdf)
    chipName <- getName(cdf)

    # Load in the mapping file
    mapDir <- getPath(cdf)
    cdfFiles <- dir(mapDir, pattern=chipType) # Find the files with the chipType in their names. This needs to be changed
    mapFiles <- dir(mapDir, pattern="mapping.txt") # Find the files the mapping.txt as the suffix. Edit this for pattern matching
    file <- intersect(mapFiles, cdfFiles)
    if (length(file)!=1) return(cat("Unable to determine the correct mapping.txt file\n"))
    file <- file.path(mapDir, file) # Add the correct pathName
    # Check the structure of the mapFile. Set everything to lowercase for ease of text matching
    mapCols <- tolower(as.character(read.table(file, header=FALSE, sep="\t", nrows=1, stringsAsFactors=FALSE)))
    # The columns that are required are:
    colReqs <- c("unit", "group", "chr", "chr strand", "chr start")
    if (length(which(!colReqs %in% mapCols))!=0) return(cat(sprintf("The mapping is missing the %s column\n",colReqs[which(!colReqs %in% mapCols)[1]])))

    map <- read.table(file, header=TRUE, sep="\t")
    colnames(map) <- tolower(colnames(map)) # Set them to lower case for simplicity
    unitNames <- getUnitNames(cdf, units)
    # Get rid of the extra rows in the map file & drop the missing levels of any factors
    map <- map[map$unit %in% unitNames,]
    map$unit <- map$unit[,drop=TRUE]
    map$group <- map$group[,drop=TRUE]
    map$chr <- map$chr[,drop=TRUE]
    map$chr.strand <- map$chr.strand[,drop=TRUE]
    map$chr.start <- map$chr.start[,drop=TRUE]
    nGroups <- length(levels(map$group))

    # Sort the map into order of units, then groups, then chr.start
    map <- map[order(map$unit, map$group, map$chr.start),]

    # Set to the defaults for the scores & colours
    itemRgb <- "Off"
    rgbCol <- matrix(0, nrow=nGroups, ncol=3, dimnames=list(NULL, c("red","green","blue")))
    ucscScores <- rep(1000, nGroups)
    # Check the scoreCel file & get the scores if they are being used
    if (!is.null(scoreCel)) {
        if (length(grep("AffymetrixCelFile", class(scoreCel)))==0) return(cat("Invalid FileType: The 'scoreCel' file must be an AffymetrixCelFile\n"))
        scoreCdf <- getCdf(scoreCel)
        # Ensure that the scoreCel is based on the monocell version of the nmain cdf
        if (getFullName(scoreCdf) != paste(getFullName(cdf), "monocell", sep=",")) return(cat("The 'scoreCel' file must be the same chipType as the cdf\n"))
        itemRgb <- "On" # Turn on the colouring if the scores check out
        # The score must be in a monocell cdf
        monoCdf <- getMonoCell(cdf)
        ugcMap <- getUnitGroupCellMap(monoCdf, units, retNames=TRUE)
        scores <- as.vector(log(extractMatrix(scoreCel, ugcMap$cell)))
        # Scale the scores to between 0 & 1
        maxScore <- max(abs(scores))
        if (maxScore==Inf)  {
            maxScore <- max(abs(scores)[which(abs(scores)!=Inf)])+1
            scores[which(scores==Inf)] <- maxScore
            scores[which(scores==-Inf)] <- -maxScore
        }

        # The ucsc browser also uses scores between 0 & 1000 so set them.
        # In reality, only the zero scores will be affected if using up & down colours
        ucscScores[which(scores==0)] <- 0  # Set the undetectable exons to 0
        scores[which(abs(scores)<thresh)] <- 0 # Set any scores to zero that are below the threshold for colouring
        scores <- scores/maxScore

        # Set the matrix of colours
        allCol <- colours()
        if (is.null(upCol)) upCol <- "black"
        if (!is.character(upCol) || length(grep(upCol[1], allCol))==0) {
            cat("Invalid 'upCol' argument. Setting to the default as green\n")
            upCol <- "green"
            flush.console()
        }
        if (!is.character(downCol) || length(grep(downCol[1], allCol))==0) {
            cat("Invalid 'downCol' argument. Setting to the default as red\n")
            downCol <- "red"
            flush.console()
        }

        upRgb <- matrix(col2rgb(upCol), ncol=1, nrow=3, dimnames=list(c("red","green","blue"),NULL))
        if (is.null(downCol)) downCol <- "black"
        downRgb <- matrix(col2rgb(downCol), ncol=1, nrow=3, dimnames=list(c("red","green","blue"),NULL))
        rgbCol <- matrix(0, nrow=nGroups, ncol=3, dimnames=list(NULL, c("red","green","blue")))
        rgbCol[which(scores>0),] <- matrix(as.integer(scores[which(scores>0)]%*%t(upRgb)),ncol=3) # Set the colours for the up-scores to be highlighted
        rgbCol[which(scores<0),] <- matrix(as.integer(-scores[which(scores<0)]%*%t(downRgb)),ncol=3)
    }
    rgbCol <- apply(rgbCol, FUN=paste, MARGIN=1, collapse=",") # Turn it into a vector

    # Get the structure of the groups to probes
    groupNames <- as.character(unique(map$group))
    groupRows <- lapply(levels(map$group), FUN=grep, map$group)
    names(groupRows) <- groupNames

    # Get the structure of the units to groups
    unitGroups <- lapply(levels(map$unit), FUN=grep, levels(map$group))
    names(unitGroups) <- levels(map$unit)

    # Set the chrom vectors
    chrom <- rep(paste("chr",as.character(map$chr[match(levels(map$unit), map$unit)]), sep=""), times=unlist(lapply(unitGroups, FUN=length)))
    strand <- as.character(map$chr.strand[match(levels(map$group), map$group)])
    chromStart <- map$chr.start[match(levels(map$group), map$group)]
    chromEnd <- map$chr.start[c(match(levels(map$group), map$group)[-1]-1, nrow(map))]+24

    # The blocks (probes within groups)
    blockCount <- unlist(lapply(groupRows, FUN=length))
    sizeFun <- function(times, x) {
        return(rep(x, times))
    }
    blockSizes <- unlist(lapply(lapply(blockCount, FUN=sizeFun, 24), FUN=paste, collapse=","))
    blockFun <- function(rows, starts) {
        return(paste(starts[rows]-starts[rows[1]], collapse=","))
    }
    blockStarts <- unlist(lapply(groupRows, FUN=blockFun, starts=map$chr.start))
    thickStart <- chromStart # This is identical to the chromStart values
    thickEnd <- chromEnd # Again, this is identical to the chromEnd values

    # Collect the bed file data as a data.frame
    bedOut <- data.frame(chrom, chromStart, chromEnd, groupNames, ucscScores, strand, thickStart, thickEnd, rgbCol, blockCount, blockSizes, blockStarts)

    if (singleBed==TRUE) {

        # Now write the file
        if (is.null(scoreCel)) {
            outFileName <-  paste(chipType,"bed",sep=".") # The filename
        }
        else {
            outFileName <- paste(getParentName(scoreCel), getName(scoreCel), "bed", sep=".")
        }

        browserTracks <- sprintf("browser full %s", paste(ucscTracks,collapse=" ")) # The tracks to display when opened
        trackName <- sprintf("track name=\"%s\" description=\"Affy %s\" visibility=2 itemRgb=\"%s\" useScore=1", chipType, chipType, itemRgb)
        cat(c(browserTracks, trackName), file=outFileName, sep="\n", append=FALSE) # Write the file with just the header
        write.table(bedOut, file=outFileName, append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE) # Now append the bed formatted info to the header
    }
    else {

        # If writing separate bed files for each unit
        for (i in 1:nUnits) {

            if (is.null(scoreCel)) {
                outFileName <- paste(unitNames[i], chipType, "bed", sep=".")
            }
            else {
                outFileName <- paste(unitNames[i], getParentName(scoreCel), getName(scoreCel), "bed", sep=".")
            }

            browserTracks <- sprintf("browser full %s", paste(ucscTracks,collapse=" ")) # The tracks to display when opened
            trackName <- sprintf("track name=\"%s\" description=\"Affy %s\" visibility=2 itemRgb=\"%s\" useScore=1", unitNames[i], chipType, itemRgb)
            cat(c(browserTracks, trackName), file=outFileName, sep="\n", append=FALSE) # Write the file with just the header
            write.table(bedOut[unitGroups[[i]],], file=outFileName, append=TRUE, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE) # Now append the bed formatted info to the header

        }

    }

    return(invisible(TRUE))

}

# cdf <- AffymetrixCdfFile$byChipType(chipType="HuEx-1_0-st-v2", tags="U-ENSG,G-ENSE,v14")
# scoreCel <- AffymetrixCelFile$fromFile(filename="TrCont - ThCont,B.CEL", path=file.path("bmeaData","contrastData", "TS_Exon,QN,phiLogFC", "HuEx-1_0-st-v2"))
# units2bed(cdf, 4:6, scoreCel=scoreCel, thresh=log(0.999/0.001))
