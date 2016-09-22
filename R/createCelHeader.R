#' @title Create a CEL header from the supplied inputs
#' 
#' @description Create a list in the same format as a CEL header, 
#' in preparation for creating a CEL file
#' 
#' @param filename the name of the CEL file that the header will be written to
#' @param cdf the cdf for the CEL file
#' @param version the version number of the CEL file to be created
#' @param cols the number of columns for the CEL file. 
#' Ignored if a \code{cdf} is supplied
#' @param rows the number of rows for the CEL file. 
#' Ignored if a \code{cdf} is supplied
#' @param algorithm a text string to be written to the \code{$algorithm} component
#' @param parameters a text string to be written to the \code{$parameters} component
#' @param chiptype a text string defining the chiptype to be  written. 
#' Ignored if a \code{cdf} is supplied
#' @param header a text string to be written to the \code{$header} component. 
#' Must be in the correct format
#' @param datheader a text string to be written to the \code{$datheader} component. 
#' Must be in the correct format
#' @param librarypackage a text string to be written to the \code{$librarypackage} component. 
#' Defaults to ""
#' @param cellmargin a text string to be written to the \code{cellmargin} component
#' @param noutliers the value to be written to the \code{$noutliers} component
#' @param nmasked defaults to zero
#' @param tags An optional text string which will be placed before the filename in the \code{$datheader} component
#' 
#' @details 
#' This function will output a list in the correct format for generation of a blank CEL file using \code{createCel}.
#' 
#' Many of the input values are redundant, but must be included to pass checks during the creation of CEL files.
#' 
#' If a cdf is specified, any values for \code{cols}, \code{rows} or \code{chiptype} will be ignored, and those on the cdf will be used.
#' 
#' @return 
#' Returns a list with the following components:
#' \itemize{
#' \item{$filename}{ the filename of the CEL file to be created}
#' \item{$version}{ the version of the CEL file}
#' \item{$cols}{ the number of rows to be generated on the CEL file}
#' \item{$rows}{ the number of rows to be generated on the CEL file}
#' \item{$total}{ the total number of cells. Will always equal \code{$rows}*\code{$cols}}
#' \item{$algorithm}{ a typical algorithm which is available for scanners}
#' \item{$parameters}{ a typical list of parameters as would be found on a CEL file}
#' \item{$chiptype}{ the type of CEL file to be written}
#' \item{$header}{ the header in the correct format}
#' \item{$datheader}{ a typical header output from a .DAT file}
#' \item{$librarypackage}{ can be any text string}
#' \item{$cellmargin}{ an integer}
#' \item{$noutliers}{ an integer}
#' \item{$nmasked}{ an integer}
#' }
#' 
#' @seealso
#' \code{\link{createCel}} \code{\link{AffymetrixCdfFile}}
#' 
#' @import aroma.affymetrix
#' 
#' @export
createCelHeader <- function(filename=NULL, cdf=NULL, version=4, cols=NULL, rows=NULL, algorithm="Percentile", parameters=NULL, chiptype=NULL, header=NULL, datheader=NULL, librarypackage="", cellmargin=4, noutliers=0, nmasked=0, tags=NULL) {

    # This creates a CEL header in the correct format. There is a fudge or two that will have no bearing on the actual files created

    if (is.null(filename)) stop("The filename must be supplied\n")

    if (!is.null(cdf)) {
        if (length(grep("AffymetrixCdfFile",class(cdf)))==0) stop("The supplied cdf must be a AffymetrixCdfFile\n")
        header <- getHeader(cdf)
        # The following will override any supplied information
        cols <- header$cols
        rows <- header$rows
        chiptype <- header$chiptype
        chipname <- getName(cdf) # The root name of the chipType, e.g. HuEx-1_0-st-v2
    }
    else {
        if (is.null(cols) || is.null(rows)) stop("If no cdf is supplied, the number of rows & columns must be specified\n")
        if (is.null(chiptype)) stop("If no cdf is supplied, the chiptype must be specified explicitly\n")
        chipname <- chiptype

    }

    total <- cols*rows
    if (is.null(parameters)){
        # This is a blatant fudge...
        parameters <- sprintf("Percentile:75;CellMargin:%i;OutlierHigh:1.500;OutlierLow:1.004;AlgVersion:6.0;FixedCellSize:TRUE;FullFeatureWidth:7;FullFeatureHeight:7;IgnoreOutliersInShiftRows:FALSE;FeatureExtraction:TRUE;PoolWidthExtenstion:1;PoolHeightExtension:1;UseSubgrids:TRUE;RandomizePixels:FALSE;ErrorBasis:StdvMean;StdMult:1.000000;NumDATSubgrids:169",as.integer(cellmargin))
    }

    # Specify the datheader in the correct format:
    if (is.null(datheader)) {
        datheader <- sprintf("[0..37616]  %s:CLS=19288RWS=19288XIN=0  YIN=0  VE=30        2.0 %s System  M10   \024  \024 %s.1sq \024  \024  \024  \024  \024 570 \024 25533.753906 \024 3.500000 \024 0.7000 \024 3", ifelse(is.null(tags),paste(chipname,basename(filename),sep="/"),paste(tags,chipname,basename(filename),sep="/")), format(Sys.time(),format="%d/%m/%y %H:%M:%S"),chipname)
    }

    # Specify the header in the correct format:
    # The GridCorners are fudges here, but they are not particulary relevant anyway...
    header <- sprintf("Cols=%i\nRows=%i\nTotalX=%i\nTotalY=%i\nTotal=%i\nOffsetX=0\nOffsetY=0\nGridCornerUL=501 503\nGridCornerUR=18775 463\nGridCornerLR=18787 18665\nGridCornerLL=510 18704\nAxis-invertX=0\nAxisInvertY=0\nswapXY=0\nDatHeader=%s\nAlgorithm=%s\nAlgorithmParameters=%s\n", cols, rows, cols, rows, total, datheader, algorithm, parameters)

    # Spit out the header:
    list(filename=filename,
         version=version, 
         cols=cols, 
         rows=rows, 
         total=total, 
         algorithm=algorithm, 
         parameters=parameters, 
         chiptype=chiptype, 
         header=header, 
         datheader=datheader, 
         librarypackage=librarypackage, 
         cellmargin=cellmargin, 
         noutliers=noutliers, 
         nmasked=nmasked)

}
