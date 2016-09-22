#' @title Make a list object where each element is an AffymetrixCelSet
#' 
#' @description Used for storing values for parameters and contrasts under the BMEA model
#' 
#' @return A list of AffymetrixCelSet objects
#' 
#' @param celSet The parent celSet from which the stored parameters are derived
#' @param type The type of values stored. 
#' Can take the values "backgroundPriors", "contrastData" or"modelData"
#' @param tags Specifies which tags have been added to the parentName.
#' See Details for more information
#' @param path Where to look for the data
#' 
#' @details 
#' For each type of celSetList, only a limited number of tags are valid:
#' \itemize{
#' \item{type=backgroundPriors}{tags can be QN, MAT, GC, RMA, GC-RMA or the like}
#' \item{type=contrastData}{tags can only be logFC of phiLogFC}
#' \item{type=modelData}{tags can be any of the model parameter names (i.e. S, sigmaS, c, mu, sigmaMu, p, sigmaP, phi)}
#' }
#' 
#' By default, all files found with the correct specification will be included in the list object
#' 
#' The cdf will be assigned automatically, 
#' i.e. BG lists will use the same as the celSet, 
#' contrasts will use the monocell CDF & model parameters will use the full (S & p) or monocell (sigmaS, c, mu, sigmaMu, sigmaP, phi) respectively
#' 
#' @import aroma.affymetrix
#' 
#' @export
AffymetrixCelSetList <- function(celSet, type, tags, path=NULL) {
  
  ###################################################
  # Check that the type is one of the 3 permissable #
  ###################################################
  if (missing(type)) {stop("Error: type must be specified\n")}
  else {
    types <- c("contrastData","modelData","backgroundPriors")
    ind <- grep(type, types)
    if(length(ind)!=1) return(cat("Error: Invalid type\n"))
    type <- types[ind]
  }
  
  #########################################################
  # Check the class of the celSet & setup the MonocellCdf #
  #########################################################
  if (length(grep("AffymetrixCelSet",class(celSet)))==0) stop("Error: Supplied celSet is not a celSet.\n")
  parentName <- getParentName(celSet)
  cdf <- getCdf(celSet)
  chipType <- getChipType(cdf)
  cdfPath <- getPath(cdf)
  rootCdfType <- getName(cdf)
  # Get the monocell cdf
  if (!isMonocellCdf(cdf)) {
    monoName <- AffymetrixCdfFile$findByChipType(chipType, tags="monocell")
    if (!is.null(monoName)){
      monoCdf <- AffymetrixCdfFile(monoName)
    }
    else{
      monoCdf <- createMonoCell(cdf)
    }
  }
  else {
    monoCdf <- cdf
  }
  
  ##################
  # Check the tags #
  ##################
  if (missing(tags)) stop("Error: tags must be specified\n")
  # Add them to the parentName if supplied:
  if (!is.null(tags)) {
    parentName <- paste(parentName, tags, sep=",")
  }
  
  ################
  # Set the path #
  ################
  if (type!="backgroundPriors") {
    if (is.null(path)) {
      celPath <- file.path("bmeaData", type, parentName, rootCdfType)
    }
    else {
      celPath <- file.path(path, "bmeaData", type, parentName, rootCdfType)
    }
    # As a quirk of dumb programming, the lambda & delta for the background priors are in separate folders.
    # For contrast & model cel files, they are all in the same directory
    
    if (!file.exists(celPath)) return(cat(sprintf("Error: Invalid directory '%s'\n",celPath)))
    celNames <- dir(celPath)
    celNames <- celNames[grep(".cel", tolower(celNames))] # Make sure that only celFiles are in the list
    if (length(celNames)==0) return(cat(sprintf("No CEL files found in '%s'\n", celPath)))
  }
  else {
    if (is.null(path)) { # If no path has been specified
      if (!file.exists(file.path("bmeaData","backgroundPriors"))) { # Check the bmeaData directory
        
        if (file.exists("backgroundPriors")) celPath <- file.path("backgroundPriors", parentName)
        else return(cat("Cannot find the backgroundPriors directory\n"))
      }
      else celPath <-  file.path("bmeaData","backgroundPriors", parentName)
    }
    else {
      celPath <- file.path(path, "backgroundPriors", parentName)
    }
    if (!file.exists(celPath)) {
      # If the backgroundPriors are in separate folders, the tags 'delta' & 'lambda' may be on the folder names
      # The last folder name will also be the chipType, not the rootCdfType
      celPath <- file.path(paste(celPath, c("lambda","delta"), sep=","),chipType)
      # Check that this path exists
      if (length(which(file.exists(celPath)))!=2) {
        return(cat(sprintf("Error: Invalid directory '%s'\n",celPath)))
      }
    }
  }
  
  
  ####################################################
  # Figure out the keyValues to return in the celSet #
  ####################################################
  if (type=="contrastData") {
    # If returning a contrast celSet, each contrast becomes a separate element
    # Get a vector with the contrast names from the celFiles:
    listNames <- unique(matrix(unlist(strsplit(celNames, split=",")), ncol=2, byrow=TRUE)[,1])
    keyVals <- unique(matrix(unlist(strsplit(celNames, split=",")), ncol=2, byrow=TRUE)[,2])
    keyVals <- substr(keyVals, start=1, stop=nchar(keyVals)-4) # Remove the '.CEL' suffix
    # These keyVals must contain 'mean', 'sd', 'maxP' & 'B' as the minimum
    ord <- match(c("mean","sd","maxP","B"), keyVals)
    if (length(which(is.na(ord)))!=0) return(cat("Missing key value for the specified contrast.\n"))
    keyVals <- keyVals[c(ord[1:2], seq(1, length(keyVals), by=1)[-ord], ord[3:4])] # This will now order them correctly
  }
  if (type=="modelData") {
    # If returning a modelData celSet, each element is a separate parameter
    listNames <- unique(matrix(unlist(strsplit(celNames, split=",")), ncol=2, byrow=TRUE)[,1])
    keyVals <- unique(matrix(unlist(strsplit(celNames, split=",")), ncol=2, byrow=TRUE)[,2])
    keyVals <- substr(keyVals, start=1, stop=nchar(keyVals)-4) # Remove the '.CEL' suffix
    # These keyVals must contain 'mean', 'sd' & 'rHat' as the minimum
    ord <- match(c("mean","sd","rHat"), keyVals)
    if (length(which(is.na(ord)))!=0) return(cat("Missing key value for the specified contrast.\n"))
    keyVals <- keyVals[c(ord[1:2], seq(1, length(keyVals), by=1)[-ord], ord[3])] # This will now order them correctly
  }
  if (type=="backgroundPriors") {
    listNames <- c("lambda","delta")
    keyVals <- getNames(celSet)
  }
  listLen <- length(listNames)
  
  ##################################
  # Now form the celSets as a list #
  ##################################
  out <- vector("list",listLen)
  names(out) <- listNames
  # Set the appropriate cdf as the cdf:
  if (type=="contrastData") cdf <- monoCdf
  if (type=="backgroundPriors") cdf <- cdf
  if (type=="modelData") {
    if (length(which(tags=="S"))!=0 || length(which(tags=="p"))!=0) cdf <- cdf
    else cdf <- monoCdf
  }
  for (i in 1:listLen) {
    if (type=="backgroundPriors") {
      tempFiles <- as.list(paste(file.path(celPath[i],keyVals),"CEL",sep="."))
      for (j in 1:length(tempFiles)) {
        tempFiles[[j]] <- AffymetrixCelFile(tempFiles[[j]], cdf=cdf)
      }
    }
    else {
      tempFiles <- as.list(paste(file.path(celPath,paste(listNames[i], keyVals, sep=",")),".CEL",sep=""))
      for (j in 1:length(tempFiles)) {
        tempFiles[[j]] <- AffymetrixCelFile(tempFiles[[j]], cdf=cdf)
      }
    }
    out[[i]] <- AffymetrixCelSet(tempFiles)
  }
  
  # Give it a class attribute & reuturn the list
  class(out) <- "AffymetrixCelSetList"
  return(out)
  
}