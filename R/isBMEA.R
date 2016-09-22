#' @title Check the class & structure of BMEA objects
#' 
#' @description Check the class & structure of BMEA objects. Generally for internal usage
#' 
#' @details If an object is required for downstream analysis, this functions check the class attribute & the structure of the object to ensure forward compatability.
#' 
#' @param x the object to be tested
#' 
#' @return logical
#' 
#' @export
isBMEA <- function(x) {

    ifelse(length(grep("BMEA",class(x)))!=0,TRUE,FALSE)

}
#' @rdname isBMEA
#' @export
isBMEA.Batch <- function(x) {

    # Tests for the class BMEA.Batch
    if (length(grep("BMEA.Batch",class(x)))==0) return(FALSE)

    # Check that it has the correct components, i.e.: "summaries","logFC","phiLogFC","units","paramToSave","sims"
    if (length(which(c("celSet","summaries","logFC","phiLogFC","conditions","units","paramToSave","sims") %in% names(x)))!=8) {
        cat("The supplied object did not have the correct components.\n")
        return(FALSE)
    }
    else{
        return(TRUE)
    }

}
#' @rdname isBMEA
#' @export
isBMEA.MCMC <- function(x) {

    # Check that it has the correct class
    if (length(grep("BMEA.MCMC",class(x)))==0) return(FALSE)
    # If the class is correct, check that it has the correct structure
    if (length(which(c("mu[1]","phi[1,1]") %in% dimnames(x)[1]))!=2) {
        cat("The mu and/or phi parameters are no correct.\n")
        return(FALSE)
    }
    else {
        return(TRUE)
    }

}