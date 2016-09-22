#' @title Define the MAT model matrix for an given set of sequences
#' 
#' @description Defines the MAT model matrix for a set of sequences.
#' 
#' @details Calculates the model matrix for the MAT background model for a given set of genomic sequences. 
#' Allows fitting of the MAT model for the purposes of background signal estimation.
#' 
#' setMatVector is designed for internal use only.
#' 
#' @param seq genomic sequences. \code{setMatVector} only takes a single string, 
#' whilst \code{setMatMatrix} takes a vector of strings.
#' Each sequence must be a string of 25 characters containing only the letters A, C, G & T
#'  
#' The function setMatVector is designed for one individual sequence.
#' 
#' @return 
#' \code{setMatMatrix} Returns a matrix with nrows=length(seq) & ncol = 80.
#' 
#' \code{setMatVector} nReturns a vector of length 80
#' 
#' @references 
#'   Kapur, K., Xing, Y., Ouyang, Z., Wong, WH. (2007) 
#'   \emph{Exon arrays provide accurate assessments of gene expression} 
#'   Genome Biol. 8(5):R82
#'   
#'   Johnson, W.E., Li, W., Meyer, C.A., Gottardo, R., Carroll, J.S., Brown, M., Liu, X.S. (2006) 
#'   \emph{ Model-based analysis of tiling-arrays for ChIP-chip.}
#'   Proc Natl Acad Sci USA 103:12457-12462
#' 
#' @examples 
#' ## Generate a random sequence
#' seq <- sample(c("A","C","G","T"), 25, replace=TRUE)
#' seq <- paste(seq,collapse="")
#' mat.vec <- setMatVector(seq) ## Returns a vector of length 80
#' 
#' ## Generate a random sequence
#' seq <- replicate(2, sample(c("A","C","G","T"), 25, replace=TRUE), simplify = FALSE)
#' seq <- sapply(seq, paste, collapse = "")
#' mat.matrix <- setMatMatrix(seq) ## Returns a 2x80 matrix
#' 
#' @export
setMatMatrix <- function(seq) {

    if(!is.character(seq)) stop("Invalid Object: 'seq' must be a character vector\n")
    if(length(which(nchar(seq)!=25))!=0) stop("Invalid Sequence Data: All strings must contain 25 characters\n")

    out <- lapply(seq, FUN=setMatVector)
    out <- matrix(unlist(out), ncol=80, byrow=TRUE)
    colnames(out) <-  c("nT", paste("A", 1:25, sep=""), paste("C", 1:25, sep=""), paste("G", 1:25, sep=""), "nAsq", "nCsq", "nGsq", "nTsq")
    out

}
#' @rdname setMatMatrix
#' @export
setMatVector <- function(seq) {
  
  # seq must be a single character string of length 25
  # Returns a vector of length 80 which corresponds to the MAT parameter structure
  
  if (nchar(seq)!=25) stop("Object must be a character string of length 25\n")
  seq <- strsplit(seq,split="")[[1]] # Break the string into 25 single characters
  
  # Defines the model vector for each 25bp probe sequence
  nA <- length(which(seq=="A")) # The number of A
  nC <- length(which(seq=="C")) # The number of C
  nG <- length(which(seq=="G")) # The number of G
  nT <- length(which(seq=="T")) # The number of T
  
  if (nA + nC + nG + nT !=25) stop("Genomic Sequence contains letters which are not A/C/G/T.\n")
  
  iA <- iC <- iG <- rep(0,25) # 3 zero vectors of length 25
  iA[which(seq=="A")] <- 1      # Change to 1 if the base at position n is A
  iC[which(seq=="C")] <- 1      # Change to 1 if the base at position n is C
  iG[which(seq=="G")] <- 1      # Change to 1 if the base at position n is G
  
  out <- c(nT,iA,iC,iG,nA^2,nC^2,nG^2,nT^2)
  return(out)
  
}