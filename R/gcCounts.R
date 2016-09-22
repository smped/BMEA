#' @title Get the number of GC bases from a batch of sequences.
#' 
#' @description Counts the number of G & C bases in a character vector of sequences
#' 
#' @details Counts the number of G & C bases in a batch of a sequences
#' 
#' @param seq a vector of 25bp sequences
#' 
#' @return Returns a \code{vector} the same length as the input \code{vector} of sequences
#' 
#' @export
gcCounts <- function(seq) {

    if (sum( nchar(seq) != 25) != 0) stop("Some of the supplied sequences are not character strings of 25(bp).\n")

    gcCount.single <- function(seq) {

        # Counts a single sequence
        seq <- unlist(strsplit(toupper(seq),split=""))
        count <- length(which(seq=="C")) + length(which(seq=="G"))
        count
    }

    vapply(seq, gcCount.single, numeric(1))

}

