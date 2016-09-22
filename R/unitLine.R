#' @title Add the Unit Line to a plot
#' 
#' @description A useful shortcut
#' 
#' @details This will simply draw the line y=x on any existing plot
#' 
#' @param col The colour of the line
#' @param ... Used to pass any other parameters to the function \code{abline}
#' 
#' @seealso \code{\link{abline}}
#' 
#' @examples
#' x <- rnorm(10)
#' y <- rnorm(10)
#' plot(x, y)
#' unitLine()
#' 
#' @export
unitLine <- function(col="blue",...) {
    abline(a=0, b=1, col=col, ...)
}
