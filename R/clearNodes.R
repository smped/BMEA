#' @title Remove all files/folders used on each node
#' 
#' @description After running the BMEA process in parallel, 
#' this function removes all the node-specific files
#' 
#' @param nodePaths a list where each element contains the folder 
#' where the node-specific CEL & log files are stored
#' 
#' @details NB: This must only be run after the nodes have been successfully merged,
#' as all information saved in the node folders will be deleted permanently
#' 
#' @import R.utils
#' 
#' @export
clearNodes <- function(nodePaths){

    for (i in 1:length(nodePaths)){
        if (file.exists(file.path(nodePaths[[i]]))) {
          removeDirectory(nodePaths[[i]], recursive=TRUE)
          message(sprintf("Deleted folder '%s' & all subsequent files.\n",nodePaths[[i]]))
          flush.console()
        }
      else {
        message("Could not locate", file.path(nodePaths[[i]]))
      }
    }

    invisible(TRUE)
}
