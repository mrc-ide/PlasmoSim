#------------------------------------------------
#' @title Import file
#'
#' @description Import file from the inst/extdata folder of this package
#' 
#' @param name name of file.
#'
#' @export

plasmosim_file <- function(name) {
  
  # load file from inst/extdata folder
  name_full <- system.file("extdata/", name, package = 'PlasmoSim', mustWork = TRUE)
  ret <- readRDS(name_full)
  
  # return
  return(ret)
}

# -----------------------------------
# takes matrix as input, converts to list format for use within Rcpp code
#' @noRd
matrix_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# takes list format returned from Rcpp and converts to matrix
#' @noRd
rcpp_to_matrix <- function(x) {
  ret <- matrix(unlist(x), nrow = length(x), byrow = TRUE)
  return(ret)
}

#------------------------------------------------
# update progress bar
#' @noRd
update_progress <- function(pb_list, name) {
  knitrProgressBar::update_progress(pb_list$pb_env[[name]])
}
