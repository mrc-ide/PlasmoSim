#------------------------------------------------
#' @title PlasmoSim
#'
#' @description A basic Plasmodium simulator. Contains functions to simulate
#'   epidemiological and genetic data from a simple model of Plasmodium
#'   falciparum transmission.
#'
#' @docType package
#' @name goodegg
NULL

#------------------------------------------------
# link to Rcpp
#' @useDynLib PlasmoSim, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#------------------------------------------------
# unload dll when package is unloaded
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("PlasmoSim", libpath)
}
