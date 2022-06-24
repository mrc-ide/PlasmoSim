
#include "misc_v14.h"

#include <Rcpp.h>
#include <vector>
#include <chrono>

//------------------------------------------------
// simulate from simple individual-based model
// [[Rcpp::export]]
Rcpp::List sim_falciparum_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress);
