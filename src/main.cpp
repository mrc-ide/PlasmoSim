#include "main.h"
#include "Parameters.h"
#include "Dispatcher.h"
#include "probability_v11.h"

using namespace std;

//------------------------------------------------
// simulate from simple individual-based model
Rcpp::List sim_falciparum_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract model parameters into separate class
  Parameters parameters(args);
  //parameters.print_summary();
  
  // R functions
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // create simulation dispatcher object
  Dispatcher dispatcher(parameters, update_progress, args_progress);
  
  // carry out simulation
  dispatcher.simulate();
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  print("simulation completed in", time_span.count(), "seconds\n");
  
  // return list
  return Rcpp::List::create(Rcpp::Named("daily_values") = dispatcher.daily_values.arr,
                            Rcpp::Named("sample_demes") = dispatcher.sample_demes.arr,
                            Rcpp::Named("sample_IDs") = dispatcher.sample_IDs.arr,
                            Rcpp::Named("sample_positive") = dispatcher.sample_positive.arr,
                            Rcpp::Named("sample_haplotypes") = dispatcher.sample_haplotypes.arr);
}
