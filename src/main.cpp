#include "main.h"
#include "Parameters.h"
#include "Dispatcher.h"
#include "probability_v17.h"

using namespace std;

//------------------------------------------------
// simulate from simple individual-based model
Rcpp::List sim_falciparum_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  
  // TODO - delete
  Rcpp::List test_list = args["test_list"];
  vector<int> tmp1 = test_list[0];
  vector<int> tmp2 = test_list[1];
  print_vector(tmp1);
  print_vector(tmp2);
  Rcpp::stop("foo1");
  
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
  //dispatcher.print_infection_map();
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  print("simulation completed in", time_span.count(), "seconds\n");
  
  // return list
  return Rcpp::List::create(Rcpp::Named("daily_values") = dispatcher.daily_values.arr,
                            Rcpp::Named("sample_output") = dispatcher.sample_output);
}
