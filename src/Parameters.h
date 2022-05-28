
#pragma once

#include <Rcpp.h>
#include <tuple>

//------------------------------------------------
// parameters of individual-based simulation model
class Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // genetic parameters
  int L;
  double lambda_oocysts;
  double lambda_products;
  double recomb_prob;
  
  // epidemiological parameters
  double a;
  double mu;
  int u;
  int v;
  int g;
  std::vector<double> prob_infection;
  std::vector<double> duration_infection;
  double infectivity;
  int max_infections;
  
  // deme parameters
  int H;
  std::vector<int> seed_infections;
  std::vector<int> M_vec;
  int n_demes;
  
  // migration
  std::vector<std::tuple<int, int, double>> mig_list;
  int n_mig_list;
  
  // demography
  std::vector<double> life_table;
  std::vector<double> age_death;
  std::vector<double> age_stable;
  int n_age;
  
  // run parameters
  int max_time;
  bool report_progress;
  
  // sampling parameters
  std::vector<std::vector<std::pair<int, int>>> sample_list;
  
  // misc parameters
  double prob_v_death;  // daily probability of mosquito death
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters() {};
  Parameters(const Rcpp::List &args);
  
  // methods
  void print_summary();
};
