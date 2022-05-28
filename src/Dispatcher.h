
#pragma once

#include "Parameters.h"
#include "Host.h"
#include "Mosquito.h"
#include "misc_v14.h"
#include "Sampler_v5.h"
#include "array_v1.h"

#include <set>
#include <map>

//------------------------------------------------
// class defining individual-based simulation model
class Dispatcher  {
  
public:
  
  // PUBLIC OBJECTS
  
  // pointers to inputs
  Parameters* param_ptr;
  Rcpp::Function* update_progress_ptr;
  Rcpp::List* args_progress_ptr;
  
  // make local copies of some parameters
  int n_demes;
  int max_time;
  double a;
  int v;
  double prob_v_death;
  int H;
  double infectivity;
  double mu;
  
  // objects for sampling from probability distributions
  Sampler sampler_age_stable;
  Sampler sampler_age_death;
  Sampler sampler_duration_infection;
  
  // ID of next infection and next haplotype
  int next_infection_ID;
  int next_haplo_ID;
  
  // counts of host types
  int H_total;
  std::vector<int> Sh;
  std::vector<int> Eh;
  std::vector<int> Ih;
  
  // population of human hosts
  std::vector<Host> host_pop;
  int next_host_ID;
  
  // store the integer index of hosts in each deme
  array_2d_int host_index;
  std::vector<std::set<int>> host_infective_index;
  
  // vector for randomly changing the order in which migration is applied
  std::vector<int> mig_order;
  
  // counts of mosquito types
  int M_total;
  std::vector<int> Sv;
  std::vector<int> Ev;
  std::vector<int> Iv;
  
  // number of mosquitoes at various stages
  std::vector<int> n_Ev_death_new;
  array_2d_int n_Ev_death;
  array_2d_int n_Ev_to_Iv;
  
  // population of mosquitoes
  std::vector<Mosquito> mosq_pop;
  
  // store integer index of mosquitoes at various stages
  array_2d_int Sv_index;
  array_3d_int Ev_death;
  array_3d_int Ev_to_Iv;
  array_2d_int Iv_index;
  
  // map of infection_ID pedigree
  std::map<int, std::vector<int>> infection_map;
  
  // objects for storing results
  array_3d_double daily_values;
  Rcpp::List sample_output;
  
  // misc
  std::vector<double> EIR;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher(Parameters &parameters, Rcpp::Function &update_progress, Rcpp::List &args_progress);
  
  // methods
  void simulate();
  void print_infection_map();
  
};
