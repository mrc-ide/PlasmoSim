
#pragma once

#include "Parameters.h"
#include "array.h"

#include <vector>

//------------------------------------------------
// forward-declare Host class
class Host;

//------------------------------------------------
// class defining mosquito
class Mosquito {
  
public:
  
  // PUBLIC OBJECTS
  
  // copy over some parameters for convenience
  int L;
  int max_infections;
  double lambda_oocysts;
  double recomb_prob;
  
  // number of oocysts
  int n_oocyst;
  
  // genetic data
  array_3d_int zygotes;
  array_3d_int products;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Mosquito() {};
  
  // methods
  void init(Parameters* param_ptr);
  void death();
  void denovo_infection(int haplo_ID);
  void new_infection(Host* host_ptr);
  std::vector<std::vector<int>> get_products(int n_products);
  std::vector<int> recombine(std::vector<std::vector<int>> &zygote);
  
};
