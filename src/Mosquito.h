
#pragma once

#include "Parameters.h"
#include "array_v1.h"

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
  
  // store zygote genotypes. Length is equal to the number of oocysts, and each
  // element is a matrix with one row if clonal or two rows (representing two
  // parents) if not clonal.
  array_3d_int zygotes;
  
  // store recombinant genotype products. Length is equal to the number of
  // oocysts, and each element is a matrix with four rows correseponding to the
  // four recombinant products. These are created on-the-fly as they are needed,
  // and so are sometimes left empty. For example, if parents are clonal there
  // is no need to generate these products as we can take directly from the
  // zygote.
  array_3d_int products;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Mosquito() {};
  
  // methods
  void init(Parameters* param_ptr);
  void death();
  void denovo_infection(int &haplo_ID);
  void new_infection(Host* host_ptr);
  std::vector<std::vector<int>> get_products(int n_products);
  std::vector<int> recombine(std::vector<std::vector<int>> &zygote);
  
};
