
#pragma once

#include "Parameters.h"

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
  
  // haplotypes
  //std::vector<int> haplotype1;
  //std::vector<int> haplotype2;
  //int n_haplotypes;
  //std::vector<std::vector<int>> products;
  //int n_products;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Mosquito() {};
  
  // methods
  void init(Parameters* param_ptr);
  void denovo_infection();
  void new_infection(Host* host_ptr);
  void death();
  //std::vector<int> get_product();
  
};
