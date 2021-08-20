
#include "Mosquito.h"
//#include "Host.h"
#include "probability_v11.h"
#include "misc_v10.h"

using namespace std;

//------------------------------------------------
// constructors
void Mosquito::init(Parameters* param_ptr) {
  L = param_ptr->L;
  max_infections = param_ptr->max_infections;
  //n_haplotypes = 0;
  //n_products = 0;
}

//------------------------------------------------
// initialise with single random haplotype
void Mosquito::denovo_infection() {
  /*
  // two possible methods for denovo creation
  int method = 2;
  
  if (method == 1) {  // coin flip at every locus
    haplotype1 = vector<int>(L);
    for (int i = 0; i < L; ++i) {
      haplotype1[i] = rbernoulli1(0.5) ? 1 : 0;
    }
  } else if (method == 2) {  // one coin flip over all loci
    haplotype1 = vector<int>(L, rbernoulli1(0.5));
  } else {
    Rcpp::stop("error in denovo_infection(): method outside range");
  }
  
  n_haplotypes = 1;
  */
}

//------------------------------------------------
// copy over genotype from infectious host
void Mosquito::new_infection(Host* host_ptr) {
  /*
  // choose two haplotypes at random to copy from host
  int rnd1 = sample2(1, host_ptr->n_infective_haplotypes_total);
  int rnd2 = sample2(1, host_ptr->n_infective_haplotypes_total);
  
  // copy first haplotype from host
  int tmp1 = 0;
  for (int i = 0; i < max_infections; ++i) {
    tmp1 += host_ptr->n_infective_haplotypes[i];
    if (tmp1 >= rnd1) {
      haplotype1 = host_ptr->haplotypes[i][tmp1-rnd1];
      n_haplotypes++;
      break;
    }
  }
  
  // if second copy is identical to first then do not bother copying over
  if (rnd1 == rnd2) {
    return;
  }
  
  // otherwise copy second haplotype from host
  tmp1 = 0;
  for (int i = 0; i < max_infections; ++i) {
    tmp1 += host_ptr->n_infective_haplotypes[i];
    if (tmp1 >= rnd2) {
      haplotype2 = host_ptr->haplotypes[i][tmp1-rnd2];
      n_haplotypes++;
      break;
    }
  }
  
  return;
  */
}

//------------------------------------------------
// death
void Mosquito::death() {
  //haplotype1.clear();
  //haplotype2.clear();
  //n_haplotypes = 0;
  //products.clear();
  //n_products = 0;
}


