
#include "Mosquito.h"
#include "Host.h"
#include "probability_v17.h"
#include "misc_v14.h"

using namespace std;

//------------------------------------------------
// constructors
void Mosquito::init(Parameters* param_ptr) {
  
  // copy over some parameters for convenience
  L = param_ptr->L;
  max_infections = param_ptr->max_infections;
  lambda_oocysts = param_ptr->lambda_oocysts;
  recomb_prob = param_ptr->recomb_prob;
  
  // initialise number of oocysts
  n_oocyst = 0;
  
}

//------------------------------------------------
// death
void Mosquito::death() {
  
  // reset number of oocysts
  n_oocyst = 0;
  
  // clear genetic data
  zygotes.arr.clear();
  products.arr.clear();
  
}

//------------------------------------------------
// initialise genetic data de novo
void Mosquito::denovo_infection(int &haplo_ID) {
  
  // initialise with single unique clonal zygote
  n_oocyst = 1;
  zygotes.arr.emplace_back(vector<vector<int>>(1, vector<int>(L, haplo_ID)));
  
  // increment haplo_ID
  haplo_ID++;
}

//------------------------------------------------
// copy over genotype from infectious host
void Mosquito::new_infection(Host* host_ptr) {
  
  // should not be possible to have infection on already-infected mosquito
  if (n_oocyst != 0) {
    Rcpp::stop("new infection on mosquito with n_oocyst != 0");
  }
  
  // draw number of oocysts
  n_oocyst = rztpois1(lambda_oocysts);
  
  // special case: if host carries a single haplotype then all oocysts will be
  // identical, therefore may as well only initialise one
  if (host_ptr->get_n_active_sexual() == 1) {
    for (int j = 0; j < max_infections; ++j) {
      if (host_ptr->infection_active[j]) {
        if (host_ptr->haplotypes[j].size() == 1) {
          n_oocyst = 1;
        }
        break;
      }
    }
  }
  
  // initialise slots in zygotes
  zygotes = array_3d_int(n_oocyst);
  
  // for each zygote, choose two parents independently
  for (int i = 0; i < n_oocyst; ++i) {
    
    // choose infection slots
    int slot1 = host_ptr->draw_active_slot();
    int slot2 = host_ptr->draw_active_slot();
    
    // choose haplotypes within slots
    int hap1 = sample2(0, host_ptr->haplotypes[slot1].size() - 1);
    int hap2 = sample2(0, host_ptr->haplotypes[slot2].size() - 1);
    
    // if identical parents then zygote is clonal, hence copy over single
    // haplotype
    if ((slot1 == slot2) && (hap1 == hap2)) {
      
      // copy over haplotype
      zygotes[i].emplace_back(host_ptr->haplotypes[slot1][hap1]);
      
    } else {  // otherwise copy over two parental haplotypes
      
      // copy over haplotypes
      zygotes[i].emplace_back(host_ptr->haplotypes[slot1][hap1]);
      zygotes[i].emplace_back(host_ptr->haplotypes[slot2][hap2]);
    }
    
  }
  
  // initialise slots in products
  products = array_3d_int(n_oocyst, 4);
  
}

//------------------------------------------------
// draw a series of haplotypes by sampling n_products times from the recombinant
// oocyst products. This is done in a smart way so that recombinant products are
// only generated as and when they are needed (i.e. not for clonal zygotes).
vector<vector<int>> Mosquito::get_products(int n_products) {
  
  // special case: if a single oocyst containing a single haplotype
  // (representing clonal parents) then copy over just this one haplotype. This
  // is the only possible outcome irrespective of the value of n_products.
  if ((zygotes.arr.size() == 1) && (zygotes[0].size() == 1)) {
    return zygotes[0];
  }
  
  // create matrix of number of times each oocyst product is sampled
  array_2d_int sample_counts(n_oocyst, 4);
  int samp_remaining = n_products;
  int prod_remaining = n_oocyst * 4;
  for (int i = 0; i < n_oocyst; ++i) {
    for (int j = 0; j < 4; ++j) {
      int n = rbinom1(samp_remaining, 1.0 / double(prod_remaining));
      sample_counts[i][j] = n;
      samp_remaining -= n;
      prod_remaining--;
    }
  }
  
  // create new oocyst products as needed and store sampled haplotypes in ret
  // object
  array_2d_int ret;
  for (int i = 0; i < n_oocyst; ++i) {
    
    // skip if this oocyst never sampled
    if (sum(sample_counts[i]) == 0) {
      continue;
    }
    
    // if zygote is clonal then return parental haplotype with no
    // recombination
    if (zygotes[i].size() == 1) {
      ret.arr.emplace_back(zygotes[i][0]);
      
    } else {  // if zygote is not clonal
      
      // loop through all 4 products
      for (int j = 0; j < 4; ++j) {
        
        // skip if this product never sampled
        if (sample_counts[i][j] == 0) {
          continue;
        }
        
        // if product does not already exist then create it by recombination
        if (products[i][j].size() == 0) {
          products[i][j] = recombine(zygotes[i]);
        }
        
        // add to return object
        ret.arr.emplace_back(products[i][j]);
        
      }
    }
  }
  
  return ret.arr;
}

//------------------------------------------------
// generate recombinant product from two parental haplotypes in zygote
vector<int> Mosquito::recombine(vector<vector<int>> &zygote) {
  
  // walk along sequence swapping between parents according to recombination
  // probability
  bool w = rbernoulli1(0.5);
  vector<int> ret = zygote[0];
  for (int i = 0; i < L; ++i) {
    if (rbernoulli1(recomb_prob)) {
      w = !w;
    }
    if (w) {
      ret[i] = zygote[1][i];
    }
  }
  return ret;
}


