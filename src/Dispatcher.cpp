
#include "Dispatcher.h"
#include "probability_v17.h"

#include <tuple>

using namespace std;

//------------------------------------------------
// constructor
Dispatcher::Dispatcher(Parameters &parameters, Rcpp::Function &update_progress, Rcpp::List &args_progress) {
  
  // pointers to inputs
  param_ptr = &parameters;
  update_progress_ptr = &update_progress;
  args_progress_ptr = &args_progress;
  
  // make local copies of some parameters for convenience
  n_demes = param_ptr->n_demes;
  max_time = param_ptr->max_time;
  a = param_ptr->a;
  v = param_ptr->v;
  prob_v_death = param_ptr->prob_v_death;
  H = param_ptr->H;
  infectivity = param_ptr->infectivity;
  mu = param_ptr->mu;

  // objects for sampling from probability distributions
  sampler_age_stable = Sampler(param_ptr->age_stable, 1000);
  sampler_age_death = Sampler(param_ptr->age_death, 1000);
  sampler_duration_infection = Sampler(param_ptr->duration_infection, 1000);
  
  // ID of next infection and next haplotype. Initialise at 0
  next_infection_ID = 0;
  next_haplo_ID = 0;
  
  // counts of host types. Start with H susceptibles in each deme
  H_total = n_demes*H;
  Sh = vector<int>(n_demes, H);
  Eh = vector<int>(n_demes);
  Ih = vector<int>(n_demes);
  
  // initialise single population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single static population
  // we can simply change the "deme" attribute of a host to represent migration
  host_pop = vector<Host>(H_total);
  next_host_ID = 0;
  
  // for each deme, store the integer index of all hosts in that deme, and the
  // integer index of infective hosts only
  host_index = array_2d_int(n_demes);
  for (int k = 0; k < n_demes; ++k) {
    host_index[k] = seq_int(k*H, (k + 1)*H - 1);
    reshuffle(host_index[k]);
  }
  host_infective_index = vector<set<int>>(n_demes);
  
  // initialise hosts
  for (int k = 0; k < n_demes; ++k) {
    for (int i = 0; i < H; ++i) {
      int this_host = host_index[k][i];
      host_pop[this_host].init(this_host, next_host_ID, k,
                               parameters,
                               sampler_age_stable,
                               sampler_age_death,
                               sampler_duration_infection,
                               host_infective_index);
    }
  }
  
  //seed initial infections
  for (int k = 0; k < n_demes; ++k) {
    for (int i = 0; i < param_ptr->seed_infections[k]; i++) {
      host_pop[host_index[k][i]].denovo_infection(next_haplo_ID, param_ptr->seed_vec[i], 0);
    }
  }
  // This works and produces seed infections across two demes of the number we assign in seed_infections
  // e.g. seed_infections = (3,2) and seed_vec randomly (as default) = (1,4,1) this produces
  // first deme: host1 = 1 infection 0, host2 = 4 infections 1,2,3,4, host3 = 1 infection 5
  // second deme: host1 = 1 infecton 6, host2 = 4 infections 7,8,9,10
  
  // This is trying to use the list and causes an abort. Changes are:
  // 1. In the loop instruction ++k to k++
  // 2. Creating vector for each list element
  // 3. Updating this vector as the input for the inner loop (replacing seed_vec[i])
  
//  for (int k = 0; k < n_demes; k++) {
//    vector<int> seed_list_subset = seed_list[k];       // When isolating this line for debugging it compiles but an error is returned "Object creates without any names"
//    for (int i = 0; i < param_ptr->seed_infections[k]; i++) {
//      host_pop[host_index[k][i]].denovo_infection(next_haplo_ID,seed_list_subset[i],0);
//    }
//  }
  
//  print_vector(seed_list_subset);
//  Rcpp::stop("debug");
  
  // vector for randomly changing the order in which migration is applied
  mig_order = seq_int(0, parameters.n_mig_list - 1);
  
  // counts of mosquito types
  M_total = sum(param_ptr->M_vec);
  Sv = param_ptr->M_vec;
  Ev = vector<int>(n_demes);
  Iv = vector<int>(n_demes);
  
  // number of mosquitoes at various stages
  n_Ev_death_new = vector<int>(v);
  n_Ev_death = array_2d_int(n_demes, v);
  n_Ev_to_Iv = array_2d_int(n_demes, v);
  
  // population of mosquitoes
  mosq_pop = vector<Mosquito>(M_total);
  
  // initialise mosquitoes
  for (int i = 0; i < M_total; ++i) {
    mosq_pop[i].init(param_ptr);
  }
  
  // store integer index of mosquitoes at various stages
  Sv_index = array_2d_int(n_demes);
  int M_cumul = 0;
  for (int k = 0; k < n_demes; ++k) {
    Sv_index[k] = seq_int(M_cumul, M_cumul + param_ptr->M_vec[k] - 1);
    M_cumul += param_ptr->M_vec[k];
  }
  Ev_death = array_3d_int(n_demes, v);
  Ev_to_Iv = array_3d_int(n_demes, v);
  Iv_index = array_2d_int(n_demes);
  
  // objects for storing results
  daily_values = array_3d_double(n_demes, max_time);
  
  // misc
  EIR = vector<double>(n_demes);
  
}

//------------------------------------------------
// run simulation
void Dispatcher::simulate() {
  
  // start message
  if (param_ptr->report_progress) {
    print("Running simulation");
  }
  
  // initialise ring buffer over extrinsic incubation period
  int ringtime = 0;
  
  // loop through daily time steps
  for (int t = 0; t < max_time; ++t) {
    
    // report progress
    if (param_ptr->report_progress) {
      (*update_progress_ptr)(*args_progress_ptr, "pb", t + 1, max_time);
    }
    
    // skip all dynamics if t = 0. Ensures that user-defined starting conditions
    // are the first values stored in output
    if (t > 0) {
      
      // update ring buffer index
      ringtime = (ringtime == v - 1) ? 0 : ringtime + 1;
      
      
      //-------- MIGRATION --------
      
      reshuffle(mig_order);
      for (int i = 0; i < param_ptr->n_mig_list; ++i) {
        
        // draw number of migrants and skip over if zero
        double prob_migration = get<2>(param_ptr->mig_list[mig_order[i]]);
        if (prob_migration == 0) {
          continue;
        }
        int n_migrants = rbinom1(H, prob_migration);
        if (n_migrants == 0) {
          continue;
        }
        
        // get demes to swap between
        int deme1 = get<0>(param_ptr->mig_list[mig_order[i]]);
        int deme2 = get<1>(param_ptr->mig_list[mig_order[i]]);
        
        // loop through number of migrants
        for (int j = 0; j < n_migrants; ++j) {
          
          // draw migrant index in both demes
          int rnd1 = sample2(0, H - 1);
          int rnd2 = sample2(0, H - 1);
          int index1 = host_index[deme1][rnd1];
          int index2 = host_index[deme2][rnd2];
          
          // update host properties
          host_pop[index1].deme = deme2;
          host_pop[index2].deme = deme1;
          
          // swap host indices
          host_index[deme1][rnd1] = index2;
          host_index[deme2][rnd2] = index1;
          
          // move infectives as needed
          if (host_pop[index1].get_n_active_sexual() > 0) {
            host_infective_index[deme1].erase(index1);
            host_infective_index[deme2].insert(index1);
          }
          if (host_pop[index2].get_n_active_sexual() > 0) {
            host_infective_index[deme2].erase(index2);
            host_infective_index[deme1].insert(index2);
          }
          
        }
      }  // end migration i loop
      
      
      // loop through demes
      for (int k = 0; k < n_demes; ++k) {
        
        
        //-------- MOSQUITO EVENTS --------
        
        
        // carry out Ev death
        int Ev_death_size = Ev_death[k][ringtime].size();
        if (Ev_death_size > 0) {
          push_back_multiple(Sv_index[k], Ev_death[k][ringtime]);
          Ev_death[k][ringtime].clear();
          Ev[k] -= Ev_death_size;
          Sv[k] += Ev_death_size;
        }
        
        // carry out Iv death
        int Iv_death = rbinom1(Iv[k], prob_v_death);
        for (int i = 0; i < Iv_death; ++i) {
          int rnd1 = sample2(0, Iv[k] - 1 - i);
          int this_mosq = Iv_index[k][rnd1];
          mosq_pop[this_mosq].death();
          Sv_index[k].push_back(this_mosq);
          quick_erase(Iv_index[k], rnd1);
        }
        Sv[k] += Iv_death;
        Iv[k] -= Iv_death;
        
        // move Ev into Iv
        int Ev_to_Iv_size = Ev_to_Iv[k][ringtime].size();
        if (Ev_to_Iv_size > 0) {
          push_back_multiple(Iv_index[k], Ev_to_Iv[k][ringtime]);
          Ev_to_Iv[k][ringtime].clear();
          Ev[k] -= Ev_to_Iv_size;
          Iv[k] += Ev_to_Iv_size;
        }
        
        // draw number of new bites on infective hosts. Note, the rate is already
        // scaled by infectivity, meaning bites that do not lead to mosquito
        // infection are ignored automatically
        int n_host_infective = host_infective_index[k].size();                  // number of infective hosts in this deme
        double rate_v_infected = a*infectivity*n_host_infective / double(H);    // rate of mosquito biting infective host and becoming infected
        double prob_v_infected_or_death = 1 - exp(-(rate_v_infected + mu));     // probability of mosquito becoming infected or dying (competing hazards)
        double prob_v_infected = rate_v_infected / (rate_v_infected + mu);      // relative probability of mosquito becoming infected vs. dying
        int v_infected_or_death = rbinom1(Sv[k], prob_v_infected_or_death);     // number of susceptible mosquitoes becoming infected or dying
        int v_infected = rbinom1(v_infected_or_death, prob_v_infected);         // number of susceptible mosquitoes becoming infected
        
        // loop through infective bites
        int death_in_lag = 0;
        for (int i = 0; i < v_infected; ++i) {
          
          // choose mosquito at random from susceptibles
          int rnd1 = sample2(0, Sv[k] - 1);
          int this_mosq = Sv_index[k][rnd1];
          
          // drop from Sv
          quick_erase(Sv_index[k], rnd1);
          
          // the majority of new mosquito infections will die in lag phase. If so
          // then no need to store genotype, instead simply add to Ev_death object
          // to schedule move back to susceptible state (equivalent to death) at
          // future time point.
          int v_time_death = rgeom1(prob_v_death) + 1;
          if (v_time_death <= v) {
            death_in_lag++;
            
            // schedule death
            Ev_death[k][(ringtime + v_time_death) % v].push_back(this_mosq);
            
          } else {
            
            // choose host at random from infectives by using iterator
            int rnd2 = sample2(0, host_infective_index[k].size() - 1);
            set<int>::iterator it = host_infective_index[k].begin();
            advance(it, rnd2);
            int this_host = *it;
            
            // infect mosquito
            mosq_pop[this_mosq].new_infection(&host_pop[this_host]);
            
            // schedule move to Iv
            Ev_to_Iv[k][ringtime].push_back(this_mosq);
          }
          
          // update deme counts
          Sv[k]--;
          Ev[k]++;
          
        } // end loop through infective bites
        
        
        //-------- HUMAN EVENTS --------
        
        
        // get number of new infectious bites on humans
        EIR[k] = a*Iv[k] / double(H);
        double prob_h_infectious_bite = 1 - exp(-EIR[k]);            // probability of new infectious bite per host
        int h_infectious_bite = rbinom1(H, prob_h_infectious_bite);  // total number of new infectious bites
        
        // apply new infectious bites
        for (int i = 0; i < h_infectious_bite; ++i) {
          
          // choose host at random
          int rnd1 = sample2(0, H - 1);
          int this_host = host_index[k][rnd1];
          
          // determine whether infectious bite is successful
          if (rbernoulli1(host_pop[this_host].get_prob_infection())) {
            
            // infect from random mosquito
            int rnd2 = sample2(0, Iv[k] - 1);
            int this_mosq = Iv_index[k][rnd2];
            host_pop[this_host].new_infection(mosq_pop[this_mosq], t);
          }
          
        }  // end loop through infectious bites
        
        // zero host counts for this deme
        Sh[k] = 0;
        Eh[k] = 0;
        Ih[k] = 0;
        
        // loop through all hosts in this deme
        for (int i = 0; i < H; ++i) {
          int this_host = host_index[k][i];
          
          // check for scheduled events on this day
          host_pop[this_host].update_events(next_host_ID, t);
          
          // update host counts
          State_host this_host_state = host_pop[this_host].get_host_state();
          if (this_host_state == Host_Sh) {
            Sh[k]++;
          } else if (this_host_state == Host_Eh) {
            Eh[k]++;
          } else if (this_host_state == Host_Ih) {
            Ih[k]++;
          }
          
        }  // end loop through hosts
        
      } // end loop through demes
      
      
    }  // end of if statement for skipping first timestep
    
    
    //-------- STORE RESULTS --------
    
    // store daily values
    for (int k = 0; k < n_demes; ++k) {
      daily_values[k][t] = {double(Sh[k]), double(Eh[k]), double(Ih[k]),
                            double(Sv[k]), double(Ev[k]), double(Iv[k]),
                            EIR[k]};
    }
    
    // store individual-based output
    int n_sample_demes = param_ptr->sample_list[t].size(); // number of demes to be sampled in this timestep
    if (n_sample_demes > 0) {
      
      for (int i = 0; i < n_sample_demes; ++i) {
        int this_deme = param_ptr->sample_list[t][i].first - 1;
        int this_n = param_ptr->sample_list[t][i].second;
        
        // sample hosts from this deme
        vector<int> rand_vec = sample4(this_n, 0, H - 1);
        for (int j = 0; j < this_n; ++j) {
          int this_index = host_index[this_deme][rand_vec[j]];
          
          Rcpp::List this_sample;
          this_sample["t"] = t + 1;
          this_sample["deme"] = this_deme + 1;
          this_sample["sample_ID"] = host_pop[this_index].ID;
          this_sample["positive"] = (host_pop[this_index].get_host_state() == Host_Ih);
          this_sample["haplotypes"] = host_pop[this_index].get_bloodstage_haplotypes();
          
          sample_output.push_back(this_sample);
        }
      }
    }
    
  } // end time loop
  
}

