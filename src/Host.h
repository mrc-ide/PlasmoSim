
#pragma once

#include "Parameters.h"
#include "Mosquito.h"
#include "array.h"
#include "Sampler_v4.h"

#include <vector>
#include <map>

//------------------------------------------------
// enumerate possible asexual and sexual infection status
enum Status_asexual {Inactive_asexual, Liverstage_asexual, Bloodstage_asexual};
enum Status_sexual {Inactive_sexual, Active_sexual};

//------------------------------------------------
// enumerate possible host states
enum State_host {Host_Sh, Host_Eh, Host_Ih};

//------------------------------------------------
// class defining host
class Host {
  
public:
  
  // PUBLIC OBJECTS
  
  // identifiers
  int index;      // where in the population (vector of hosts) this host resides
  int ID;         // unique ID, incremented upon death
  int home_deme;  // deme into which this host was born
  int deme;       // deme in which this host currently resides
  
  // pointers to external objects
  Parameters* param_ptr;
  Sampler* sampler_age_stable_ptr;
  Sampler* sampler_age_death_ptr;
  Sampler* sampler_duration_infection_ptr;
  std::vector<std::set<int>>* host_infective_index_ptr;
  
  // copy over some parameters for convenience
  int L;
  int max_infections;
  int n_age;
  int max_time;
  int u;
  int g;
  double lambda_products;
  
  // dates of birth and death
  int birth_day;
  int death_day;
  
  // cumulative number of infections to date
  int cumul_inf;
  
  // infection-level objects
  std::vector<bool> infection_active;
  std::vector<Status_asexual> infection_status_asexual;
  std::vector<Status_sexual> infection_status_sexual;
  std::vector<int> time_Eh_to_Ih;
  std::vector<int> time_Ih_to_Sh;
  std::vector<int> time_begin_infective;
  std::vector<int> time_end_infective;
  
  // counts of the number of infections in each stage
  int n_liverstage_asexual;
  int n_bloodstage_asexual;
  int n_active_sexual;
  
  // time of next event
  int time_next_event;
  
  // haplotypes
  array_3d_int haplotypes;  // each infection slot can hold multiple haplotypes. Each haplotype is a vector of integers
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  
  // other methods
  void init(int index, int &ID, int deme,
            Parameters &parameters,
            Sampler &sampler_age_stable,
            Sampler &sampler_age_death,
            Sampler &sampler_duration_infection,
            std::vector<std::set<int>> &host_infective_index);
  void draw_starting_age();
  void death(int &ID, int t);
  
  void denovo_infection(int &haplo_ID);
  void new_infection(Mosquito &mosq, int t);
  void update_events(int &ID, int t);
  void update_time_next_event();
  
  void Eh_to_Ih(int this_slot);
  void Ih_to_Sh(int this_slot);
  void begin_infective(int this_slot);
  void end_infective(int this_slot);
  
  // getters and setters
  int get_n_infections();
  State_host get_host_state();
  double get_prob_infection();
  std::vector<std::vector<int>> get_bloodstage_haplotypes();
  
};
