
#include "Host.h"
#include "probability_v17.h"
#include "misc_v14.h"

using namespace std;

//------------------------------------------------
// initialise host
void Host::init(int index, int &ID, int deme,
                Parameters &parameters,
                Sampler &sampler_age_stable,
                Sampler &sampler_age_death,
                Sampler &sampler_duration_infection,
                vector<set<int>> &host_infective_index) {
  
  // identifiers
  this->index = index;
  this->ID = ID++;
  home_deme = deme;
  this->deme = deme;
  
  // pointers
  param_ptr = &parameters;
  sampler_age_stable_ptr = &sampler_age_stable;
  sampler_age_death_ptr = &sampler_age_death;
  sampler_duration_infection_ptr = &sampler_duration_infection;
  host_infective_index_ptr = &host_infective_index;
  
  // copy over some parameters for convenience
  L = param_ptr->L;
  max_infections = param_ptr->max_infections;
  n_age = param_ptr->n_age;
  max_time = param_ptr->max_time;
  u = param_ptr->u;
  g = param_ptr->g;
  lambda_products = param_ptr->lambda_products;
  
  // draw birth and death days from stable demography distribution
  draw_starting_age();
  
  // in the unlikely even that die on day zero, delay death by one day
  if (death_day == 0) {
    death_day++;
  }
  
  // initialise cumulative number of infections to date
  cumul_inf = 0;
  
  // initialise infection objects
  infection_active = vector<bool>(max_infections, false);
  infection_status_asexual = vector<Status_asexual>(max_infections, Inactive_asexual);
  infection_status_sexual = vector<Status_sexual>(max_infections, Inactive_sexual);
  
  // initialise times to be greater than maximum simulation time (meaning they
  // will never happen)
  time_Eh_to_Ih = vector<int>(max_infections, max_time + 1);
  time_Ih_to_Sh = vector<int>(max_infections, max_time + 1);
  time_begin_infective = vector<int>(max_infections, max_time + 1);
  time_end_infective = vector<int>(max_infections, max_time + 1);
  
  // initialise time of next event
  time_next_event = death_day;
  
  // initiliase haplotypes
  haplotypes = array_3d_int(max_infections);
  
}

//------------------------------------------------
// draw birth and death days given a known stable demography distribution, and
// schedule death for future. Only needed at start of simulation, after which
// hosts that die are instantly re-born and have death drawn from entire life
// table
void Host::draw_starting_age() {
  
  // draw age from stable demography distribution
  int age_years = sampler_age_stable_ptr->draw();
  int extra_days = sample2(0, 364);
  int age_days = age_years*365 + extra_days;
  
  // draw duration of life from demography distribution looking forward from
  // current age. This is tricky, as we must account for the fact that if we are
  // already part way into an age group then we have a reduced probability of
  // dying within that age group.
  //
  // math derivation:
  // assume constant rate of death r over 1-year period. Total probability of
  // dying this year is therefore p = 1 - exp(-r). We know p from life table,
  // from which we can derive r = -log(1-p). If we are already a proportion x
  // through this year, then the probability of dying in the remaining time is
  // Pr(die) = 1 - exp(-r(1-x)). Sustituting in r and simplifying we get Pr(die)
  // = 1 - (1-p)^(1-x). In code, x can never equal exactly 1 as there must
  // always be at least one day remaining in the current year. Notice that as x
  // approaches 1 the chance of dying this year tends to 0, UNLESS p = 1 in
  // which case death is certain this year.
  int life_days = 0;
  double prop_year_remaining = 1.0 - (extra_days + 1) / 365.0; // (x in the above derivation. Ranges from 0 to 364/365)
  double prob_die_this_year = 1.0 - pow(1.0 - param_ptr->life_table[age_years], 1.0 - prop_year_remaining);
  if (rbernoulli1(prob_die_this_year)) {
    life_days = age_years*365 + sample2(extra_days, 364);
  } else {
    // if we do not die in the current year of age then loop through all
    // remaining years, performing Bernoulli draw from probability of dying in
    // that year
    for (size_t i = (age_years + 1); i < param_ptr->life_table.size(); ++i) {
      if (rbernoulli1(param_ptr->life_table[i])) {
        life_days = i*365 + sample2(0, 364);
        break;
      }
    }
  }
  
  // calculate final birth and death days from age at time 0 and life_days
  birth_day = -age_days;
  death_day = life_days - age_days;
  
}

//------------------------------------------------
// death
void Host::death(int &ID, int t) {
  
  // drop from infective list if necessary
  (*host_infective_index_ptr)[deme].erase(index);
  
  // set new unique ID and increment
  this->ID = ID++;
  
  // make current deme home deme
  home_deme = deme;
  
  // set date of birth as today
  this->birth_day = t;
  
  // draw life duration from demography distribution
  int life_years = sampler_age_death_ptr->draw();
  int life_days = life_years*365 + sample2(1, 364);  // makes simplifying assumption that cannot die on day of birth
  death_day = birth_day + life_days;
  
  // reset cumulative number of infections
  cumul_inf = 0;
  
  // reset infection-level objects
  fill(infection_active.begin(), infection_active.end(), false);
  fill(infection_status_asexual.begin(), infection_status_asexual.end(), Inactive_asexual);
  fill(infection_status_sexual.begin(), infection_status_sexual.end(), Inactive_sexual);
  fill(time_Eh_to_Ih.begin(), time_Eh_to_Ih.end(), max_time + 1);
  fill(time_Ih_to_Sh.begin(), time_Ih_to_Sh.end(), max_time + 1);
  fill(time_begin_infective.begin(), time_begin_infective.end(), max_time + 1);
  fill(time_end_infective.begin(), time_end_infective.end(), max_time + 1);
  
  // reset time of next event
  time_next_event = death_day;
  
  // reset haplotypes
  haplotypes = array_3d_int(max_infections);
  
}

//------------------------------------------------
// de-novo infection
void Host::denovo_infection(int &haplo_ID, int n_haplos, int t) {
  
  // find next infection slot and apply changes to represent new infection
  int this_slot = update_infection_slots(t);
  
  // give each haplotype a unique ID
  haplotypes.arr[this_slot] = vector<vector<int>>(n_haplos);
  for (int i = 0; i < n_haplos; ++i) {
    haplotypes.arr[this_slot][i] = vector<int>(L, haplo_ID);
    haplo_ID++;
  }
  
}

//------------------------------------------------
// new infection from mosquito
void Host::new_infection(Mosquito &mosq, int t) {
  
  // return if already at max_infections
  if (get_n_infections() == max_infections) {
    return;
  }
  
  // find next infection slot and apply changes to represent new infection
  int this_slot = update_infection_slots(t);
  
  // get haplotypes from mosquito
  update_haplotypes_from_mosquito(this_slot, mosq);
  
}

//------------------------------------------------
// update infection slots to represent new infection. Return the chosen slot
int Host::update_infection_slots(int t) {
  
  // update cumulative infections
  cumul_inf++;
  
  // get next free infection slot
  int this_slot = get_free_infection_slot();
  
  // draw duration of infection
  int duration_infection = sampler_duration_infection_ptr->draw() + 1;
  
  // add new infection
  infection_active[this_slot] = true;
  infection_status_asexual[this_slot] = Liverstage_asexual;
  time_Eh_to_Ih[this_slot] = t + u;                               // begin bloodstage
  time_Ih_to_Sh[this_slot] = t + u + duration_infection;          // end bloodstage
  time_begin_infective[this_slot] = t + u + g;                    // begin infective
  time_end_infective[this_slot] = t + u + g + duration_infection; // end infective
  
  // update time of next event
  if ((t + u) < time_next_event) {
    time_next_event = t + u;
  }
  
  return this_slot;
}

//------------------------------------------------
// update haplotypes array by drawing recombinant products from mosquito
void Host::update_haplotypes_from_mosquito(int this_slot, Mosquito &mosq) {
  
  // draw number of recombinant products obtained from mosquito
  int n_products = rztpois1(lambda_products);
  
  // copy over products of recombination
  haplotypes.arr[this_slot] = mosq.get_products(n_products);
  
}

//------------------------------------------------
// enact all events scheduled for this day
void Host::update_events(int &ID, int t) {
  
  // loop through potentially multiple times until all events enacted
  while (time_next_event == t) {
    
    // check for death
    if (death_day == t) {
      death(ID, t);
    }
    
    // check for events on individual infection objects
    for (int i = 0; i < max_infections; ++i) {
      if (time_Eh_to_Ih[i] == t) {
        Eh_to_Ih(i);
      }
      if (time_Ih_to_Sh[i] == t) {
        Ih_to_Sh(i);
      }
      if (time_begin_infective[i] == t) {
        begin_infective(i);
      }
      if (time_end_infective[i] == t) {
        end_infective(i);
      }
    }
    
    // update time to next event
    update_time_next_event();
    
  }
  
}

//------------------------------------------------
// update time to next events
void Host::update_time_next_event() {
  time_next_event = death_day;
  for (int i = 0; i < max_infections; ++i) {
    time_next_event = (time_Eh_to_Ih[i] < time_next_event) ? time_Eh_to_Ih[i] : time_next_event;
    time_next_event = (time_Ih_to_Sh[i] < time_next_event) ? time_Ih_to_Sh[i] : time_next_event;
    time_next_event = (time_begin_infective[i] < time_next_event) ? time_begin_infective[i] : time_next_event;
    time_next_event = (time_end_infective[i] < time_next_event) ? time_end_infective[i] : time_next_event;
  }
}

//------------------------------------------------
// move from Eh state to Ih
void Host::Eh_to_Ih(int this_slot) {
  
  // update infection status
  infection_status_asexual[this_slot] = Bloodstage_asexual;
  
  // reset time to next event
  time_Eh_to_Ih[this_slot] = max_time + 1;
  
}

//------------------------------------------------
// move from Ih state to Sh
void Host::Ih_to_Sh(int this_slot) {
  
  // update infection status
  infection_status_asexual[this_slot] = Inactive_asexual;
  
  // reset time to next event
  time_Ih_to_Sh[this_slot] = max_time + 1;
  
}

//------------------------------------------------
// begin infective period
void Host::begin_infective(int this_slot) {
  
  // update infection status
  infection_status_sexual[this_slot] = Active_sexual;
  
  // reset time to next event
  time_begin_infective[this_slot] = max_time + 1;
  
  // add to infectives set (will do nothing if already in set)
  (*host_infective_index_ptr)[deme].insert(index);
  
}

//------------------------------------------------
// end infective period
void Host::end_infective(int this_slot) {
  
  // update infection status
  infection_status_sexual[this_slot] = Inactive_sexual;
  infection_active[this_slot] = false;
  
  // reset time to next event
  time_end_infective[this_slot] = max_time + 1;
  
  // if no longer infective then drop from infectives set
  if (get_n_active_sexual() == 0) {
    (*host_infective_index_ptr)[deme].erase(index);
  }
  
  // clear heplotypes
  haplotypes[this_slot].clear();
  
}

//------------------------------------------------
// sample from active slots with equal probability
int Host::draw_active_slot() {
  
  // error if no active infections to sample
  int n_active_sexual = get_n_active_sexual();
  if (n_active_sexual == 0) {
    Rcpp::stop("no active infections to sample");
  }
  
  // sample from active slots only
  int ret = 0;
  int n = sample2(1, n_active_sexual);
  for (int j = 0; j < max_infections; ++j) {
    n -= infection_active[j];
    if (n == 0) {
      ret = j;
      break;
    }
  }
  
  return ret;
}

// #################################
// #                               #
// #      GETTERS AND SETTERS      #
// #                               #
// #################################

//------------------------------------------------
// get total number of infections. This counts infections in either asexual or
// sexual states
int Host::get_free_infection_slot() {
  
  int ret = 0;
  for (int i = 0; i < max_infections; ++i) {
    if (!infection_active[i]) {
      break;
    }
    ret++;
  }
  if (ret == max_infections) {
    Rcpp::stop("could not find free infection slot");
  }
  
  return ret;
}

//------------------------------------------------
// get total number of infections. This counts infections in either asexual or
// sexual states
int Host::get_n_infections() {
  return sum_bool(infection_active);
}

//------------------------------------------------
// get number of infections that have asexual liverstage parasites present
int Host::get_n_liverstage_asexual() {
  int ret = 0;
  for (int i = 0; i < max_infections; ++i) {
    if (infection_status_asexual[i] == Liverstage_asexual) {
      ret++;
    }
  }
  return ret;
}

//------------------------------------------------
// get number of infections that have asexual bloodstage parasites present
int Host::get_n_bloodstage_asexual() {
  int ret = 0;
  for (int i = 0; i < max_infections; ++i) {
    if (infection_status_asexual[i] == Bloodstage_asexual) {
      ret++;
    }
  }
  return ret;
}

//------------------------------------------------
// get number of infections for which gametocytes are present
int Host::get_n_active_sexual() {
  int ret = 0;
  for (int i = 0; i < max_infections; ++i) {
    if (infection_status_sexual[i] == Active_sexual) {
      ret++;
    }
  }
  return ret;
}

//------------------------------------------------
// get host state
State_host Host::get_host_state() {
  if (get_n_bloodstage_asexual() == 0) {
    if (get_n_liverstage_asexual() == 0) {
      return Host_Sh;
    } else {
      return Host_Eh;
    }
  }
  return Host_Ih;
}

//------------------------------------------------
// get current probability of infection
double Host::get_prob_infection() {
  int p = int(param_ptr->prob_infection.size());
  int i = (cumul_inf > (p - 1)) ? p - 1 : cumul_inf;
  return param_ptr->prob_infection[i];
}

//------------------------------------------------
// return array of all haplotypes that correspond to bloodstage infections
vector<vector<int>> Host::get_bloodstage_haplotypes() {
  vector<vector<int>> ret;
  for (int i = 0; i < max_infections; ++i) {
    if (infection_status_asexual[i] == Bloodstage_asexual) {
      for (int j = 0; j < haplotypes[i].size(); ++j) {
        ret.push_back(haplotypes[i][j]);
      }
    }
  }
  return ret;
}
