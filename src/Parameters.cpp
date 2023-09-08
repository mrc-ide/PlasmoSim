
#include "Parameters.h"
#include "misc_v14.h"
#include "array_v1.h"

using namespace std;

//------------------------------------------------
// constructor
Parameters::Parameters(const Rcpp::List &args) {
  
  // genetic parameters
  L = rcpp_to_int(args["L"]);
  lambda_oocysts = rcpp_to_double(args["lambda_oocysts"]);
  lambda_products = rcpp_to_double(args["lambda_products"]);
  recomb_prob = rcpp_to_double(args["recomb_prob"]);
  
  // epidemiological parameters
  a = rcpp_to_double(args["a"]);
  mu = rcpp_to_double(args["mu"]);
  u = rcpp_to_int(args["u"]);
  v = rcpp_to_int(args["v"]);
  g = rcpp_to_int(args["g"]);
  prob_infection = rcpp_to_vector_double(args["prob_infection"]);
  duration_infection = rcpp_to_vector_double(args["duration_infection"]);
  infectivity = rcpp_to_double(args["infectivity"]);
  max_infections = rcpp_to_int(args["max_infections"]);
  
  // deme parameters
  H = rcpp_to_int(args["H"]);
  seed_infections = rcpp_to_vector_int(args["seed_infections"]);
  M_vec = rcpp_to_vector_int(args["M"]);
  n_demes = int(M_vec.size());
  
  // get migration list from matrix
  vector<vector<double>> mig_matrix = rcpp_to_matrix_double(args["mig_matrix"]);
  for (int i = 0; i < (n_demes - 1); ++i) {
    for (int j = (i + 1); j < n_demes; ++j) {
      if (mig_matrix[i][j] > 0) {
        mig_list.push_back(make_tuple(i, j, mig_matrix[i][j]));
      }
    }
  }
  n_mig_list = mig_list.size();
  
  // demography
  life_table = rcpp_to_vector_double(args["life_table"]);
  age_death = rcpp_to_vector_double(args["age_death"]);
  age_stable = rcpp_to_vector_double(args["age_stable"]);
  n_age = int(age_stable.size());
  
  // run parameters
  max_time = rcpp_to_int(args["max_time"]);
  report_progress = rcpp_to_bool(args["report_progress"]);
  
  // sampling parameters. sample_list contains an element for each day of
  // simulation. Each element is a vector over all demes that need to be sampled
  // on this day. Each element of the vector is a pair that gives the index of
  // the deme and the number of hosts to sample.
  Rcpp::List sample_dataframe = args["sample_dataframe"];
  vector<int> sample_df_deme = rcpp_to_vector_int(sample_dataframe["deme"]);
  vector<int> sample_df_time = rcpp_to_vector_int(sample_dataframe["time"]);
  vector<int> sample_df_n = rcpp_to_vector_int(sample_dataframe["n"]);
  sample_list = vector<vector<pair<int, int>>>(max_time);
  for (size_t i = 0; i < sample_df_deme.size(); ++i) {
    int t = sample_df_time[i] - 1;
    if (t < max_time) {
      sample_list[t].push_back({sample_df_deme[i], sample_df_n[i]});
    }
  }
  
  // misc parameters
  prob_v_death = 1 - exp(-mu);  // daily probability of mosquito death
  
}

//------------------------------------------------
// print summary
void Parameters::print_summary() {
  
  // epidemiological parameters
  print("-- epidemiological parameters --");
  print("a: ", a);
  print("mu: ", mu);
  print("u: ", u);
  print("v: ", v);
  print("g: ", g);
  print("prob_infection: ");
  print_vector(prob_infection);
  print("duration_infection: ");
  print_vector(duration_infection);
  print("infectivity: ", infectivity);
  print("max_infections: ", max_infections);
  print("");
  
  // deme parameters
  print("-- deme parameters --");
  print("H: ", H);
  print("seed_infections: ");
  print_vector(seed_infections);
  print("M_vec: ");
  print_vector(M_vec);
  print("n_demes: ", n_demes);
  print("");
  
  // genetic parameters
  print("-- genetic parameters --");
  print("L: ", L);
  print("lambda_oocysts: ", lambda_oocysts);
  print("lambda_products: ", lambda_products);
  
  // demography
  print("-- demography --");
  print("life_table: ");
  print_vector(life_table);
  print("age_death: ");
  print_vector(age_death);
  print("age_stable: ");
  print_vector(age_stable);
  print("n_age: ", n_age);
  print("");
  
  // migration
  print("-- migration --");
  for (int i = 0; i < n_mig_list; ++i) {
    print(get<0>(mig_list[i]), get<1>(mig_list[i]), get<2>(mig_list[i]));
  }
  print("");
  
  // run parameters
  print("-- run parameters --");
  print("max_time: ", max_time);
  print("");
  
  // misc parameters
  print("-- misc parameters --");
  print("prob_v_death: ", prob_v_death);
  print("");
}
