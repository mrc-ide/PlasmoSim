
#------------------------------------------------
#' @title Simulate genetic data from simple P. falciparum model
#'
#' @description Simulate genetic data from a simple model of P. falciparum
#'   epidemiology and genetics.
#'
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on
#'   humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. mu = -log(p) unless otherwise
#'   specified.
#' @param u intrinsic incubation period. The number of days from infection to
#'   blood-stage infection in a human host.
#' @param v extrinsic incubation period. The number of days from infection to
#'   becoming infectious in a mosquito.
#' @param g lag time between human blood-stage infection and production of
#'   gametocytes.
#' @param prob_infection probability a human becomes infected after being bitten
#'   by an infected mosquito.
#' @param duration_infection vector specifying probability distribution of time
#'   (in days) of a malaria episode.
#' @param infectivity probability a mosquito becomes infected after biting an
#'   infective human host.
#' @param max_infections maximum number of infections that an individual
#'   can hold simultaneously.
#' @param H human population size, which is assumed to be the same in every
#'   deme.
#' @param seed_infections vector specifying the initial number of infected
#'   humans in each deme.
#' @param M vector specifying mosquito population size (strictly the number of
#'   adult female mosquitoes) in each deme.
#' @param mig_matrix migration matrix specifing the daily probability of
#'   migrating from each deme to each other deme. Migration must be equal in
#'   both directions, meaning this matrix must be symmetric.
#' @param L number of loci. The maximum number of loci is 1000, as at higher
#'   numbers haplotypes begin to exceed integer representation (2^L).
#' @param mean_oocysts the average number of viable oocysts generated from
#'   gametocytes upon biting an infective host. The actual number of oocysts is
#'   generated from a zero-truncated Poisson distribution with this mean.
#' @param mean_products parasite genotypes are passed from mosquito to host by
#'   sampling N times with replacement from the available oocysts products (the
#'   available number of products is 4 times the number of oocysts). N is drawn
#'   independently for each infection from a zero-truncated Poisson distribution
#'   with mean given by \code{mean_products}. Hence, large values of this
#'   parameter increase the chance of co-transmission of multiple genotypes,
#'   while small values increase the chance of picking up just a single
#'   genotype.
#' @param recomb_prob the probability of a recombination breakpoint between any
#'   sequential pair of loci. Assumed to be the same for all loci.
#' @param max_time number of days in the simulation.
#' @param sample_dataframe a dataframe specifying outputs from the model. Must
#'   contain the following three columns:
#'   \enumerate{
#'     \item \code{deme}: which numbered deme to sample from.
#'     \item \code{time}: the day on which samples are taken.
#'     \item \code{n}: the number of hosts to randomly sample from the population.
#'   }
#' @param report_progress if \code{TRUE} then a progress bar is printed to the
#'   console during simuation.
#'
#' @importFrom utils txtProgressBar
#' @importFrom stats dgeom setNames optim
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

sim_falciparum <- function(a = 0.3,
                           p = 0.9,
                           mu = -log(p),
                           u = 12,
                           v = 10,
                           g = 10,
                           prob_infection = 0.1,
                           duration_infection = dgeom(1:500, 1/100),
                           infectivity = 0.1,
                           max_infections = 5,
                           H = 1000,
                           seed_infections = 100,
                           seed_vec = sample.int(5, seed_infections, replace = TRUE),
                           seed_list = lapply(seed_infections,sample.int,n=5,replace=TRUE),
                           M = 1000,
                           mig_matrix = diag(length(M)),
                           L = 24,
                           mean_oocysts = 2.0,
                           mean_products = 5,
                           recomb_prob = 0.1,
                           max_time = max(sample_dataframe$time),
                           sample_dataframe = data.frame(deme = 1, time = 365, n = 100),
                           report_progress = TRUE) {
  
  # avoid no visible binding warning
  S <- NULL

  # check inputs
  assert_single_bounded(a)
  assert_single_bounded(p)
  assert_single_pos(mu)
  assert_single_pos_int(u, zero_allowed = FALSE)
  assert_single_pos_int(v, zero_allowed = FALSE)
  assert_single_pos_int(g, zero_allowed = FALSE)
  assert_vector_bounded(prob_infection)
  assert_pos(prob_infection[1], zero_allowed = FALSE)
  assert_vector_pos(duration_infection, zero_allowed = TRUE)
  assert_single_bounded(infectivity)
  assert_single_pos_int(max_infections, zero_allowed = FALSE)
  assert_single_pos_int(H, zero_allowed = FALSE)
  assert_pos_int(seed_infections, zero_allowed = TRUE)
  assert_leq(seed_infections, H)
  assert_vector_pos(seed_vec, zero_allowed = FALSE)
  assert_list(seed_list)
  assert_eq(lengths(seed_list),seed_infections)
  assert_pos(unlist(seed_list), zero_allowed = TRUE)
  assert_same_length(M, seed_list)
  assert_pos_int(M, zero_allowed = FALSE)
  assert_same_length(M, seed_infections)
  n_demes <- length(M)
  assert_symmetric_matrix(mig_matrix)
  assert_dim(mig_matrix, c(n_demes, n_demes))
  assert_bounded(mig_matrix)
  assert_eq(rowSums(mig_matrix), rep(1,n_demes))
  assert_single_pos_int(L, zero_allowed = FALSE)
  assert_leq(L, 1000)
  assert_single_pos(mean_oocysts)
  assert_gr(mean_oocysts, 1.0)
  assert_single_pos(mean_products)
  assert_gr(mean_products, 1.0)
  assert_single_bounded(recomb_prob)
  assert_single_pos_int(max_time, zero_allowed = FALSE)
  assert_dataframe(sample_dataframe)
  assert_in(c("deme", "time", "n"), names(sample_dataframe))
  assert_pos_int(sample_dataframe$deme, zero_allowed = FALSE)
  assert_leq(sample_dataframe$deme, n_demes, message = sprintf("sample_dataframe demes must be consistent with the number of demes implied by the migration matrix, i.e. between 1 and %s", n_demes))
  assert_pos_int(sample_dataframe$time, zero_allowed = FALSE)
  assert_pos_int(sample_dataframe$n, zero_allowed = FALSE)
  assert_single_logical(report_progress)
  
  # normalise infection duration distribution
  duration_infection <- duration_infection / sum(duration_infection)
  
  # if any prob_infection is zero then this automatically defines the value of
  # max_infections
  if (any(prob_infection == 0)) {
    max_infections <- min(max_infections, which(prob_infection == 0)[1] - 1)
  }
  
  # solve zero-truncated Poisson distribution for lambda_oocysts parameter to
  # match user-specified mean
  lambda_loss <- function(lambda, mu) {
    abs(lambda / (1 - exp(-lambda)) - mu)
  }
  lambda_oocysts <- optim(1.0, lambda_loss, method = "Brent", lower = 0.0, upper = 11.0, mu = mean_oocysts)$par
  assert_non_NA(lambda_oocysts)
  assert_non_null(lambda_oocysts)
  assert_le(lambda_oocysts, 10, message = "mean_oocysts too high")
  
  # solve zero-truncated Poisson distribution for lambda_products parameter to
  # match user-specified mean
  lambda_products <- optim(1.0, lambda_loss, method = "Brent", lower = 0.0, upper = 101.0, mu = mean_products)$par
  assert_non_NA(lambda_products)
  assert_non_null(lambda_products)
  assert_le(lambda_products, 100, message = "mean_products too high")
  
  # read in Mali demography distribution
  mali_demog <- plasmosim_file("mali_demog.rds")
  
  
  # ---------------------------------------------
  # set up arguments for input into C++
  
  # create function list
  args_functions <- list(update_progress = update_progress)
  
  # create progress bars
  pb <- txtProgressBar(0, max_time, initial = NA, style = 3)
  args_progress <- list(pb = pb)
  
  # create argument list
  args <- list(a = a,
               mu = mu,
               u = u,
               v = v,
               g = g,
               prob_infection = prob_infection,
               duration_infection = duration_infection,
               infectivity = infectivity,
               max_infections = max_infections,
               H = H,
               seed_infections = seed_infections,
               seed_vec = seed_vec,
               seed_list = seed_list,
               M = M,
               mig_matrix = matrix_to_rcpp(mig_matrix),
               L = L,
               lambda_oocysts = lambda_oocysts,
               lambda_products = lambda_products,
               recomb_prob = recomb_prob,
               life_table = mali_demog$life_table,
               age_death = mali_demog$age_death,
               age_stable = mali_demog$age_stable,
               max_time = max_time,
               sample_dataframe = sample_dataframe,
               report_progress = report_progress)
  
  # ---------------------------------------------
  # run efficient C++ function
  
  output_raw <- sim_falciparum_cpp(args, args_functions, args_progress)
  
  #return(output_raw)
  # ---------------------------------------------
  # process raw output
  
  message("processing output")
  
  # get daily values
  daily_values <- mapply(function(i) {
    output_raw$daily_values[[i]] %>%
      rcpp_to_matrix() %>%
      as.data.frame() %>%
      setNames(c("S", "E", "I", "Sv", "Ev", "Iv", "EIR")) %>%
      dplyr::mutate(time = seq_along(S),
                    deme = i,
                    .before = 1)
  }, seq_along(output_raw$daily_values), SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()
  
  # process basic individual-level data
  indlevel <- mapply(function(x) {
    data.frame(time = x$t,
               deme = x$deme,
               sample_ID = x$sample_ID,
               positive = x$positive)
  }, output_raw$sample_output, SIMPLIFY = FALSE) %>%
    dplyr::bind_rows()
  
  # append haplotypes as list
  indlevel$haplotypes <- mapply(function(x) {
    l <- length(x$haplotypes)
    if (l > 0) {
      ret <- matrix(unlist(x$haplotypes), nrow = l, byrow = TRUE)
    } else {
      ret <- NULL
    }
    return(ret)
  }, output_raw$sample_output)
  
  # get unique hash of each haplotype
  indlevel$haplo_ID <- mapply(function(x) {
    if (is.null(x)) {
      ret <- NULL
    } else {
      haplo_string <- paste(x, collapse = ".")
      ret <- openssl::md5(haplo_string)
    }
    return(ret)
  }, indlevel$haplotypes)
  
  # reorder columns
  indlevel <- indlevel %>% 
    dplyr::relocate(.data$haplo_ID, .before = .data$haplotypes)
  
  # return list
  ret <- list(daily_values = daily_values,
              indlevel = indlevel,
              final_haplos = output_raw$final_haplos)
  
  return(ret)
}

#------------------------------------------------
#' @title Get proportion identical between two haplotype matrices
#'
#' @description Compare two sets of haplotypes (matrices), and return the
#'   proportion of identical sites over all pairwise comparisons. Values can be
#'   any numeric value; for example if values represent ancestry then this
#'   function returns the average identity by descent, or if values represent
#'   alleles then it returns the average identity by state.
#'
#' @param mat1,mat2 matrices representing sets of haplotypes to compare.
#'   Haplotypes are in rows and loci are in columns.
#'   
#' @export

get_haplotype_identity <- function(mat1, mat2) {
  
  # check inputs
  assert_matrix_numeric(mat1)
  assert_matrix_numeric(mat2)
  assert_eq(ncol(mat1), ncol(mat2), message = "haplotype matrices must have the same number of columns (loci)")
  
  # get number of matches and number of comparisons between all pairs of
  # haplotypes
  numer <- sum(apply(mat1, 1, function(x) {
    sum(apply(mat2, 1, function(y) x == y))
  }))
  denom <- ncol(mat1)*nrow(mat1)*nrow(mat2)
  
  # return proportion identical
  return(numer / denom)
}

#------------------------------------------------
#' @title Get pairwise genetic identity matrix
#'
#' @description Calculates pairwise genetic identity between all samples. If
#'   \code{deme_level = TRUE} this is averaged over all individuals within a
#'   deme, to produce average pairwise relatedness within and between demes.
#'   Each time point in the sample is considered independently, and output as a
#'   list.
#'
#' @param sim_output simulation output from \code{sim_falciparum()}.
#' @param deme_level if \code{TRUE} then return pairwise identity at the deme
#'   level, averaged over all individuals within each deme. Otherwise return at
#'   the individual level.
#'   
#' @export

get_identity_matrix <- function(sim_output, deme_level = FALSE) {
  
  # avoid no visible binding note
  positive <- deme <- NULL
  
  # check inputs
  assert_list_named(sim_output)
  assert_in("indlevel", names(sim_output))
  assert_dataframe(sim_output$indlevel)
  assert_in(c("time", "deme", "sample_ID", "positive", "haplotypes"), names(sim_output$indlevel))
  assert_vector_pos_int(sim_output$indlevel$time)
  assert_vector_pos_int(sim_output$indlevel$deme)
  assert_vector_pos_int(sim_output$indlevel$sample_ID)
  assert_vector_logical(sim_output$indlevel$positive)
  assert_list(sim_output$indlevel$haplotypes)
  
  # split individual-level output by sampling time
  tsplit <- split(sim_output$indlevel, f = sim_output$indlevel$time)
  
  # loop over sampling times
  indlevel_list <- mapply(function(x) {
    
    # subset to positives
    x <- subset(x, positive == 1)
    if (nrow(x) == 0) {
      return(NULL)
    }
    
    # get identity in all pairwise comparisons
    ret_mat <- mapply(function(i) {
      mapply(function(j) {
        get_haplotype_identity(x$haplotypes[[i]], x$haplotypes[[j]])
      }, seq_along(x$haplotypes))
    }, seq_along(x$haplotypes))
    
    # set row and column names
    colnames(ret_mat) <-  rownames(ret_mat) <- x$sample_ID
    
    ret_mat
  }, tsplit, SIMPLIFY = FALSE)
  
  # if not summarising at deme level then return as is
  if (!deme_level) {
    return(indlevel_list)
  }
  
  # produce pairwise matrix summarised at deme level
  demelevel_list <- mapply(function(t_i) {
    
    # subset to positives
    x <- subset(tsplit[[t_i]], positive == 1)
    if (nrow(x) == 0) {
      return(NULL)
    }
    
    # get dataframe of demes and IDs
    deme_IDs <- subset(x, select = c("deme", "sample_ID"))
    
    # get unique demes
    u <- unique(x$deme)
    
    # get pairwise similarity averaged over all individuals within a deme
    ret_mat <- mapply(function(i) {
      mapply(function(j) {
        IDs_i <- subset(deme_IDs, deme == u[i])$sample_ID
        IDs_j <- subset(deme_IDs, deme == u[j])$sample_ID
        w1 <- match(IDs_i, as.numeric(colnames(indlevel_list[[t_i]])))
        w2 <- match(IDs_j, as.numeric(colnames(indlevel_list[[t_i]])))
        mean(indlevel_list[[t_i]][w1, w2])
      }, seq_along(u))
    }, seq_along(u))
    
    # set row and column names
    colnames(ret_mat) <-  rownames(ret_mat) <- u
    
    ret_mat
  }, seq_along(tsplit), SIMPLIFY = FALSE)
  
  return(demelevel_list)
}
