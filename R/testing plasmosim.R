library(dplyr)
library(PlasmoSim)
sim1 <- sim_falciparum(sample_dataframe = data.frame(deme = 1, time = 13, n = 1000))
sim1$indlevel$haplotypes[[1]]

# MIGRATION
# define x- and y-locations of demes
demes_x <- seq(0, 100, 10)
demes_y <- seq(0, 100, 10)

# get coordinates into dataframe
deme_df <- expand.grid(demes_x, demes_y)
names(deme_df) <- c("x", "y")

# add deme-specific properties
deme_df <- deme_df %>%
  dplyr::mutate(deme = seq_len(nrow(deme_df)),
                M = 500 + x, # mosquito population size increasing from left to right
                seed_infections = ifelse(x == 0 & y == 0, 50, 0))

####### NEW START ####### 
seed_vec <- list()
for (i in 1:length(deme_df$seed_infections)) {
  seed_vec[[i]] = sample.int(5, deme_df$seed_infections[i], replace = TRUE)
}
deme_df$seed_vec <- seed_vec
####### NEW END ####### 

head(deme_df)

# define migration matrix based on distance
d <- deme_df %>%
  dplyr::select(x, y) %>%
  dist() %>%
  as.matrix()
mig_matrix <- 0.01 * (d <= 10)

# ensure migration probabilities sum to 1 over rows
diag(mig_matrix) <- 0
diag(mig_matrix) <- 1 - rowSums(mig_matrix)

# define output dataframe
sample_df <- deme_df %>%
  dplyr::filter(deme %in% c(13, 26, 35, 109)) %>%
  dplyr::mutate(time = 30 * 365,
                n = 100)

sample_df

# simulate
set.seed(1)
sim_mig <- sim_falciparum(H = 100,
                          M = deme_df$M,
                          seed_infections = deme_df$seed_infections,
                          mig_matrix = mig_matrix,
                          L = 24,
                          sample_dataframe = sample_df)
