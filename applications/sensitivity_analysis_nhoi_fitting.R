library(dplyr)
source("helper_functions.R")

# This file contains code used to perform a sensitivity analysis for the 
# no-highest-order interaction assumption for the Kosovo data set

#### Load posterior samples ####
load("nhoi_fits/pi_tilde_samples.RData")
load("nhoi_fits/N_samples.RData")

#### Number of files and number of observed casualties ####
K <- 4
n <- 4400

#### Hyperparameters for negative-binomial priors for N ####
M <- c(8000, 12000, 10000)
a <- c(43, 9, 1.6)

#### Range of values for sensitivity parameter ####
xi <- c(1 / 2, 2 / 3, 3 / 2, 2)

#### Conduct sensitivity analysis using LCMCR fit ####
N_lcmcr$interaction <- 1
set.seed(61)
for(j in 1:length(xi)){
    pi_0_temp <- get_pi_0_hoi(pi_tilde_lcmcr, K, xi = xi[j])
    for(i in 1:length(M)){
        temp_accept <- rejection_sampler_nbinom(pi_0_temp, M[i], a[i], n)
        temp_N <- n + n_0_sampler_nbinom(pi_0_temp[temp_accept], M[i], a[i], n)
        N_lcmcr <- rbind(N_lcmcr, data.frame(N = temp_N, prior = i, 
                                             interaction = xi[j]))
    }
    N_lcmcr <- rbind(N_lcmcr, 
                     data.frame(N = n + 
                                    n_0_sampler_isp(pi_0_temp[pi_0_temp > 0], 
                                                    n),
                                prior = length(M) + 1, interaction = xi[j]))
}

#### Save the relevant output ####
save(N_lcmcr, file = "nhoi_fits/N_samples_sensitivity.RData")
