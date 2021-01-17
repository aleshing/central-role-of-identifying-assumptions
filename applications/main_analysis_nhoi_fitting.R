library(MCMCpack)
library(cmdstanr)
library(LCMCR)
library(conting)
library(posterior)
library(dplyr)
source("helper_functions.R")

# This file contains code used to analyze the Kosovo data set using the 
# no-highest-order interaction assumption

#### Read in Kosovo data from the LCMCR package ####
data(kosovo_aggregate)
K <- ncol(kosovo_aggregate)
n <- nrow(kosovo_aggregate)
kosovo_cont <- kosovo_aggregate %>% group_by(EXH, ABA, OSCE, HRW) %>% 
    summarise(Freq = n()) %>% ungroup() %>% as.data.frame()

#### Hyperparameters for negative-binomial priors for N ####
M <- c(8000, 12000, 10000)
a <- c(43, 9, 1.6)

#### Dirichlet prior for observed cell probabilities ####
set.seed(60)
num_dir <- 100000
alphas <- rep(1, nrow(kosovo_cont))
pi_tilde_dir <- rdirichlet(num_dir, alpha = kosovo_cont$Freq + alphas)
pi_0_dir <- get_pi_0_hoi(pi_tilde_dir, K, xi = 1)
set.seed(61)
N_dir <- data.frame(N = numeric(), prior = numeric())
for(i in 1:length(M)){
    temp_accept <-  rejection_sampler_nbinom(pi_0_dir, M[i], a[i], n)
    temp_N <- n + n_0_sampler_nbinom(pi_0_dir[temp_accept], M[i], a[i], n)
    N_dir <- rbind(N_dir, data.frame(N = temp_N, prior = i))
}
N_dir <- rbind(N_dir, data.frame(N = n + n_0_sampler_isp(pi_0_dir, n),
                                 prior = length(M) + 1))

#### Log-linear prior for observed cell probabilities, using Stan ####
set.seed(60)
log_linear_model <- cmdstan_model('log_linear_model.stan')
X <- glm(Freq ~ (EXH + ABA + OSCE + HRW) ^ 3, data = kosovo_cont,
         family = poisson, x = TRUE)$x[, -1]
ll_sd <- 5
stan_data <- list(K = K, twoK = 2 ^ K, ns = kosovo_cont$Freq, n = n, 
                  p = ncol(X), X = X, ll_sd = ll_sd)
ll_fit <- log_linear_model$sample(data = stan_data, seed = 60, 
                                  num_chains = 4, num_cores = 2, 
                                  max_depth = 10, adapt_delta = 0.8,
                                  num_warmup = 1000, num_samples = 25000)
pi_tilde_ll <- ll_fit$draws() %>% 
    subset_draws(variable = "pi_tilde", regex = TRUE) %>% 
    merge_chains() %>% as_draws_df() %>% select(-c(.chain, .iteration, .draw)) 
pi_0_ll <- get_pi_0_hoi(pi_tilde_ll, K, xi = 1)
set.seed(61)
N_ll <- data.frame(N = numeric(), prior = numeric())
for(i in 1:length(M)){
    temp_accept <-  rejection_sampler_nbinom(pi_0_ll, M[i], a[i], n)
    temp_N <- n + n_0_sampler_nbinom(pi_0_ll[temp_accept], M[i], a[i], n)
    N_ll <- rbind(N_ll, data.frame(N = temp_N, prior = i))
}
N_ll <- rbind(N_ll, data.frame(N = n + n_0_sampler_isp(pi_0_ll, n),
                                 prior = length(M) + 1))

#### Latent class model prior for observed cell probabilities, using LCMCR ####
set.seed(60)
num_lcmcr_burn <- 50000
num_lcmcr <- 100000
lcmcr_fit <- fit_lcmcr(kosovo_cont, K, burn_iter = num_lcmcr_burn, 
                       samp_iter = num_lcmcr, K_star = 10, seed = 60)
pi_tilde_lcmcr <- get_pi_tilde_lcmcr(lcmcr_fit, K)
pi_0_lcmcr <- get_pi_0_hoi(pi_tilde_lcmcr, K, xi = 1)
set.seed(61)
N_lcmcr <- data.frame(N = numeric(), prior = numeric())
for(i in 1:length(M)){
    temp_accept <-  rejection_sampler_nbinom(pi_0_lcmcr, M[i], a[i], n)
    temp_N <- n + n_0_sampler_nbinom(pi_0_lcmcr[temp_accept], M[i], a[i], n)
    N_lcmcr <- rbind(N_lcmcr, data.frame(N = temp_N, prior = i))
}
N_lcmcr <- rbind(N_lcmcr, data.frame(N = n + n_0_sampler_isp(pi_0_lcmcr, n),
                                     prior = length(M) + 1))

#### Log-linear prior for observed cell probabilities, using conting ####
set.seed(60)
conting_data <- as.data.frame(rbind(c(rep(0, K), NA), kosovo_cont))
num_conting_burn <- 50000
num_conting <- 100000
# This model takes ~ 45 minutes to fit
conting_fit <- bict(Freq ~ (EXH + ABA + OSCE + HRW) ^ 3, data = conting_data, 
                    n.sample = num_conting_burn + num_conting, prior = "SBH", 
                    progress = TRUE)
pi_tilde_conting <- get_pi_tilde_conting(conting_fit, K, num_conting_burn)
pi_0_conting <- get_pi_0_hoi(pi_tilde_conting, K, xi = 1)
set.seed(61)
N_conting <- data.frame(N = numeric(), prior = numeric())
for(i in 1:length(M)){
    temp_accept <-  rejection_sampler_nbinom(pi_0_conting, M[i], a[i], n)
    temp_N <- n + n_0_sampler_nbinom(pi_0_conting[temp_accept], M[i], a[i], n)
    N_conting <- rbind(N_conting, data.frame(N = temp_N, prior = i))
}
N_conting <- rbind(N_conting, 
                   data.frame(N = n + n_0_sampler_isp(pi_0_conting, n),
                              prior = length(M) + 1))

#### Save the relevant output ####
save(pi_tilde_dir, pi_tilde_ll, pi_tilde_lcmcr, pi_tilde_conting, 
     file = "nhoi_fits/pi_tilde_samples.RData")
save(N_dir, N_ll, N_lcmcr, N_conting, file = "nhoi_fits/N_samples.RData")

