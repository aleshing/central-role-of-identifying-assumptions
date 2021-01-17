library(LCMCR)
source("helper_functions.R")

# Example 2 using lists 3 and 4 from simulation of Manrique-Vallier (2016)

#### Set latent class model parameters ####

# Number of lists
K <- 2
# Number of classes
J <- 2
# Class probabilities 
nu <- c(0.9, 0.1)
# Class specific observation probabilities 
q <- matrix(c(0.099, 0.132, 0.759, 0.99), ncol = K, byrow = TRUE)
# Cell probabilities
pis <- exp(get_log_pis(nu, q, K))
# Missing cell probabilities
pi_0 <- pis[1]
# Conditional cell probabilities
pi_tilde <- pis[2:4] / (1- pi_0)

#### Run simulation ####

# Simulation settings
num_rep <- 200
num_dpmm_burn <- 50000
num_dpmm <- 200000
Ns <- c(2000, 10000, 100000)
results_N <- vector("list", length(Ns))
results_pi_0 <- vector("list", length(Ns))
# Run Simulation
for(i in 1:length(Ns)){
    results_N[[i]] <- matrix(0, nrow = num_rep, ncol = 6)
    results_pi_0[[i]] <- matrix(0, nrow = num_rep, ncol = 6)
    for(j in 1:num_rep){
        print(paste0("Iteration ", j, ", N = ", Ns[i]))
        data <- generate_lcm_data(Ns[i], nu, q, K, seed = j + i)
        data_obs <- data[rowSums(data) > 0,]
        fit_dpmm <- fit_lcmcr(data_obs, K, burn_iter = num_dpmm_burn, 
                              samp_iter = num_dpmm, K_star = 2, seed = 64, 
                              tab = FALSE, a = 0.25, b = 0.25)
        pi_0_dpmm <- fit_dpmm$Get_Trace('prob_zero')
        N_dpmm <- nrow(data_obs) + fit_dpmm$Get_Trace('n0')
        results_N[[i]][j, ] <- c(mean(N_dpmm), 
                                 quantile(N_dpmm, 
                                          prob = c(0.5, 0.025, 0.975,
                                                   0.25, 0.75)))
        results_pi_0[[i]][j, ] <- c(mean(pi_0_dpmm), 
                                    quantile(pi_0_dpmm, 
                                             prob = c(0.5, 0.025, 0.975,
                                                      0.25, 0.75)))
        if(j %% 10 == 0){
            save(results_N, results_pi_0, file = "example_2_simulation.RData")
        }
    }
}

#### Summarize Simulation ####

for(i in 1:length(Ns)){
    print(paste0("Example 2: Population Size " , Ns[i], " Results"))
    print(paste0("Mean Posterior Median: ", mean(results_pi_0[[i]][, 2])))
    print(paste0("95% Credible Interval Coverage: ", 
                 mean(results_pi_0[[i]][, 4] > pi_0 & 
                          results_pi_0[[i]][, 3] < pi_0)))
    print(paste0("Mean 95% Credible Interval Width: ", 
                 mean(results_pi_0[[i]][, 4] - results_pi_0[[i]][, 3])))
    print(paste0("50% Credible Interval Coverage: ", 
                 mean(results_pi_0[[i]][, 6] > pi_0 & 
                          results_pi_0[[i]][, 5] < pi_0)))
    print(paste0("Mean 50% Credible Interval Width: ", 
                 mean(results_pi_0[[i]][, 6] - results_pi_0[[i]][, 5])))
}
