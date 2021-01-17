# This file contains helper functions used to run the simulations

# Fit a Bayesian latent class model in LCMCR
fit_lcmcr <- function(cont, K, burn_iter = 10000, samp_iter = 10000,
                      K_star = 10, seed = 42, tabular = TRUE, a = 0.25, 
                      b = 0.25){
    if(tabular){
        colnames(cont)[K + 1] <- "Freq"
        fit <- lcmCR(captures = cont, tabular = TRUE, 
                     in_list_label = '1', not_in_list_label = '0', K = K_star, 
                     a_alpha = a, b_alpha = b, seed = seed,
                     buffer_size = samp_iter, thinning = 1)
    }
    else{
        fit <- lcmCR(captures = cont, tabular = FALSE, in_list_label = '1', 
                     not_in_list_label = '0', K = K_star, a_alpha = a, 
                     b_alpha = b, seed = seed, buffer_size = samp_iter, 
                     thinning = 1)
    }
    
    fit$Update(burn_iter)
    fit$Set_Trace(c('n0', 'prob_zero', 'log_nuK', 'log_lambdaJK2'))
    fit$Activate_Tracing()
    fit$Update(samp_iter)
    return(fit)
}

# Calculate pis, i.e. the cell probabilities, for a given set of latent class 
# model parameters, on the log scale
get_log_pis <- function(nu, q, K){
    X <- matrix(as.numeric(intToBits(0:(2 ^ K - 1))), nrow = 2 ^ K, 
                byrow = TRUE)[, K:1]
    J <- length(nu)
    log_pis <- rep(0, nrow(X))
    for(i in 1:nrow(X)){
        temp_prod <- rep(0, J)
        for(j in 1:J){
            temp_prod[j]<- sum(X[i, ] * log(q[j, ])) +
                sum((1 - X[i, ]) * log(1 - q[j, ])) + log(nu[j])
        }
        max_temp_prod <- max(temp_prod)
        log_pis[i] <- max_temp_prod + 
            log(sum(exp(temp_prod - max_temp_prod)))
    }
    return(log_pis)
}

# Generate data from a J class latent class model
generate_lcm_data <- function(N, nu, q, K, seed = 42){ 
    set.seed(seed)
    
    J <- length(nu)
    data <- matrix(0, nrow = N, ncol = K)
    zs <- t(rmultinom(N, size = 1, prob = nu))
    
    for(i in 1:N){
        for(j in 1:J){
            if(zs[i, j] == 1){
                for(k in 1:K){
                    data[i, k] <- rbinom(1, size = 1, prob = q[j, k])
                }
            }
        }
    }
    return(data)
}