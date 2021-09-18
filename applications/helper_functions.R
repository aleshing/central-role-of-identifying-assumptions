# This file contains helper functions used to run analyses

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

# Get samples of observed cell probabilities from LCMCR output
get_pi_tilde_lcmcr <- function(fit, K){
    X <- matrix(as.numeric(intToBits(0:(2 ^ K - 1))), nrow = 2 ^ K, 
                byrow = TRUE)[, K:1]
    samp_iter <- fit$Get_Status()$buffer_size
    fit_log_nu <- fit$Get_Trace('log_nuK')
    fit_log_lambda <- fit$Get_Trace('log_lambdaJK2')
    fit_prob_zero <- fit$Get_Trace('prob_zero')
    fit_probs <- matrix(0, nrow = samp_iter, ncol = (2 ^ K) - 1)
    K_star <- fit$Get_Param("K")
    
    for(j in 3:(2 ^ K)){
        temp <- matrix(0, nrow = samp_iter, ncol = K_star)
        for(k in 1:K_star){
            temp[, k] <- fit_log_nu[, k] 
            for(l in 1:K){
                temp[, k] <- temp[, k] + X[j, l] * fit_log_lambda[, l, k, 2] +
                    (1 - X[j, l]) * fit_log_lambda[, l, k, 1]
            }
        }
        temp_max <- apply(temp, 1, max)
        fit_probs[, j - 1] <- exp(temp_max + log(rowSums(exp(temp - temp_max))))
    }
    fit_probs[, 1] <- 1 - (fit_prob_zero + rowSums(fit_probs))
    pi_tilde <- fit_probs / rowSums(fit_probs)
    return(pi_tilde)
}

# Get samples of observed cell probabilities from conting output
get_pi_tilde_conting <- function(fit, K, burn){
    ll_params <- fit$BETA[-(1:burn), ]
    X <- fit$maximal.mod$x
    num_samp <- nrow(ll_params)
    pi_tilde <- matrix(0, nrow = num_samp, ncol = (2 ^ K) - 1)
    for(i in 1:num_samp){
        pi_tilde[i, ] <- exp(ll_params[i, ] %*% t(X))[-1]
        pi_tilde[i, ] <- pi_tilde[i, ] / sum(pi_tilde[i, ])
    }
    return(pi_tilde)
}

# Given a matrix of posterior samples for the observed cell probabilities,
# pi_tilde_mat, where each row is a sample, calculate the missing cell 
# probabilities using a sensitivity parameter of xi
# xi = 1 corresponds to the no-highest-order interaction assumption
get_pi_0_hoi <- function(pi_tilde_mat, K, xi = exp(0)){
    X <- matrix(as.numeric(intToBits(0:(2 ^ K - 1))), nrow = 2 ^ K, 
                byrow = TRUE)[, K:1]
    temp <- rowSums(2 - X[-1, ]) %% 2
    evens <- which(temp == 0)
    odds <- which(temp == 1)
    if(K == 2){
        log_ratio <- rowSums(log(pi_tilde_mat[, odds])) - 
            log(pi_tilde_mat[, evens])
    }
    else{
        log_ratio <- rowSums(log(pi_tilde_mat[, odds])) - 
            rowSums(log(pi_tilde_mat[, evens]))
    }
    
    ratio <- exp(log_ratio)
    
    return(ratio / (xi + ratio))
}

# Given a matrix of posterior samples for the observed cell probabilities,
# pi_tilde_mat, where each row is a sample, calculate the missing cell 
# probabilities using an odds ratio of xi for the marginal 2x2 table
# of list_1 and list_2
# xi = 1 corresponds to the 2-list marginal no-highest-order interaction 
# assumption
get_pi_0_marg_hoi <- function(pi_tilde_mat, K, list_1, list_2, xi = exp(0)){
    X <- matrix(as.numeric(intToBits(0:(2 ^ K - 1))), nrow = 2 ^ K, 
                byrow = TRUE)[, K:1]
    
    X_marg <- matrix(as.numeric(intToBits(0:(2 ^ 2 - 1))), nrow = 2 ^ 2, 
                     byrow = TRUE)[, 2:1]
    pi_tilde_marg <- matrix(0, nrow = nrow(pi_tilde_mat), ncol = 4)
    
    X <- X[2:(2 ^ K), ]
    for(i in 1:nrow(X)){
        temp <- X[i, sort(c(list_1, list_2))]
        for(j in 1:nrow(X_marg)){
            if(temp[1] == X_marg[j, 1] & temp[2] == X_marg[j, 2]){
                pi_tilde_marg[, j] <- pi_tilde_marg[, j] + 
                    pi_tilde_mat[, i]
            }
        }
    }
    
    temp <- exp(log(pi_tilde_marg[, 2]) + log(pi_tilde_marg[, 3])  -
                    log(pi_tilde_marg[, 4]))
    numerator <- temp - xi * pi_tilde_marg[, 1]
    
    return(numerator / (xi + numerator))
}

# Given samples from the missing cell, pi_0, from the posterior for the observed 
# cell probabilities under just the conditional likelihood, run a rejection 
# sampler to obtain samples from the full posterior for pi_0 under the prior 
# N ~ NB(a, M/(M + a))
# Note this function outputs an accept/ reject decision for each sample,
# not the accepted samples
rejection_sampler_nbinom <- function(pi_0, M, a, n){
    # log(runif(length(pi_0))) < 
    #     dnbinom(n, size = a, prob = a / ((1 - pi_0) * M + a), log = TRUE) -
    #     dnbinom(n, size = a, prob = a / (n + a), log = TRUE)
    log(runif(length(pi_0))) < 
        log(pi_0 > 0) +
        (dnbinom(n, size = a, prob = a / ((1 - pi_0) * M + a), log = TRUE) -
             dnbinom(n, size = a, prob = a / (n + a), log = TRUE))
}

# Given samples from the missing cell, pi_0, under the full posterior for pi_0 
# under the prior N ~ NB(a, M/(M + a)), sample from the posterior for n_0
n_0_sampler_nbinom <- function(pi_0, M, a, n){
    rnbinom(length(pi_0), size = n + a, prob = 1 - (M / (M + a)) * pi_0) 
}

# Given samples from the missing cell, pi_0, under the full posterior for pi_0 
# under the improper scale prior for N, sample from the posterior for n_0
n_0_sampler_isp <- function(pi_0, n){
    rnbinom(length(pi_0), size = n, prob = 1 - pi_0) 
}

# Given a contingency table of data, cont, estimate the population size in a 
# frequentist framework using a highest order interaction of xi
# The last row of cont should correspond to the inclusion pattern (1, ..., 1)
# and the last column of cont should be named Freq
get_N_hoi_freq <- function(cont, xi = exp(0)){
    K <- ncol(cont) - 1
    
    # When specifying the offset, need to be careful and make sure you're using 
    # right parameterization since xi doesn't exactly correspond to 
    # exp(interaction) (i.e. sometimes it's exp(-interaction))
    cont$offset <- c(rep(0, (2 ^ K) - 2), ((-1) ^ (1 - K %% 2)) * log(xi))
    
    col_names <- colnames(cont)[1:K]
    if(K > 2){
        form <- paste0("Freq ~ (", paste(col_names, collapse = " + "),
                       ") ^ ", (K - 1))
    }
    else{
        form <- paste0("Freq ~ ", paste(col_names, collapse = " + "))
    }
    fit <- glm(form, data = cont, family = poisson, offset = offset)
    N_est <- unname(sum(cont$Freq) + exp(fit$coef[1]))
    var_intercept <- unname(summary(fit)$cov.unscaled[1, 1])
    # See Rivest and Levesque (2001) equation (5) for justification of se
    N_se  <- unname(sqrt(exp(fit$coef[1]) + 
                             (exp(2 * fit$coef[1])) * var_intercept))
    
    return(c(est = N_est, se = N_se, lower = N_est - 1.96 * N_se, 
             upper = N_est + 1.96 * N_se))
}