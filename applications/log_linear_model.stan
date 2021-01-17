data 
{
    // Data
    int<lower=0> K;               // Number of lists
    int<lower=0> twoK;            // 2 ^ K
    vector<lower=0>[twoK - 1] ns; // Observed cell counts
    int<lower=0> n;               // Observed sample size, sum(ns)
    int<lower=K> p;               // Number of log-linear parameters
    matrix[twoK - 1, p] X;        // Design matrix for log-linear parameters
    
    // Prior hyperparmaters
    real<lower=0> ll_sd;          // Standard deviation of prior for log-linear 
                                  // parameters
}
parameters 
{
    vector[p] lambda_raw;         // Auxiliary variable for log-linear 
                                  // parameters
}
transformed parameters{
    vector[twoK - 1] log_mu;      
    real log_mu_sum;          
    vector[p] lambda;             // Log-linear parameters
    
    lambda = lambda_raw * ll_sd;
    log_mu = X * lambda;
    log_mu_sum = log_sum_exp(log_mu);
}
model 
{
    lambda_raw ~ std_normal();
    target += dot_product(ns, log_mu) - n * log_mu_sum;
}
generated quantities{
    vector[twoK - 1] pi_tilde;    // Observed cell probabilities
    {
        real mu_sum = exp(log_mu_sum);
        pi_tilde = exp(log_mu) / mu_sum;
    }
}
