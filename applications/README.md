# Code to Reproduce Applications

This directory contains code to reproduce the Kosovo data analyses presented in 
Section 5 and Appendix E. We describe how to reproduce these analyses below. 
Before doing so we note that the `R` script `helper_functions.R` contains helper 
functions used to run the analyses, and the `stan` file `log_linear_model.stan`
contains `stan` code used to fit a log-linear model with normal priors for the
log-linear parameters. 

To reproduce the main analyses of the Kosovo data set, one follows these 
steps:
1. Set your working directory in `R` to the current directory
2. Run the `R` script `main_analysis_marg_nhoi_fitting.R`. This produces 
posterior samples for the observed cell probabilities and the population size in a 
Bayesian framework under the 2-list marginal no-highest-order interaction assumption, 
which are  saved in the directory `marg_nhoi_fits`.
3. Run the `R` script `main_analysis_marg_nhoi_freq_fitting.R`. This produces 
estimates of the population size in a frequentist framework under the 2-list marginal 
no-highest-order interaction assumption, which are  saved in the directory 
`marg_nhoi_fits`.
4. Run the `R` script `main_analysis_nhoi_fitting.R`. This produces posterior 
samples for the observed cell probabilities and the population size in a 
Bayesian framework under the  no-highest-order interaction assumption, which are
saved in the directory `nhoi_fits`.
5. Run the `R` script `main_analysis_nhoi_freq_fitting.R`. This produces 
estimates of the population size in a frequentist framework under the  no-highest-order 
interaction assumption, which are saved in the directory `nhoi_fits`.
6. Run the `R` script `main_analysis_marg_nhoi_plotting.R`. This produces 
Figure 1 and Table 17 presented in Appndix E.2. The plot is saved in the folder 
`plots`.
7. Run the `R` script `main_analysis_nhoi_plotting.R`. This produces 
Figure 2 and Table 18 presented in Appndix E.2. The plot is saved in the folder 
`plots`.

To reproduce the sensitivity analyses of the Kosovo data set, one follows these 
steps:
1. Run the `R` script `sensitivity_analysis_marg_nhoi_fitting.R`. This produces 
posterior samples for the observed cell probabilities and the population size  in a 
Bayesian framework under the sensitivity analysis for the 2-list marginal no-highest-order 
interaction assumption, which are saved in the directory `marg_nhoi_fits`.
2. Run the `R` script `sensitivity_analysis_marg_nhoi_freq_fitting.R`. This 
produces estimates of the population size in a frequentist framework under the sensitivity 
analysis for the 2-list marginal no-highest-order interaction assumption, which are saved 
in the directory `marg_nhoi_fits`.
3. Run the `R` script `sensitivity_analysis_nhoi_fitting.R`. This produces 
posterior samples for the observed cell probabilities and the population size in a 
Bayesian framework under the sensitivity analysis for the no-highest-order interaction 
assumption, which are saved in the directory `nhoi_fits`.
4. Run the `R` script `sensitivity_analysis_nhoi_freq_fitting.R`. This 
produces estimates of the population size in a frequentist framework under the 
sensitivity analysis for the no-highest-order interaction assumption, which are saved in 
the directory `nhoi_fits`.
5. Run the `R` script `sensitivity_analysis_marg_nhoi_plotting.R`. This produces 
Table 4 presented in Section 5.4. 
6. Run the `R` script `sensitivity_analysis_nhoi_plotting.R`. This produces 
Table 19 presented in Appendix E.3. 
