# Code to Reproduce Applications

This directory contains code to reproduce the Kosovo data analyses presented in 
Section 6 and Appendix D. We describe how to reproduce these analyses below. 
Before doing so we note that the `R` script `helper_functions.R` contains helper 
functions used to run the analyses, and the `stan` file `log_linear_model.stan`
contains `stan` code used to fit a log-linear model with normal priors for the
log-linear parameters. 

To reproduce the analysis of the Kosovo data set presented in Section 6 using 
the 2-list marginal no-highest-order interaction assumption, one follows these 
steps:
1. Set your working directory in `R` to the current directory
2. Run the `R` script `main_analysis_marg_nhoi_fitting.R`. This produces 
posterior samples for the observed cell probabilities and the population size 
under the 2-list marginal no-highest-order interaction assumption, which are 
saved in the directory `marg_nhoi_fits`.
3. Run the `R` script `main_analysis_marg_nhoi_plotting.R`. This produces the 
plot and table presented in Section 6.4. The plot is saved in the folder 
`plots`.
4. Run the `R` script `sensitivity_analysis_marg_nhoi_fitting.R`. This produces 
posterior samples for the observed cell probabilities and the population size 
under the sensitivity analysis for the 2-list marginal no-highest-order 
interaction assumption, which are saved in the directory `marg_nhoi_fits`.
3. Run the `R` script `sensitivity_analysis_marg_nhoi_plotting.R`. This produces 
the plot and table presented in Section 6.5. The plot is saved in the folder 
`plots`.

To reproduce the analysis of the Kosovo data set presented in Appendix D using 
the no-highest-order interaction assumption, one follows these steps:
1. Set your working directory in `R` to the current directory
2. Run the `R` script `main_analysis_nhoi_fitting.R`. This produces posterior 
samples for the observed cell probabilities and the population size under the 
no-highest-order interaction assumption, which are saved in the directory 
`nhoi_fits`.
3. Run the `R` script `main_analysis_nhoi_plotting.R`. This produces the 
plot and table presented in Appendix D.2. The plot is saved in the folder 
`plots`.
4. Run the `R` script `sensitivity_analysis_nhoi_fitting.R`. This produces 
posterior samples for the observed cell probabilities and the population size 
under the sensitivity analysis for the no-highest-order interaction assumption, 
which are saved in the directory `nhoi_fits`.
3. Run the `R` script `sensitivity_analysis_nhoi_plotting.R`. This produces 
the plot and table presented in Appendix D.3. The plot is saved in the folder 
`plots`.