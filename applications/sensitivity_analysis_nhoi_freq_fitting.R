library(LCMCR)
library(Rcapture)
library(dplyr)
source("helper_functions.R")

# This file contains code used to perform a sensitivity analysis for the 
# no-highest-order interaction assumption for the Kosovo data set, in a 
# frequentist framework

#### Read in Kosovo data from the LCMCR package ####
data(kosovo_aggregate)
K <- ncol(kosovo_aggregate)
n <- nrow(kosovo_aggregate)
kosovo_cont <- kosovo_aggregate %>% group_by(EXH, ABA, OSCE, HRW) %>% 
    summarise(Freq = n()) %>% ungroup() %>% as.data.frame() %>%
    mutate(EXH = as.numeric(as.character(EXH)),
           ABA = as.numeric(as.character(ABA)),
           OSCE = as.numeric(as.character(OSCE)),
           HRW = as.numeric(as.character(HRW)))

#### Range of values for sensitivity parameter ####
xi <- c(1 / 2, 2 / 3, 1, 3 / 2, 2)
N_freq <- data.frame(interaction = NA, est = NA, se = NA, lower = NA, 
                     upper = NA)

#### Fit model with HOI given in xi ####
for(i in 1:length(xi)){
    fit <- c(interaction = xi[i], get_N_hoi_freq(kosovo_cont, xi = xi[i]))
    N_freq <- rbind(N_freq, fit)
}
N_freq <- N_freq[-1, ]

#### Save the relevant output ####
save(N_freq, file = "nhoi_fits/Freq_fit_sensitivity.RData")