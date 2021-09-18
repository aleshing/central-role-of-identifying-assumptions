library(LCMCR)
library(Rcapture)
library(dplyr)
source("helper_functions.R")

# This file contains code used to perform a sensitivity analysis for the 2-list 
# marginal no-highest-order interaction assumption for the Kosovo data set, in a
# frequentist framework

#### Read in Kosovo data from the LCMCR package ####
data(kosovo_aggregate)
kosovo_cont <- kosovo_aggregate %>% group_by(EXH, ABA, OSCE, HRW) %>% 
    summarise(Freq = n()) %>% ungroup() %>% as.data.frame() %>%
    group_by(ABA, HRW) %>% summarise(Freq = sum(Freq)) %>% 
    ungroup() %>% as.data.frame() %>%
    mutate(ABA = as.numeric(as.character(ABA)),
           HRW = as.numeric(as.character(HRW))) %>%
    filter(ABA == 1 | HRW == 1)

K <- ncol(kosovo_cont) - 1
n <- sum(kosovo_cont$Freq)

#### Range of values for sensitivity parameter ####
xi <- c(1, 0.9, 0.8, 0.7)
N_freq <- data.frame(interaction = NA, est = NA, se = NA, lower = NA, 
                     upper = NA)

#### Fit model with HOI given in xi ####
for(i in 1:length(xi)){
    fit <- c(interaction = xi[i], get_N_hoi_freq(kosovo_cont, xi = xi[i]))
    N_freq <- rbind(N_freq, fit)
}
N_freq <- N_freq[-1, ]

#### Save the relevant output ####
save(N_freq, file = "marg_nhoi_fits/Freq_fit_sensitivity.RData")