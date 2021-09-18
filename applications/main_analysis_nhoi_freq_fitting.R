library(LCMCR)
library(Rcapture)
library(dplyr)
source("helper_functions.R")

# This file contains code used to analyze the Kosovo data set using the 
# no-highest-order interaction assumption, in a frequentist framework

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

#### Fit model with NHOI assumption ####
rcapture_fit <- closedpMS.t(kosovo_cont, dfreq = TRUE)
N_freq <- c(est = rcapture_fit$result[1, 1], 
            se = rcapture_fit$results[1, 2],
            lower = rcapture_fit$results[1, 1] - 
                1.96 * rcapture_fit$results[1, 2],
            upper = rcapture_fit$results[1, 1] + 
                1.96 * rcapture_fit$results[1, 2])

#### Save the relevant output ####
save(rcapture_fit, N_freq, file = "nhoi_fits/Freq_fit.RData")