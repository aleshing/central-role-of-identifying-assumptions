library(dplyr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(xtable)
source("helper_functions.R")

# This file contains code used to summarize the sensitivity analysis for the 
# 2-list marginal no-highest-order interaction assumption for the Kosovo data 
# set

#### Load posterior samples ####
load("marg_nhoi_fits/N_samples_sensitivity.RData")

#### Load frequentist fit ####
load("marg_nhoi_fits/Freq_fit_sensitivity.RData")

#### Hyperparameters for negative-binomial priors for N ####
M <- c(8000, 12000, 10000)
a <- c(43, 9, 1.6)

#### Pre-processing for plots ####
N_lcmcr <- N_lcmcr %>% mutate(interaction = as.factor(interaction))
xis <- c("0.7", "0.8", "0.9", "1")
priors <- c("NB(8000, 43)", "NB(12000, 9)", "NB(10000, 1.6)", 
            "Improper Scale Prior")
N_lcmcr$which <- "Posterior"
for(i in 1:length(M)){
    temp <- data.frame(N = n_0_sampler_nbinom(rep(1, 100000), M = M[i],  
                                              a = a[i], n = 0),  
                       prior = i, which = "Prior")
    for(j in 1:length(xis)){
        temp$interaction <- xis[j]
        N_lcmcr  <- rbind(N_lcmcr , temp)
    }
}
N_lcmcr$prior <- as.factor(N_lcmcr$prior)
priors <- c("NB(8000, 43)", "NB(12000, 9)", "NB(10000, 1.6)", 
            "Improper Scale Prior")
levels(N_lcmcr$prior) <- priors
N_lcmcr <- N_lcmcr %>% 
    mutate(prior = fct_relevel(prior, rev(priors)))
N_lcmcr$interaction <- factor(N_lcmcr$interaction, levels = xis, 
                              ordered = TRUE, 
                              labels = c(expression(paste(xi , " = ", "0.7")),
                                         expression(paste(xi , " = ", "0.8")),
                                         expression(paste(xi , " = ", "0.9")),
                                         expression(paste(xi , " = ", "1"))))

#### Plot 1: Priors and Posteriors Plotted Together ####
N_lcmcr %>% 
    mutate(which = fct_relevel(as.factor(which), c("Prior", "Posterior"))) %>% 
    ggplot(aes(x = N, fill = which)) + geom_density(alpha = 0.5) + 
    theme_bw() + 
    facet_grid(interaction ~ prior,
               labeller = labeller(interaction = label_parsed)) +
    xlim(0, 35000) + theme(legend.position = "top") + 
    labs(y = "Density", fill = "Prior or Posterior?") 
ggsave("plots/sensitivity_analysis_marg_nhoi_plot1.pdf", height = 7, width = 8, 
       units = "in")

#### Plot 2: Results for both improper scale and negative-binomial priors ####
priors <- c("Improper Scale Prior", "Negative-Binomial", "NB(12000, 9)", 
            "NB(8000, 43)")
levels(N_lcmcr$prior) <- priors
N_lcmcr %>% 
    mutate(which = fct_relevel(as.factor(which), c("Prior", "Posterior"))) %>% 
    filter(which == "Posterior", prior %in% c("Negative-Binomial", 
                                              "Improper Scale Prior")) %>%
    ggplot(aes(x = N)) + 
    geom_density(aes(linetype = prior), alpha = 0.5) + 
    theme_bw() + 
    facet_grid(cols = vars(interaction),
               labeller = labeller(interaction = label_parsed)) + 
    xlim(6500, 17000) + theme(legend.position = "top") + 
    labs(y = "Density", linetype = "Prior for N") 
ggsave("plots/sensitivity_analysis_marg_nhoi_plot2.pdf", height = 4, width = 8, 
       units = "in")

#### Plot 3: Results for just negative-binomial prior ####
N_lcmcr %>% 
    mutate(which = fct_relevel(as.factor(which), c("Prior", "Posterior"))) %>% 
    filter(which == "Posterior", prior %in% c("Negative-Binomial")) %>%
    ggplot(aes(x = N)) + 
    geom_density(alpha = 0.5) + 
    theme_bw() + 
    facet_grid(cols = vars(interaction),
               labeller = labeller(interaction = label_parsed)) + 
    xlim(6500, 17000) + 
    labs(y = "Density", linetype = "Prior for N") 
ggsave("plots/sensitivity_analysis_marg_nhoi_plot3.pdf", height = 4, width = 8, 
       units = "in")

#### Plot 4: Plot Used in Paper ####
levels(N_lcmcr$interaction) <- xis
N_lcmcr$interaction <- fct_rev(N_lcmcr$interaction)
N_lcmcr %>% 
    mutate(which = fct_relevel(as.factor(which), c("Prior", "Posterior"))) %>% 
    filter(which == "Posterior", prior %in% c("Negative-Binomial")) %>%
    ggplot(aes(x = N)) + 
    geom_density(aes(linetype = interaction), alpha = 0.5) + 
    theme_bw() + 
    xlim(5000, 20000) + #  theme(legend.position = "top") +
    labs(y = "Density", linetype = expression(xi)) 
ggsave("plots/sensitivity_analysis_marg_nhoi_plot4.pdf", height = 4, width = 6, 
       units = "in")

#### Table Used in Paper ####
for_xtable <- N_lcmcr %>% filter(which == "Posterior",
                                 prior %in% c("Negative-Binomial")) %>%
    group_by(prior, interaction) %>%
    summarise(mean = mean(N), lower = quantile(N, probs = 0.025),
              median = median(N), upper = quantile(N, probs = 0.975)) %>%
    mutate(field = paste0(floor(mean), " [", floor(lower), ", ", floor(upper),
                          "]")) %>%
    ungroup() %>% select(field) %>% unlist() %>% matrix(nrow = 1, ncol = 4,
                                                        byrow = TRUE)
colnames(for_xtable) <- rev(xis)
xtable(for_xtable)

for_xtable_freq <- N_freq %>%
    mutate(field = paste0(floor(est), " [", floor(lower), ", ", floor(upper),
                          "]")) %>%
    ungroup() %>% select(field) %>% unlist() %>% matrix(nrow = 1, ncol = 4, 
                                                        byrow = TRUE) 
colnames(for_xtable_freq) <- rev(xis)

for_xtable_comb <- cbind(c("Frequentist", "Bayesian"), 
                         rbind(for_xtable_freq, for_xtable))
xtable(for_xtable_comb)

