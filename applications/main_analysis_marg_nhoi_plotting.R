library(dplyr)
library(ggplot2)
library(ggpubr)
library(forcats)
library(xtable)
source("helper_functions.R")

# This file contains code used to summarize the analysis of the Kosovo data set 
# using the 2-list marginal no-highest-order interaction assumption

#### Load posterior samples ####
load("marg_nhoi_fits/N_samples.RData")

#### Hyperparameters for negative-binomial priors for N ####
M <- c(8000, 12000, 10000)
a <- c(43, 9, 1.6)

#### Pre-processing for plots ####
models <- c("Dirichlet", "Log-Linear", "LCMCR", "Conting")
N_samples <- rbind(cbind(N_dir, model = models[1]),
                   cbind(N_ll, model = models[2]),
                   cbind(N_lcmcr, model = models[3]),
                   cbind(N_conting, model = models[4]))
N_samples$which <- "Posterior"
for(i in 1:length(M)){
    temp <- data.frame(N = n_0_sampler_nbinom(rep(1, 100000), M = M[i],  
                                              a = a[i], n = 0),  
                       prior = i)
    for(j in 1:length(models)){
        temp$model <- models[j]
        temp$which <- "Prior"
        N_samples <- rbind(N_samples, temp)
    }
}
N_samples$prior <- as.factor(N_samples$prior)
priors <- c("NB(8000, 43)", "NB(12000, 9)", "NB(10000, 1.6)", 
            "Improper Scale Prior")
levels(N_samples$prior) <- priors
N_samples <- N_samples %>% 
    mutate(prior = fct_relevel(prior, rev(priors)))

#### Plot 1: Priors and Posteriors Plotted Together ####
N_samples %>% 
    mutate(which = fct_relevel(as.factor(which), c("Prior", "Posterior"))) %>% 
    # filter(which == "Posterior") %>%
    ggplot(aes(x = N, fill = which)) + geom_density(alpha = 0.5) + 
    theme_bw() + facet_grid(rows = vars(model), cols = vars(prior)) + 
    xlim(0, 30000) + theme(legend.position = "top") + 
    labs(y = "Density", fill = "Prior or Posterior?") 
ggsave("plots/main_analysis_marg_nhoi_plot1.pdf", height = 7, width = 8, 
       units = "in")

#### Plot 2: Plot Used in Paper ####
priors <- c("Improper Scale Prior", "Negative-Binomial", "NB(12000, 9)", 
            "NB(8000, 43)")
levels(N_samples$prior) <- priors
N_samples %>% 
    mutate(which = fct_relevel(as.factor(which), c("Prior", "Posterior"))) %>% 
    filter(which == "Posterior", prior %in% c("Negative-Binomial", 
                                              "Improper Scale Prior")) %>%
    ggplot(aes(x = N)) + 
    geom_density(aes(linetype = prior), alpha = 0.5) + 
    theme_bw() + facet_grid(cols = vars(model)) + 
    xlim(6000, 13000) + theme(legend.position = "top") + 
    labs(y = "Density", linetype = "Prior for N") 
ggsave("plots/main_analysis_marg_nhoi_plot2.pdf", height = 4, width = 8, 
       units = "in")

#### Table Used in Paper ####
for_xtable <- N_samples %>% filter(which == "Posterior", 
                                   prior %in% c("Negative-Binomial", 
                                                "Improper Scale Prior")) %>% 
    group_by(prior, model) %>% 
    summarise(mean = mean(N), lower = quantile(N, probs = 0.025),
              median = median(N), upper = quantile(N, probs = 0.975)) %>%
    mutate(field = paste0(floor(mean), " [", floor(lower), ", ", floor(upper),
                          "]")) %>% 
    ungroup() %>% select(field) %>% unlist() %>% matrix(nrow = 4, ncol = 2, 
                                                        byrow = FALSE) 
rownames(for_xtable) <- models
colnames(for_xtable) <- c("Improper Scale Prior", "Negative-Binomial")
xtable(for_xtable)