## ---------------------------
##
## Script name: Correlations 
##
## Purpose of script: Calculating correlation between CWM species based and pollen based
##
## Author: Annegreet Veeken
##
## Date created: 08-02-2022
## ---------------------------
##
## Notes:
##   - re-run when leafC and leafN are available
##
## ---------------------------

## Load libraries
library(tidyverse)
library(rjags) 
library(runjags)
library(mcmc)
library(coda)
library(gridExtra)

## Load data and prepare
dfCWM <- readRDS("RDS_files/03_Posterior_CWM.rds") 

# Make df per trait, make tidy
lCWM <- dfCWM %>%
  dplyr::select(-id) %>% 
  group_by(trait) %>% 
  group_split(.keep = TRUE) %>% 
  purrr::map(., ~pivot_wider(., id_cols = c(trait, sitename, country), names_from = zone ,
                             values_from = Mean))



# dfCWM <- readRDS("RDS_files/03_Posterior_CWM.rds") %>%
#   as_tibble() %>%
#   pivot_wider(names_from = zone, values_from = c(Mean,SD))
# 
# dfCWM %>% 
#   filter(trait == "SLA") %>% 
#   ggplot(aes(x = exp(Mean_pollencorrected), y = Mean_zoneA)) +
#   geom_point(col = "darkorange") +
#   geom_smooth(method = "lm", col = "darkorange") 
# 
# 
# ggplot(dfCWM, aes(x = Mean_pollencorrected, y = Mean_zoneB)) +
#   geom_point(col = "darkorchid") +
#   geom_smooth(method = "lm", col = "darkorchid")+
#   facet_wrap(~trait)
# 
# ggplot(dfCWM, aes(x = Mean_pollencorrected, y = Mean_speciesC)) +
#   geom_point(col = "cyan4") +
#   geom_smooth(method = "lm", col = "cyan4")+
#   facet_wrap(~trait)

## Plots ----
plot_func <- function(df,traitname){
p <- ggplot(df, aes(x = pollencorrected, y =zoneA)) +
  geom_point(col = "darkorange") +
  geom_smooth(method = "lm", col = "darkorange")+
  geom_point(data = df, aes(x = pollencorrected, y = zoneB), col = "darkorchid") +
  geom_smooth(data = df, aes(x = pollencorrected, y = zoneB), method = "lm", col = "darkorchid")+
  geom_point(data = df, aes(x = pollencorrected, y = zoneC), col = "cyan4") +
  geom_smooth(data = df, aes(x = pollencorrected, y = zoneC), method = "lm", col = "cyan4")+
  scale_x_continuous("Pollen-based CWM") +
  scale_y_continuous("Vegetation CWM") +
  ggtitle(traitname) +
  theme_bw()
}

pall <- purrr::map2(lCWM[2:3], c("PlantHeight", "SLA"), ~plot_func(df = .x, traitname = .y))

trait_adj <- grid.arrange(pall[[1]], pall[[2]],
             # pall$LeafN,pall$LeafC,
             ncol = 2, nrow = 1)

plot_func <- function(df,traitname){
  p <- ggplot(df, aes(x = pollen, y = zoneA)) +
    geom_point(col = "darkorange") +
    geom_smooth(method = "lm", col = "darkorange")+
    geom_point(data = df, aes(x = pollen, y = zoneB), col = "darkorchid") +
    geom_smooth(data = df, aes(x = pollen, y = zoneB), method = "lm", col = "darkorchid")+
    geom_point(data = df, aes(x = pollen, y = zoneC), col = "cyan4") +
    geom_smooth(data = df, aes(x = pollen, y = zoneC), method = "lm", col = "cyan4")+
    scale_x_continuous("Pollen-based CWM") +
    scale_y_continuous("Vegetation CWM") +
    ggtitle(traitname) +
    theme_bw()
}

pall2 <- purrr::map2(lCWM[2:3],  c("PlantHeight", "SLA"), ~plot_func(df = .x, traitname = .y))

trait_pol <- grid.arrange(pall2[[1]], pall[[2]],
                          ncol = 2, nrow = 1)

## Regression bayesian----

# Trait <- dfCWM %>% 
#   filter(trait == "PlantHeight") %>% 
#   filter(zone %in% c("pollencorrected","zoneA")) %>% 
#   dplyr::select(sitename,zone, Mean) %>% 
#   pivot_wider(names_from = zone, values_from = Mean) %>% 
#   dplyr::select(-sitename) %>% 
#   as.matrix()
# N <- nrow(Trait)
# 
# # https://www.sumsar.net/blog/2013/08/bayesian-estimation-of-correlation/
# model_string <- "
#   model {
#     for(i in 1:N) {
#       Trait[i,1:2] ~ dmnorm(mu[], prec[ , ])
#     }
# 
#     # Constructing the covariance matrix and the corresponding precision matrix.
#     prec[1:2,1:2] <- inverse(cov[,])
#     cov[1,1] <- sigma[1] * sigma[1]
#     cov[1,2] <- sigma[1] * sigma[2] * rho
#     cov[2,1] <- sigma[1] * sigma[2] * rho
#     cov[2,2] <- sigma[2] * sigma[2]
#     
#     # Uninformative priors on all parameters which could, of course, be made more informative.
#     sigma[1] ~ dunif(0, 1000) 
#     sigma[2] ~ dunif(0, 1000)
#     rho ~ dunif(-1, 1)
#     mu[1] ~ dnorm(0, 0.001)
#     mu[2] ~ dnorm(0, 0.001)
# 
#     # Generate random draws from the estimated bivariate normal distribution
#     x_rand ~ dmnorm(mu[], prec[ , ])
# 
#     #data# N, Trait
#     #monitor# mu, rho, sigma, x_rand
#   }"
# 
# run.jags(model_string)
