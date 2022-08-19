## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Annegreet Veeken
##
## Date Created: 2022-07-24
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


## Load packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(rjags)) install.packages("rjags")
if(!require(runjags)) install.packages("runjags")
if(!require(mvtnorm)) install.packages("mvtnorm") #to generate correlated data with rmvnorm.
if(!require(car)) install.packages("car") #To plot the estimated bivariate normal distribution.
library(gridExtra)

# collate data
traits <- c("LA","SLA","PlantHeight")

pollen_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "03_CWM_estimates_") %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(pattern = "pollen") 
  
pollen_lab <- pollen_files %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_remove(., "03_CWM_estimates_") %>% 
  str_remove(., ".rds")

veg_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "03_CWM_estimates_") %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(pattern = "pollen", negate = TRUE)
veg_lab <- veg_files %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_remove(., "03_CWM_estimates_") %>% 
  str_remove(., ".rds")

folderpath.fun <- function(x)
  {paste("RDS_files/", x, sep = "/")}

dfCWM_pol <- pollen_files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(pollen_lab, ~readRDS(.) %>% 
                as.data.frame() %>% 
                rownames_to_column("parameter") %>% 
                filter(., str_detect(parameter, "cwm\\[")) %>%
               
                mutate(label = .y,
                       id = seq(1:nrow(.))
                       )
              ) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  dplyr::select(sitename,label, pollen_mean = Mean, pollen_sd = SD) %>% 
  mutate(trait = str_extract(label, pattern = paste(traits, collapse = "|")) %>% 
           str_remove_all("_"),
         sitename = str_remove(sitename, pattern = "X"),
         growthform = case_when(str_detect(label, pattern = "herb") ~ "herb",
                                str_detect(label, pattern = "trsh") ~ "trsh",
                                TRUE ~ "allpft"),
         pollination = case_when(str_detect(label, pattern = "not wind") ~ "not wind",
                                 str_detect(label, pattern = "wind") ~ "wind",
                                 TRUE ~ "allpol"),
         correction = case_when(str_detect(label, pattern = "adjustedpercent_mean") ~ "adjusted",
                                str_detect(label, pattern = "trsh|herb") ~ "adjusted",
                                str_detect(label, pattern = "wind") ~ "adjusted",
                            str_detect(label, pattern = "_percent") ~ "non",
                            str_detect(label, pattern = "_draw") ~ "draw"),
         country = str_extract(label, pattern = "Scotland|Switzerland")
  ) 

dfCWM_veg <- veg_files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(veg_lab, ~readRDS(.) %>% 
                as.data.frame() %>% 
                rownames_to_column("parameter") %>% 
                filter(., str_detect(parameter, "cwm\\[")) %>%
                mutate(label = .y,
                       id = seq(1:nrow(.))) %>% 
                mutate(across(one_of("sitename"),  as.character))) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(zone = case_when(str_detect(label,"Switzerland") ~ "BDM",
                          TRUE ~ zone)) %>% 
  dplyr::select(sitename, trait, veg_mean = Mean, veg_sd = SD,
                growthform, pollination, taxres, zone) 

 
cor_bay <- function(cwm, selectedzone, selectedtrait){
  CWM <- cwm %>% 
    filter(zone == selectedzone,
           trait == selectedtrait)
  CWM_mean <- CWM %>% 
    dplyr::select(pollen_mean, veg_mean) %>% 
    as.matrix() 
  CWM_sd <- CWM %>% 
    dplyr::select(pollen_sd, veg_sd) %>% 
    as.matrix()
  N <- nrow(CWM_mean)
  
  # https://www.sumsar.net/blog/2013/08/bayesian-estimation-of-correlation/
  model_string <- "
    model {
      for(i in 1:N) {
        CWM_mean[i,1:2] ~ dmnorm(mu[], prec[,,i])
      }
  
      # Constructing the covariance matrix and the corresponding precision matrix.
      for(i in 1:N){
        prec[1:2,1:2,i] <- inverse(cov[,,i])
        cov[1,1,i] <- CWM_sd[i,1] * CWM_sd[i,1]
        cov[1,2,i] <- CWM_sd[i,1] * CWM_sd[i,2] * rho
        cov[2,1,i] <- CWM_sd[i,1] * CWM_sd[i,2] * rho
        cov[2,2,i] <- CWM_sd[i,2] * CWM_sd[i,2]
      }
  
      # Priors
      rho ~ dunif(-1, 1)
      mu[1] ~ dnorm(0, 0.001)
      mu[2] ~ dnorm(0, 0.001)
  
      # Generate random draws from the estimated bivariate normal distribution
      for(i in 1:N){
        x_rand[i,1:2] ~ dmnorm(mu[], prec[,,i])
        }
      #data# N, CWM_mean, CWM_sd
      #monitor# rho
    }"
  
  inits_list <-  list(mu = c(min(CWM_mean)*0.1, max(CWM_mean)*10),
                    rho = cor(CWM_mean[, 1], CWM_mean[, 2]))
  
  results <- run.jags(model_string, inits = inits_list)
  results_summary <- as.data.frame(summary(results)) %>% 
    rownames_to_column("parameter") %>% 
    mutate(trait = selectedtrait,
           zone = selectedzone)
  return(results_summary)
  }

## All pollen no correction ----
cwm <- dfCWM_pol %>% 
  filter(growthform == "allpft",
       pollination == "allpol",
       correction == "non") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 
combo <- expand_grid(trait = unique(dfCWM_pol$trait), zone = unique(dfCWM_veg$zone))
rhos <- purrr::map2_dfr(combo$zone, combo$trait, ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y))
(p_scatter <- 
  ggplot(cwm, aes(x = pollen_mean, y = veg_mean, color = zone)) +
  geom_point() +
  scale_fill_manual() + 
  scale_x_continuous("Pollen") +
  scale_y_continuous("Vegetation") +
  facet_wrap(~trait, scales = "free") +
  ggtitle("Uncorrected pollen data") + 
  theme_bw() + 
  theme(legend.position =  "bottom"))

(p_rho <- ggplot(rhos, aes(x = zone, y=Mean)) + 
  geom_errorbar(aes(ymin=Lower95, ymax=Upper95), width=.1) +
  geom_point(shape = 21, color = "black", fill = "darkorange",
             size = 3) + 
  scale_x_discrete("") + 
  scale_y_continuous("", limits = c(-1,1)) +
  facet_wrap(~trait) +
  coord_flip() +
  theme_bw()
)

pol <- grid.arrange(p_scatter, p_rho)
ggsave("Figures/Pollen_uncorrected.png", pol, width = 7, height = 7)

## Pollen correction ----
cwm <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "adjusted") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 
combo <- expand_grid(trait = unique(dfCWM_pol$trait), zone = unique(dfCWM_veg$zone))
rhos <- purrr::map2_dfr(combo$zone, combo$trait, ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y))
(p_scatter <- 
    ggplot(cwm, aes(x = pollen_mean, y = veg_mean, color = zone)) +
    geom_point() +
    scale_fill_manual() + 
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    facet_wrap(~trait, scales = "free") +
    ggtitle("Corrected pollen data") + 
    theme_bw() + 
    theme(legend.position =  "bottom"))

(p_rho <- ggplot(rhos, aes(x = zone, y=Mean)) + 
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95), width=.1) +
    geom_point(shape = 21, color = "black", fill = "darkorange",
               size = 3) + 
    scale_x_discrete("") + 
    scale_y_continuous("", limits = c(-1,1)) +
    facet_wrap(~trait) +
    coord_flip() +
    theme_bw()
)

pol <- grid.arrange(p_scatter, p_rho)
ggsave("Figures/Pollen_corrected.png", pol, width = 7, height = 7)

## Pollen correction + uncertainty ----
cwm <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "draw", 
         country == "Scotland") %>% 
  mutate(sitename = rep(c("C001", "C002", "C003","C004", "C005", "C006","C007", "C008", 
                      "C009","C010", "C011", "C012","C013", "C014", "C015", "C016"),300),
         draw = str_extract(label, pattern = "[draw]")) %>%  
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 

cor_bay_draws <- function(cwm, selectedzone, selectedtrait){
  CWM <- cwm %>% 
    filter(zone == selectedzone,
           trait == selectedtrait)
  CWM_mean <- CWM %>% 
    dplyr::select(pollen_mean, veg_mean) %>% 
    as.matrix() 
  CWM_sd <- CWM %>% 
    dplyr::select(pollen_sd, veg_sd) %>% 
    as.matrix()
  N <- nrow(CWM_mean)
  
  # https://www.sumsar.net/blog/2013/08/bayesian-estimation-of-correlation/
  model_string <- "
    model {
      for(i in 1:N) {
        CWM_mean[i,1:2] ~ dmnorm(mu[], prec[,,i])
      }
  
      # Constructing the covariance matrix and the corresponding precision matrix.
      for(i in 1:N){
        prec[1:2,1:2,i] <- inverse(cov[,,i])
        cov[1,1,i] <- CWM_sd[i,1] * CWM_sd[i,1]
        cov[1,2,i] <- CWM_sd[i,1] * CWM_sd[i,2] * rho
        cov[2,1,i] <- CWM_sd[i,1] * CWM_sd[i,2] * rho
        cov[2,2,i] <- CWM_sd[i,2] * CWM_sd[i,2]
      }
  
      # Priors
      rho ~ dunif(-1, 1)
      mu[1] ~ dnorm(0, 0.001)
      mu[2] ~ dnorm(0, 0.001)
  
      # Generate random draws from the estimated bivariate normal distribution
      for(i in 1:N){
        x_rand[i,1:2] ~ dmnorm(mu[], prec[,,i])
        }
      #data# N, CWM_mean, CWM_sd
      #monitor# rho
    }"
  
  inits_list <-  list(mu = c(min(CWM_mean)*0.1, max(CWM_mean)*10),
                      rho = cor(CWM_mean[, 1], CWM_mean[, 2]))
  
  results <- run.jags(model_string, inits = inits_list)
  results_summary <- as.data.frame(summary(results)) %>% 
    rownames_to_column("parameter") %>% 
    mutate(trait = selectedtrait,
           zone = selectedzone)
  return(results_summary)
}


combo <- expand_grid(trait = unique(dfCWM_pol$trait), zone = unique(dfCWM_veg$zone))
rhos <- purrr::map2_dfr(combo$zone, combo$trait, ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y))

## Pollination mode ----
# Wind
cwm <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "wind",
         correction == "adjusted") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 
combo <- expand_grid(trait = unique(dfCWM_pol$trait), zone = unique(dfCWM_veg$zone))
rhos <- purrr::map2_dfr(combo$zone, combo$trait, ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y))
(p_scatter <- 
    ggplot(cwm, aes(x = pollen_mean, y = veg_mean, color = zone)) +
    geom_point() +
    scale_fill_manual() + 
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    facet_wrap(~trait, scales = "free") +
    ggtitle("wind pollinated taxa") + 
    theme_bw() + 
    theme(legend.position =  "bottom"))

(p_rho <- ggplot(rhos, aes(x = zone, y=Mean)) + 
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95), width=.1) +
    geom_point(shape = 21, color = "black", fill = "darkorange",
               size = 3) + 
    scale_x_discrete("") + 
    scale_y_continuous("", limits = c(-1,1)) +
    facet_wrap(~trait) +
    coord_flip() +
    theme_bw()
)

pol <- grid.arrange(p_scatter, p_rho)
ggsave("Figures/Wind_pollinated.png", pol, width = 7, height = 7)
# No wind
cwm <- dfCWM_pol %>%
  filter(growthform == "allpft",
         pollination == "not wind",
         correction == "adjusted")
cwm$pollination[cwm$pollination == "not wind"] <- "no wind"
cwm <- cwm  %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>%
  filter(taxres == "stand.spec")
combo <- expand_grid(trait = unique(dfCWM_pol$trait), zone = unique(dfCWM_veg$zone))
rhos <- purrr::map2_dfr(combo$zone, combo$trait, ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y))
(p_scatter <-
    ggplot(cwm, aes(x = pollen_mean, y = veg_mean, color = zone)) +
    geom_point() +
    scale_fill_manual() +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    facet_wrap(~trait, scales = "free") +
    ggtitle("Not wind pollinated taxa") +
    theme_bw() +
    theme(legend.position =  "bottom"))

(p_rho <- ggplot(rhos, aes(x = zone, y=Mean)) +
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95), width=.1) +
    geom_point(shape = 21, color = "black", fill = "darkorange",
               size = 3) +
    scale_x_discrete("") +
    scale_y_continuous("", limits = c(-1,1)) +
    facet_wrap(~trait) +
    coord_flip() +
    theme_bw()
)

pol <- grid.arrange(p_scatter, p_rho)
ggsave("Figures/NotWind_pollinated.png", pol, width = 7, height = 7)

## Plant functional type ----
# Trees and Shrubs
cwm <- dfCWM_pol %>% 
  filter(growthform == "trsh",
         pollination == "allpol",
         correction == "adjusted") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 
combo <- expand_grid(trait = unique(dfCWM_pol$trait), zone = unique(dfCWM_veg$zone))
rhos <- purrr::map2_dfr(combo$zone, combo$trait, ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y))
(p_scatter <- 
    ggplot(cwm, aes(x = pollen_mean, y = veg_mean, color = zone)) +
    geom_point() +
    scale_fill_manual() + 
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    facet_wrap(~trait, scales = "free") +
    ggtitle("trees and shrub taxa") + 
    theme_bw() + 
    theme(legend.position =  "bottom"))

(p_rho <- ggplot(rhos, aes(x = zone, y=Mean)) + 
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95), width=.1) +
    geom_point(shape = 21, color = "black", fill = "darkorange",
               size = 3) + 
    scale_x_discrete("") + 
    scale_y_continuous("", limits = c(-1,1)) +
    facet_wrap(~trait) +
    coord_flip() +
    theme_bw()
)

pol <- grid.arrange(p_scatter, p_rho)
ggsave("Figures/TRSH.png", pol, width = 7, height = 7)
# Herbs
cwm <- dfCWM_pol %>% 
  filter(growthform == "herb",
         pollination == "allpol",
         correction == "adjusted") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 
combo <- expand_grid(trait = unique(dfCWM_pol$trait), zone = unique(dfCWM_veg$zone))
rhos <- purrr::map2_dfr(combo$zone, combo$trait, ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y))
(p_scatter <- 
    ggplot(cwm, aes(x = pollen_mean, y = veg_mean, color = zone)) +
    geom_point() +
    scale_fill_manual() + 
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    facet_wrap(~trait, scales = "free") +
    ggtitle("herb taxa") + 
    theme_bw() + 
    theme(legend.position =  "bottom"))

(p_rho <- ggplot(rhos, aes(x = zone, y=Mean)) + 
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95), width=.1) +
    geom_point(shape = 21, color = "black", fill = "darkorange",
               size = 3) + 
    scale_x_discrete("") + 
    scale_y_continuous("", limits = c(-1,1)) +
    facet_wrap(~trait) +
    coord_flip() +
    theme_bw()
)

pol <- grid.arrange(p_scatter, p_rho)
ggsave("Figures/HERB.png", pol, width = 7, height = 7)
## Taxonomic resolution ----
# Genus
cwm <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "adjusted") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "genus") %>% 
  filter(!zone == "BDM") # temporary filter till convergence issue is fixed
combo <- expand_grid(trait = unique(dfCWM_pol$trait), zone = unique(dfCWM_veg$zone)[-1])
rhos <- purrr::map2_dfr(combo$zone, combo$trait, ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y))
(p_scatter <- 
    ggplot(cwm, aes(x = pollen_mean, y = veg_mean, color = zone)) +
    geom_point() +
    scale_fill_manual() + 
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    facet_wrap(~trait, scales = "free") +
    ggtitle("Genus level") + 
    theme_bw() + 
    theme(legend.position =  "bottom"))

(p_rho <- ggplot(rhos, aes(x = zone, y=Mean)) + 
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95), width=.1) +
    geom_point(shape = 21, color = "black", fill = "darkorange",
               size = 3) + 
    scale_x_discrete("") + 
    scale_y_continuous("", limits = c(-1,1)) +
    facet_wrap(~trait) +
    coord_flip() +
    theme_bw()
)

pol <- grid.arrange(p_scatter, p_rho)
ggsave("Figures/veg_genus.png", pol, width = 7, height = 7)
# Family
cwm <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "adjusted") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "family") %>% 
  filter(!zone == "BDM") # temporary filter till convergence issue is fixed
combo <- expand_grid(trait = unique(dfCWM_pol$trait), zone = unique(dfCWM_veg$zone)[-1])
rhos <- purrr::map2_dfr(combo$zone, combo$trait, ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y))
(p_scatter <- 
    ggplot(cwm, aes(x = pollen_mean, y = veg_mean, color = zone)) +
    geom_point() +
    scale_fill_manual() + 
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    facet_wrap(~trait, scales = "free") +
    ggtitle("Family level") + 
    theme_bw() + 
    theme(legend.position =  "bottom"))

(p_rho <- ggplot(rhos, aes(x = zone, y=Mean)) + 
    geom_errorbar(aes(ymin=Lower95, ymax=Upper95), width=.1) +
    geom_point(shape = 21, color = "black", fill = "darkorange",
               size = 3) + 
    scale_x_discrete("") + 
    scale_y_continuous("", limits = c(-1,1)) +
    facet_wrap(~trait) +
    coord_flip() +
    theme_bw()
)

pol <- grid.arrange(p_scatter, p_rho)
ggsave("Figures/veg-Family.png", pol, width = 7, height = 7)

