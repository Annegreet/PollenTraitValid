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
if(!require(gridExtra)) install.packages("gridExtra")

# load PFT summary table
# filter sites with to little vegetation cover of particular type
sum <- readRDS("RDS_files/01_Percentage_cover_pft.rds")
dropsites <- sum %>% 
  transmute(drop_herb = case_when(`Non-Woody`<=5 ~ paste(zone, sitename)),
         drop_trsh = case_when(Woody <= 5 ~ paste(zone,sitename)),
         drop_wind = case_when(wind <= 5 ~ paste(zone,sitename)),
         drop_nowind = case_when(`not wind` <= 5 ~ paste(zone,sitename)))

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
                mutate(label = .y)) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(correction = case_when(str_detect(label, pattern = "adjustedpercent_draw") ~ "correction draw",
                                str_detect(label, pattern = "adjustedpercent_mean") ~ "correction",
                                str_detect(label, pattern = "percent") ~ "no correction",
                                TRUE ~ "correction"),
         pollination = replace_na(pollination, "allpol"),
         growthform = replace_na(growthform, "allpft"),
         sitename = str_remove(sitename, "X")
         ) %>% 
  dplyr::select(country, sitename, pollen_mean = Mean, 
                pollen_sd = SD, trait, correction, pollination, growthform) 

dfCWM_veg <- veg_files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(veg_lab, ~readRDS(.) %>% 
                as.data.frame() %>% 
                rownames_to_column("parameter") %>% 
                filter(., str_detect(parameter, "cwm\\[")) %>%
                mutate(label = .y) %>% 
                mutate(across(one_of("sitename"),  as.character))) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(country = str_extract(label, pattern = "Scotland|Switzerland"),
         zone = paste(country, zone),
         pollination = recode(pollination, `no wind` = "not wind")) %>% 
  dplyr::select(country, sitename, trait, veg_mean = Mean, veg_sd = SD,
                growthform, pollination, taxres, zone) %>% 
  # remove sites with to little cover of pft (5% or less) 
  # see sites in df dropsites
  filter(!(growthform == "trsh" & zone == "zoneA"), # Too little sites with tree cover in zone A
         !(growthform == "trsh" & zone == "zoneB" & sitename %in% c("C005", "C013")),
         !(pollination == "wind" & zone == "zoneA" & sitename %in% c("665174", "C009","C015")),
         !(pollination == "not wind" & zone == "zoneA" & sitename %in% c("563166")))
#Save compiled data
saveRDS(dfCWM_pol, "RDS_files/04_CWM_estimates_pollen.rds") 
saveRDS(dfCWM_veg, "RDS_files/04_CWM_estimates_vegetation.rds") 

 
selectedzone <- combo$zone[1]
selectedtrait <- combo$trait[1]

cor_bay <- function(cwm, selectedzone, selectedtrait){
  CWM <- cwm %>% 
    filter(zone == selectedzone,
           trait == selectedtrait)
  CWM_mean <- CWM %>% 
    ungroup() %>% 
    dplyr::select(pollen_mean, veg_mean) %>% 
    as.matrix() 
  CWM_sd <- CWM %>%
    ungroup() %>%
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

trait_lab <-  c("Leaf area", "Plant height", 
                "Specific leaf area")
names(trait_lab) <- c("LA", "PlantHeight","SLA")

## No correction + correction----
# No correction
cwm <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "no correction") %>% 
  left_join(dfCWM_veg, by = c("country", "sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") %>% 
  group_by(trait, zone) %>% 
  mutate(nobs = n())
combo <- cwm %>% 
  dplyr::select(trait, zone) %>% 
  distinct()
rhos <- purrr::map2_dfr(combo$zone, combo$trait, 
                        ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y)) %>% 
  mutate(treatment = "No correction")

# Pollen correction 
cwm2 <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "correction") %>% 
  left_join(dfCWM_veg, by = c("country","sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") %>% 
  group_by(trait, zone) %>% 
  mutate(nobs = n())
combo2 <- cwm2 %>% 
  dplyr::select(trait, zone) %>% 
  distinct()
rhos2 <- purrr::map2_dfr(combo2$zone, combo2$trait, 
                         ~cor_bay(cwm2, selectedzone = .x, selectedtrait = .y)) %>% 
  mutate(treatment = "Correction") %>% 
  bind_rows(rhos)

# (p_scatter <- 
#     ggplot(cwm, aes(x = pollen_mean, y = veg_mean, color = zone)) +
#     geom_point() +
#     scale_fill_manual() + 
#     scale_x_continuous("Pollen") +
#     scale_y_continuous("Vegetation") +
#     facet_wrap(~trait, scales = "free") +
#     ggtitle("Corrected pollen data") + 
#     theme_bw() + 
#     theme(legend.position =  "bottom"))

(p_rho <- ggplot(rhos2, aes(x = zone, y = Mean, fill = treatment)) + 
    # geoms
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .1,
                  position = position_dodge(0.5)) +
    geom_point(shape = 21, size = 3, position = position_dodge(0.5)) + 
    # labels and scales
    scale_y_continuous("", limits = c(-1,1)) +
    scale_fill_manual(values = c("darkorange", "purple")) +
    labs(fill = "Treatment") +
    coord_flip() +
    # Faceting
    facet_wrap(~trait, labeller = labeller(trait = trait_lab)) +
    # Theme
    theme_bw()
)


ggsave("Figures/Pollen_corrected.png", p_rho, width = 7, height = 7)

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
         correction == "correction") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec")  %>% 
  group_by(trait, zone) %>% 
  mutate(nobs = n())
combo <- cwm %>% 
  dplyr::select(trait, zone) %>% 
  distinct()
rhos <- purrr::map2_dfr(combo$zone, combo$trait, 
                        ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y)) %>% 
  mutate(treatment = "Wind") 

# No wind
cwm2 <- dfCWM_pol %>%
  filter(growthform == "allpft",
         pollination == "not wind",
         correction == "correction")
cwm2 <- cwm2  %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>%
  filter(taxres == "stand.spec")
combo2 <- cwm2 %>% 
  dplyr::select(trait, zone) %>% 
  distinct()
rhos2 <- purrr::map2_dfr(combo2$zone, combo2$trait, 
                         ~cor_bay(cwm2, selectedzone = .x, selectedtrait = .y)) %>% 
  mutate(treatment = "Not wind") %>% 
  bind_rows(rhos)

(p_rho <- ggplot(rhos2, aes(x = zone, y = Mean, fill = treatment)) + 
    # geoms
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .1,
                  position = position_dodge(0.5)) +
    geom_point(shape = 21, size = 3, position = position_dodge(0.5)) + 
    # labels and scales
    scale_y_continuous("", limits = c(-1,1)) +
    scale_fill_manual(values = c("darkorange", "purple")) +
    labs(fill = "Treatment") +
    coord_flip() +
    # Faceting
    facet_wrap(~trait, labeller = labeller(trait = trait_lab)) +
    # Theme
    theme_bw()
)

pol <- grid.arrange(p_scatter, p_rho)
ggsave("Figures/Pollination_mode.png", p_rho, width = 7, height = 7)

## Plant functional type ----
# Trees and Shrubs
cwm <- dfCWM_pol %>% 
  filter(growthform == "trsh",
         pollination == "allpol",
         correction == "correction") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 
combo <- cwm %>% 
  dplyr::select(trait, zone) %>% 
  distinct()
rhos <- purrr::map2_dfr(combo$zone, combo$trait, 
                        ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y)) %>% 
  mutate(treatment = "Woody vegetation")

# Herbs
cwm2 <- dfCWM_pol %>% 
  filter(growthform == "herb",
         pollination == "allpol",
         correction == "correction") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 
combo2 <- cwm2 %>% 
  dplyr::select(trait, zone) %>% 
  distinct()
rhos2 <- purrr::map2_dfr(combo2$zone, combo2$trait, 
                         ~cor_bay(cwm2, selectedzone = .x, selectedtrait = .y)) %>% 
  mutate(treatment = "Non-woody vegetation")

(p_rho <- ggplot(rhos2, aes(x = zone, y = Mean, fill = treatment)) + 
    # geoms
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = .1,
                  position = position_dodge(0.5)) +
    geom_point(shape = 21, size = 3, position = position_dodge(0.5)) + 
    # labels and scales
    scale_y_continuous("", limits = c(-1,1)) +
    scale_fill_manual(values = c("darkorange", "purple")) +
    labs(fill = "Treatment") +
    coord_flip() +
    # Faceting
    facet_wrap(~trait, labeller = labeller(trait = trait_lab)) +
    # Theme
    theme_bw()
)

ggsave("Figures/PFT.png", p_rho, width = 7, height = 7)
## Taxonomic resolution ----
# Genus
cwm <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "correction") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "genus")
combo <- cwm %>% 
  dplyr::select(trait, zone) %>% 
  distinct()
rhos <- purrr::map2_dfr(combo$zone, combo$trait, 
                        ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y)) %>% 
  mutate(treatment = "Genus")

# Family
cwm2 <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "correction") %>% 
  left_join(dfCWM_veg, by = c("sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "family") 
combo2 <- cwm2 %>% 
  dplyr::select(trait, zone) %>% 
  distinct()
rhos2 <- purrr::map2_dfr(combo2$zone, combo2$trait,
                         ~cor_bay(cwm2, selectedzone = .x, selectedtrait = .y)) %>% 
  mutate(treatment = "Family") %>% 
  bind_rows(rhos)
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

ggsave("Figures/Taxonomic_resolution.png", p_rho, width = 7, height = 7)

