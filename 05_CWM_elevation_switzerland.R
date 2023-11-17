## ---------------------------
##
## Script name:  05_CWM_elevation_switzerland
##
## Purpose of script: Calculating CWM change along elevation gradient from pollen and vegetation
##
## Author: Annegreet Veeken
##
## Date Created: 2022-09-22
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   
## References:
## - Ch 17 Doing Bayesian analysis, Kruschke
## - cow_raising.R, Matt Denwood, PR statistics Bayesian course
## - http://biometry.github.io/APES//LectureNotes/StatsCafe/Linear_models_jags.html
## ---------------------------

set.seed(1)

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rjags)) install.packages("rjags")
if (!require(runjags)) install.packages("runjags")
if (!require(rlang)) install.packages("rlang")

## Load data
bdm_meta <- read.csv("Data/bdm_metadata.csv") %>% 
  mutate(sitename = as.character(aID_STAO))

# collate data
traits <- c("LA","SLA","PlantHeight")

pollen_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "03_CWM_estimates_|03_zCWM_estimates_") %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(pattern = "pollen") %>%  
  str_subset(pattern = "Switzerland") 

pollen_lab <- pollen_files %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_remove(., ".rds")

veg_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "03_CWM_estimates_|03_zCWM_estimates_") %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(pattern = "Switzerland") %>% 
  str_subset(pattern = "zoneA")

veg_lab <- veg_files %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_remove(., ".rds") 

folderpath.fun <- function(x)
{paste("RDS_files/", x, sep = "/")}

dfCWM_pol <- pollen_files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(pollen_lab, ~readRDS(.) %>% 
                mutate(label = .y) %>% 
                rownames_to_column("rowname")) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(correction = case_when(str_detect(label, pattern = "adjustedpercent_draw") ~ str_extract(label, pattern = "draw[0-9]*$"),
                                str_detect(label, pattern = "adjustedpercent_mean") ~ "correction",
                                str_detect(label, pattern = "percent") ~ "no correction",
                                TRUE ~ "no correction"),
         sitename = str_remove(sitename, "X"),
         pollination = replace_na(pollination, "allpol"),
         growthform = replace_na(growthform, "allpft"),
         standardized = ifelse(str_detect(label, "zCWM"),"yes",ifelse(str_detect(rowname, "zcwm"),"yes","no")),
         sp_list = ifelse(str_detect(label, "field"), "field", "GBIF")
  ) %>% 
  dplyr::select(country, sitename, cwm_mean = Mean, 
                cwm_sd = SD, trait, correction,  standardized, growthform, pollination, sp_list) %>% 
  left_join(bdm_meta, by = "sitename")
sitenames <- dfCWM_pol %>% pull(sitename) %>% unique

dfCWM_veg <- veg_files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(veg_lab, ~readRDS(.) %>% 
                as.data.frame() %>% 
                mutate(label = .y) %>% 
                mutate(across(one_of("sitename"),  as.character))) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(country = str_extract(label, pattern = "Scotland|Switzerland"),
         standardized = ifelse(str_detect(label, "zCWM"), "yes", "no")) %>% 
  dplyr::select(country, sitename, trait, cwm_mean = Mean, cwm_sd = SD,
                zone, standardized, growthform, pollination) %>% 
  # create treatment column
  mutate(zone = str_remove(zone, pattern = "Scotland |Switzerland ")) %>% 
  left_join(bdm_meta, by = "sitename")  %>% 
  # select sites in the pollen dataset
  filter(sitename %in% sitenames)

# saveRDS(dfCWM_pol, "RDS_files/04_zCWM_estimates_pollen_Switzerland.rds")
# saveRDS(dfCWM_veg, "RDS_files/04_zCWM_estimates_vegetation_Switzerland.rds")

selectedtrait <- "PlantHeight"

# Linear model ----
lm_bay <- function(cwm, selectedhabitat,  selectedpft, selectedpol, selectedtrait, selecteddata){
  cwm <- dfCWM_veg %>% 
    filter(habitat01 %in% selectedhabitat) %>% 
    filter(growthform == selectedpft) %>% 
    filter(pollination == selectedpol) %>% 
    filter(trait == selectedtrait) 
  
  # Model data
  CWM <- cwm %>% 
    filter(standardized == "no") %>% 
    pull(cwm_mean)
  zCWM <- cwm %>% 
    filter(standardized == "yes") %>% 
    pull(cwm_mean) 
  zSD <- cwm %>% 
    filter(standardized == "yes") %>% 
    pull(cwm_sd) 
  Elev <- cwm %>% 
    filter(standardized == "yes") %>% 
    pull(elevation)
  zElev <- (Elev - mean(Elev)) / sd(Elev) # standardize elevation
  N <- length(zCWM)

  # values to predict on
  maxmin <- cwm %>%
    filter(habitat01 %in% selectedhabitat) %>% 
    filter(growthform == selectedpft) %>%
    summarise(min = min(elevation), max = max(elevation))
  Elev_pred <- seq(from = as.numeric(maxmin[1,"min"]), to = as.numeric(maxmin[1,"max"]), by = 2)
  zElev_pred <- (Elev_pred - mean(Elev_pred)) / sd(Elev_pred)
  Npred <- length(zElev_pred)
  
  # define model
  model_string <- "model {
      for(i in 1:N) {
      # regression on standardized data
        zCWM[i] ~ dt(regression_fitted[i], 1/zsd^2, df)
        regression_fitted[i] <- zintercept +  zslope_elevation * zElev[i] 
        regression_residual[i] <- zCWM[i] - regression_fitted[i]
        zSD[i] ~ dnorm(zsd, tau_sd)
        }
    
      # priors
      zsd ~ dlnorm(0, 1/(1)^2)
      tau_sd ~ dgamma(0.001, 0.001)
      zintercept ~ dnorm(0, 1/(10)^2)
      zslope_elevation ~ dnorm(0, 1/(10)^2)
      df ~ dexp(1/30) # df of t distribution, Values of d) larger than about 30.0 make the t distribution roughly normal
  
      # transform to original scale
      slope_elevation <- zslope_elevation * sd(CWM)/sd(Elev)
      intercept <- zintercept * sd(CWM) + mean(CWM) - zslope_elevation*mean(Elev)*sd(CWM)/sd(Elev)
      
      # predict
      for(i in 1:Npred){
        zCWM_pred[i] ~ dt(zintercept +  zslope_elevation * zElev_pred[i], 1/zsd^2, df)
        CWM_pred[i] <- zCWM_pred[i] * sd(CWM) + mean(CWM) # transform to original scale
      }
  
      #data# N, zCWM, zSD, zElev, CWM, Elev, Npred, zElev_pred
      #monitor# zintercept, zslope_elevation, zsd, df, slope_elevation, intercept, CWM_pred
      #residual# regression_residual
      #fitted# regression_fitted
      #inits# zintercept, zslope_elevation, df
    }
  "
  # Inits
  zintercept <- list(chain1 = min(zCWM)*0.1, chain2 = mean(zCWM), chain3 = max(zCWM)*10)
  zslope_elevation <- list(chain1 = c(-5), chain2 = c(5), chain3 = 0)
  df <- list(chain1 = 1/1, chain2 = 1/15, chain3 = 1/30)
  
  # Run model
  (results <- run.jags(model_string, n.chains = 3, 
                        modules = "glm on", sample = 40000))
  
  # check for sufficient samples size
  results_summary <- as.data.frame(summary(results, vars = c("zintercept", "zslope_elevation", "zsd", "df"))) %>% 
    mutate(trait = selectedtrait,
           data = selecteddata)
  if (!all(results_summary[,"SSeff"] >= 10000))
    warn(paste("Effective sample size potentially too low", selectedtrait, selecteddata, paste(selectedhabitat, collapse = "_")))
  
  # check model assumptions
  # Extract residuals and fitted estimates:
  residuals <- resid(results)
  fitted <- fitted(results)
  png(paste0("Figures/CWM_elevation_model_check_", selectedtrait,"_", selecteddata,"_", paste(selectedhabitat, collapse = "_"),
             "_",selectedpft,"_",selectedpol, ".png"))
  par(mfrow = c(2,1))
  # Plot of residuals indicates a potential problem:
  qqnorm(residuals); qqline(residuals)
  # Looking at residuals vs fitted indicates increasing variance:
  plot(fitted, residuals); abline(h = 0)
  dev.off()

  # save convergence diagnostics
  plot(results, vars = c("zintercept", "zslope_elevation", "zsd", "df"), 
       file = paste0("Convergence_diagnostics/04_Covergence_lm_", selectedtrait,
                              "_", selecteddata,"_",paste(selectedhabitat, collapse = "_"),"_",selectedpft,"_",selectedpol, ".pdf" ))

  # plot model results
  lm_mcmc <- as.mcmc(results)
  # credibility around model parameters
  pred_cwm_mean <- mean(lm_mcmc[,"intercept"]) + Elev_pred * mean(lm_mcmc[,"slope_elevation"])
  pred_mean_dist <- matrix(NA, nrow = nrow(lm_mcmc), ncol = length(Elev_pred))
  for (i in 1:nrow(pred_mean_dist)) {
    pred_mean_dist[i,] <- lm_mcmc[i,"intercept"] + Elev_pred * lm_mcmc[i,"slope_elevation"]
  }
  credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.05)
  credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.95)
  # uncertainty in prediction
  uncertainty_lower <- apply(lm_mcmc[, grep("CWM_pred", colnames(lm_mcmc), value = FALSE)], MARGIN = 2, quantile, prob = 0.05)
  uncertainty_upper <- apply(lm_mcmc[, grep("CWM_pred", colnames(lm_mcmc), value = FALSE)], MARGIN = 2, quantile, prob = 0.95)

  df_pred <- tibble(elevation = Elev_pred,
                    habitat01 = paste(selectedhabitat, collapse = "_"),
                    growthform = selectedpft,
                    pollination = selectedpol,
                    data = selecteddata,
                    cwm_mean = pred_cwm_mean,
                    credible_lower = credible_lower,
                    credible_upper = credible_upper,
                    uncertainty_lower = uncertainty_lower,
                    uncertainty_upper = uncertainty_upper,
                    trait = selectedtrait)
  saveRDS(df_pred, paste0("RDS_files/05_lm_switzerland_", selectedtrait, "_", selecteddata,"_",paste(selectedhabitat, collapse = "_"),
                          "_",selectedpft,"_",selectedpol,".rds"))

  }

traits <- c("PlantHeight", "LA", "SLA")
## vegetation 
cwm_veg <- dfCWM_veg %>% 
  filter(growthform == "allpft")
purrr::map(traits, ~lm_bay(cwm = cwm_veg,
                           selectedtrait = .x,
                           selectedhabitat = "forests",
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selecteddata = "Vegetation"))
purrr::map(traits, ~lm_bay(cwm = cwm_veg,
                           selectedtrait = .x,
                           selectedhabitat = "grasslands",
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selecteddata = "Vegetation"))
purrr::map(traits, ~lm_bay(cwm = cwm_veg,
                           selectedtrait = .x,
                           selectedhabitat = c("forests","grasslands"),
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selecteddata = "Vegetation"))
# herb
purrr::map(traits, ~lm_bay(cwm = dfCWM_veg, 
                           selectedtrait = .x, 
                           selectedhabitat = c("forests","grasslands"), 
                           selectedpft = "herb",
                           selectedpol = "allpol",
                           selecteddata = "Vegetation"))
# trsh
purrr::map(traits, ~lm_bay(cwm = dfCWM_veg, 
                           selectedtrait = .x, 
                           selectedhabitat = c("forests","grasslands"), 
                           selectedpft = "trsh",
                           selectedpol = "allpol",
                           selecteddata = "Vegetation"))
# not wind pollinated
purrr::map(traits, ~lm_bay(cwm = dfCWM_veg, 
                           selectedtrait = .x, 
                           selectedhabitat = c("forests","grasslands"), 
                           selectedpft = "allpft",
                           selectedpol = "no wind",
                           selecteddata = "Vegetation"))
# wind pollinated
purrr::map(traits, ~lm_bay(cwm = dfCWM_veg, 
                           selectedtrait = .x, 
                           selectedhabitat = c("forests","grasslands"), 
                           selectedpft = "allpft",
                           selectedpol = "wind",
                           selecteddata = "Vegetation"))
# No subset
purrr::map(traits, ~lm_bay(cwm = dfCWM_veg,
                           selectedtrait = .x,
                           selectedhabitat = c("forests","grasslands"),
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selecteddata = "Vegetation"))

## pollen
cwm_pol <- dfCWM_pol %>% 
  drop_na(elevation) %>% 
  filter(correction == "no correction") %>% 
  filter(sp_list == "GBIF")
# habitat
purrr::map(traits, ~lm_bay(cwm = cwm_pol,
                           selectedtrait = .x,
                           selectedhabitat = "forests",
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selecteddata = "Pollen"))
purrr::map(traits, ~lm_bay(cwm = cwm_pol,
                           selectedtrait = .x,
                           selectedhabitat = "grasslands",
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selecteddata = "Pollen"))
# PFT
purrr::map(traits, ~lm_bay(cwm = cwm_pol, 
                           selectedtrait = .x, 
                           selectedpft = "herb",
                           selectedpol = "allpol",
                           selectedhabitat = c("forests", "grasslands"), 
                           selecteddata = "Pollen"))
purrr::map(traits, ~lm_bay(cwm = cwm_pol, 
                           selectedtrait = .x, 
                           selectedpft = "trsh",
                           selectedpol = "allpol",
                           selectedhabitat = c("forests", "grasslands"), 
                           selecteddata = "Pollen"))
# pollination mode
purrr::map(traits, ~lm_bay(cwm = cwm_pol, 
                           selectedtrait = .x, 
                           selectedpft = "allpft",
                           selectedpol = "wind",
                           selectedhabitat = c("forests", "grasslands"), 
                           selecteddata = "Pollen"))
purrr::map(traits, ~lm_bay(cwm = cwm_pol, 
                           selectedtrait = .x, 
                           selectedpft = "allpft",
                           selectedpol = "not wind",
                           selectedhabitat = c("forests", "grasslands"), 
                           selecteddata = "Pollen"))
# No subset
purrr::map(traits, ~lm_bay(cwm = cwm_pol, 
                           selectedtrait = .x, 
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selectedhabitat = c("forests", "grasslands"), 
                           selecteddata = "Pollen"))

# CWM with sp list based on field observations
cwm_pol_field <- dfCWM_pol %>%
  drop_na(elevation) %>% 
  filter(correction == "no correction") %>% 
  filter(sp_list == "field")
purrr::map(traits, ~lm_bay(cwm = cwm_pol_field,
                           selectedtrait = .x,
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selectedhabitat = c("forests", "grasslands"),
                           selecteddata = "Pollen_fieldsp"))
# corrected
cwm_pol_correction <- dfCWM_pol %>%
  drop_na(elevation) %>% 
  filter(correction == "correction") 
purrr::map(traits, ~lm_bay(cwm = cwm_pol_correction,
                           selectedtrait = .x,
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selectedhabitat = c("forests", "grasslands"),
                           selecteddata = "Pollen_correction"))
purrr::map(traits, ~lm_bay(cwm = cwm_pol_correction,
                           selectedtrait = .x,
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selectedhabitat = c("forests"),
                           selecteddata = "Pollen_correction"))
purrr::map(traits, ~lm_bay(cwm = cwm_pol_correction,
                           selectedtrait = .x,
                           selectedpft = "allpft",
                           selectedpol = "allpol",
                           selectedhabitat = "grasslands",
                           selecteddata = "Pollen_correction"))
purrr::map(traits, ~lm_bay(cwm = cwm_pol_correction,
                           selectedtrait = .x,
                           selectedpft = "trsh",
                           selectedpol = "allpol",
                           selectedhabitat = c("forests", "grasslands"),
                           selecteddata = "Pollen_correction"))
purrr::map(traits, ~lm_bay(cwm = cwm_pol_correction,
                           selectedtrait = .x,
                           selectedpft = "herb",
                           selectedpol = "allpol",
                           selectedhabitat = c("forests", "grasslands"),
                           selecteddata = "Pollen_correction"))
purrr::map(traits, ~lm_bay(cwm = cwm_pol_correction,
                           selectedtrait = .x,
                           selectedpft = "allpft",
                           selectedpol = "wind",
                           selectedhabitat = c("forests", "grasslands"),
                           selecteddata = "Pollen_correction"))
purrr::map(traits, ~lm_bay(cwm = cwm_pol_correction,
                           selectedtrait = .x,
                           selectedpft = "allpft",
                           selectedpol = "not wind",
                           selectedhabitat = c("forests", "grasslands"),
                           selecteddata = "Pollen_correction"))


# Plotting ----
# load in lm results
lm_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "05_lm") 

lm_values <- lm_files %>% 
  folderpath.fun(.) %>% 
  purrr::map2_df(.,lm_files, ~readRDS(.) %>% 
                   mutate(f = .y)) 

# plot model results
cwm_range <- lm_values %>% 
  group_by(trait) %>%
  summarise(across(c("uncertainty_lower","uncertainty_upper"), ~range(., na.rm = T))) %>%
  transmute(max_lim = max(uncertainty_upper, na.rm = T),
            min_lim = min(uncertainty_lower, na.rm = T)) %>%
  pivot_longer(-trait, names_to = "limits", values_to = "cwm_mean") %>%
  mutate(elevation = rep(c(400,2600,400,2600)), habitat01 = NA, growthform = NA) 
labs_trait <- as_labeller(c(PlantHeight = "Height~(cm)(log)",
                            LA = "Leaf~area~(cm^{2})(log)",
                            SLA ="SLA~(mm^{2}/mg)(log)"),
                          default = label_parsed)

cwm_points_veg <- dfCWM_veg[dfCWM_veg$standardized == "no" & dfCWM_veg$growthform %in% c("trsh","herb"),] %>% 
  select(elevation, growthform,trait, cwm_mean) %>% 
  mutate(data = "Vegetation")
cwm_points <- cwm_pol[cwm_pol$standardized == "no" & cwm_pol$growthform %in% c("trsh","herb"),] %>% 
  select(elevation, growthform, trait, cwm_mean) %>% 
  mutate(data = "Pollen") %>% 
  bind_rows(cwm_points_veg) %>% 
  mutate(data =  factor(data, levels = c("Vegetation","Pollen", "Pollen_corrected"))) 
(p <- lm_values %>% 
    filter(habitat01 == "forests_grasslands") %>% 
    filter(growthform == "herb"| growthform == "trsh") %>%
    mutate(data =  factor(data, levels = c("Vegetation","Pollen", "Pollen_corrected"))) %>% 
    ggplot(aes(x = elevation, y = cwm_mean, colour = growthform)) +
    geom_ribbon(aes(ymin = uncertainty_lower, ymax = uncertainty_upper, fill = growthform), 
                alpha = 0.2, colour = NA) +
    geom_ribbon(aes(ymin = credible_lower, ymax = credible_upper, fill = growthform), 
                alpha = 0.2, colour = NA) +
    geom_line() +
    geom_point(data = cwm_points,
               aes(x = elevation, y = cwm_mean), size = 0.75) +
    geom_blank(data = cwm_range) +
    scale_color_manual(values = c(trsh = "darkorange", herb = "purple"),
                       label = c(trsh = "Woody", herb = "Non-woody")) + 
    scale_fill_manual(values = c(trsh = "darkorange", herb = "purple"),
                      label = c(trsh = "Woody", herb = "Non-woody")) + 
    scale_y_continuous("") + 
    scale_x_continuous("Elevation (m.s.l)", limits = c(400,2600)) +
    facet_grid(trait~data, scales = "free",
               labeller = labeller(trait = labs_trait)) + 
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()) +
    guides(fill = "none"))
ggsave("Figures/CWM_elevation_veg_pol.png", height = 7, width = 7)
(p <- lm_values %>% 
    filter(data == "Vegetation"|data == "Pollen") %>% 
    filter(!habitat01 == "forests_grasslands") %>% 
    filter(growthform %in% c("herb", "trsh")) %>% 
    ggplot(aes(x = elevation, y = cwm_mean, colour = growthform)) +
    geom_ribbon(aes(ymin = uncertainty_lower, ymax = uncertainty_upper, fill = growthform), 
                alpha = 0.2, colour = NA) +
    geom_ribbon(aes(ymin = credible_lower, ymax = credible_upper, fill = growthform), 
                alpha = 0.2, colour = NA) +
    geom_line() +
    # geom_point(data = dfCWM_veg[dfCWM_veg$standardized == "no" & 
    #                               dfCWM_veg$growthform %in% c("trsh","herb"),], 
    #            aes(x = elevation, y = cwm_mean), size = 0.75) +
    geom_blank(data = cwm_range) +
    scale_color_manual(values = c(trsh = "darkorange", herb = "purple"),
                       label = c(trsh = "Woody", herb = "Non-woody")) + 
    scale_fill_manual(values = c(trsh = "darkorange", herb = "purple"),
                      label = c(trsh = "Woody", herb = "Non-woody")) + 
    scale_y_continuous("") + 
    scale_x_continuous("Elevation (m.s.l)", limits = c(400,2600)) +
    facet_grid(trait~data, scales = "free",
               labeller = labeller(trait = labs_trait)) + 
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()) +
    guides(fill = "none"))
ggsave("Figures/CWM_elevation_vegetation_pft.png", p, width = 7, height = 4)

(p <- lm_values %>% 
    filter(data == "Pollen" ) %>% 
    filter(growthform %in% c("herb", "trsh")) %>% 
    ggplot(aes(x = elevation, y = cwm_mean, colour = growthform)) +
    geom_ribbon(aes(ymin = uncertainty_lower, ymax = uncertainty_upper, fill = growthform), 
                alpha = 0.2, colour = NA) +
    geom_ribbon(aes(ymin = credible_lower, ymax = credible_upper, fill = growthform), 
                alpha = 0.2, colour = NA) +
    geom_line() +
    geom_point(data = cwm_pol[cwm_pol$standardized == "no" & 
                                cwm_pol$growthform %in% c("trsh","herb"),], 
               aes(x = elevation, y = cwm_mean), size = 0.75) +
    geom_blank(data = cwm_range) +
    scale_color_manual(values = c(trsh = "darkorange", herb = "purple"),
                       label = c(trsh = "Woody", herb = "Non-woody")) + 
    scale_fill_manual(values = c(trsh = "darkorange", herb = "purple"),
                      label = c(trsh = "Woody", herb = "Non-woody")) + 
    scale_y_continuous("") + 
    scale_x_continuous("Elevation (m.s.l)", limits = c(400,2600)) +
    facet_wrap(data~trait, scales = "free",
               labeller = labs_trait
    ) + 
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()) +
    guides(fill = "none"))
ggsave("Figures/CWM_elevation_pollen_pft.png", p, width = 7, height = 4)

# By habitat 
(p <- lm_values %>% 
    # filter(habitat01 %in% c("forests", "grasslands")) %>% 
    filter(data  %in% c("Pollen", "Vegetation", "Pollen_correction")) %>% 
    filter(pollination == "allpol", growthform == "allpft") %>% 
    ggplot(aes(x = elevation, y = cwm_mean, colour = habitat01)) +
    geom_ribbon(aes(ymin = uncertainty_lower, ymax = uncertainty_upper, fill = habitat01), 
                alpha = 0.2, colour = NA) +
    geom_ribbon(aes(ymin = credible_lower, ymax = credible_upper, fill = habitat01), 
                alpha = 0.2, colour = NA) +
    geom_line() +
    geom_blank(data = cwm_range) +
    # scale_color_manual(values = c(forests = "darkorange", grasslands = "purple"),
    #                    label = c(forests = "Forests", grasslands = "Grassland")) + 
    # scale_fill_manual(values = c(forests = "darkorange", grasslands = "purple"),
    #                   label = c(forests = "Forest", grasslands = "Grassland")) + 
    scale_y_continuous("") + 
    scale_x_continuous("Elevation (m.s.l)", limits = c(400,2600)) +
    facet_grid(trait~data, scales = "free"
               # labeller = labs_trait
               ) + 
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()) +
    guides(fill = "none"))
ggsave("Figures/CWM_elevation_pollen_fieldsp.png", p, width = 7, height = 4)

