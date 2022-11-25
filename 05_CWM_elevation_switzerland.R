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
bdm_meta <- read_csv("Data/bdm_metadata.csv") %>% 
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
                mutate(label = .y)) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(correction = case_when(str_detect(label, pattern = "adjustedpercent_draw") ~ str_extract(label, pattern = "draw[0-9]*$"),
                                str_detect(label, pattern = "adjustedpercent_mean") ~ "correction",
                                str_detect(label, pattern = "percent") ~ "no correction",
                                TRUE ~ "no correction"),
         pollination = replace_na(pollination, "allpol"),
         growthform = replace_na(growthform, "allpft"),
         sitename = str_remove(sitename, "X"),
         standardized = ifelse(str_detect(label, "zCWM"),"yes","no")
  ) %>% 
  dplyr::select(country, sitename, cwm_mean = Mean, 
                cwm_sd = SD, trait, correction, pollination, growthform, standardized) %>% 
  # create treatment column
  mutate(treatment = case_when((pollination == "allpol" & growthform == "allpft") ~ "Correction factor",
                               (pollination != "allpol" & growthform == "allpft") ~ "Pollination mode",
                               (growthform != "allpft" & pollination == "allpol") ~ "Growth form")) %>% 
  left_join(bdm_meta, by = "sitename") %>% 
  filter(!sitename == "unknown")

dfCWM_veg <- veg_files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(veg_lab, ~readRDS(.) %>% 
                as.data.frame() %>% 
                mutate(label = .y) %>% 
                mutate(across(one_of("sitename"),  as.character))) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(country = str_extract(label, pattern = "Scotland|Switzerland"),
         pollination = dplyr::recode(pollination, `no wind` = "not wind"),
         standardized = ifelse(str_detect(label, "zCWM"), "yes", "no")) %>% 
  dplyr::select(country, sitename, trait, cwm_mean = Mean, cwm_sd = SD,
                growthform, pollination, taxres, zone, standardized) %>% 
  # create treatment column
  mutate(treatment = case_when((pollination == "allpol" & growthform == "allpft" &
                                  taxres == "stand.spec") ~ "Correction factor",
                               taxres %in% c("family","genus") ~ "Taxonomic resolution",
                               (pollination != "allpol" & growthform == "allpft") ~ "Pollination mode",
                               (growthform != "allpft" & pollination == "allpol") ~ "Growth form"),
         zone = str_remove(zone, pattern = "Scotland |Switzerland ")) %>% 
  left_join(bdm_meta, by = "sitename")

# saveRDS(dfCWM_pol, "RDS_files/04_zCWM_estimates_pollen_Switzerland.rds")
# saveRDS(dfCWM_veg, "RDS_files/04_zCWM_estimates_vegetation_Switzerland.rds")

selectedtrait <- "PlantHeight"

# Linear model ----
lm_bay <- function(cwm, selectedhabitat,  selectedtrait, selecteddata){
  # Model data
  CWM <- cwm %>% 
    filter(habitat01 == selectedhabitat) %>% 
    filter(trait == selectedtrait) %>% 
    filter(standardized == "no") %>% 
    pull(cwm_mean)
  zCWM <- cwm %>% 
    filter(habitat01 == selectedhabitat) %>% 
    filter(trait == selectedtrait) %>% 
    filter(standardized == "yes") %>% 
    pull(cwm_mean) 
  zSD <- cwm %>% 
    filter(habitat01 == selectedhabitat) %>% 
    filter(trait == selectedtrait) %>% 
    filter(standardized == "yes") %>% 
    pull(cwm_sd) 
  Elev <- cwm %>% 
    filter(habitat01 == selectedhabitat) %>% 
    filter(trait == selectedtrait) %>% 
    filter(standardized == "yes") %>% 
    pull(elevation)
  zElev <- (Elev - mean(Elev)) / sd(Elev) # standardize elevation
  N <- length(zCWM)

  # values to predict on
  maxmin <- cwm %>%
    filter(habitat01 == selectedhabitat) %>% 
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
    warn(paste("Effective sample size potentially too low", selectedtrait, selecteddata, selectedhabitat))
  
  # check model assumptions
  # Extract residuals and fitted estimates:
  residuals <- resid(results)
  fitted <- fitted(results)
  png(paste0("Figures/CWM_elevation_model_check_", selectedtrait,"_", selecteddata,"_", selectedhabitat, ".png"))
  par(mfrow = c(2,1))
  # Plot of residuals indicates a potential problem:
  qqnorm(residuals); qqline(residuals)
  # Looking at residuals vs fitted indicates increasing variance:
  plot(fitted, residuals); abline(h = 0)
  dev.off()

  # save convergence diagnostics
  plot(results, vars = c("zintercept", "zslope_elevation", "zsd", "df"), 
       file = paste0("Convergence_diagnostics/04_Covergence_correlation_regression_", selectedtrait,
                              "_", selecteddata,"_",selectedhabitat, ".pdf" ))

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
                    habitat01 = selectedhabitat,
                    data = selecteddata,
                    cwm_mean = pred_cwm_mean,
                    credible_lower = credible_lower,
                    credible_upper = credible_upper,
                    uncertainty_lower = uncertainty_lower,
                    uncertainty_upper = uncertainty_upper,
                    trait = selectedtrait)
  saveRDS(df_pred, paste0("RDS_files/05_lm_switzerland_", selectedtrait, "_", selecteddata,"_",selectedhabitat, ".rds"))
  
 (p <- ggplot(df_pred, aes(x = elevation, y = cwm_mean)) +
    geom_ribbon(aes(ymin = credible_lower, ymax = credible_upper), alpha = 0.2, color = "grey") +
    geom_ribbon(aes(ymin = uncertainty_lower, ymax = uncertainty_upper), alpha = 0.2, color = "grey") +
    geom_line(color = "darkgreen") +
    geom_point(data = cwm[cwm$trait == selectedtrait & cwm$standardized == "no" & cwm$habitat01 == selectedhabitat,], aes(x = elevation, y = cwm_mean)) +
    scale_y_continuous("") + 
    scale_x_continuous("Elevation (m)") +
    theme_bw())
  ggsave(paste0("Figures/CWM_elevation_", selectedtrait,"_", selecteddata,"_", selectedhabitat, ".png"),p)
  }

traits <- dfCWM_pol$trait %>% unique()
cwm_veg <- dfCWM_veg %>% 
  filter(treatment == "Correction factor") %>% 
  mutate(trait = factor(trait, level = c("LA", "PlantHeight", "SLA"),
                        labels = c(LA = "Leaf~area~(mm^{2})(log)",
                                   PlantHeight = "Height~(cm)(log)",
                                   SLA = "SLA~(mm^{2}/mg)(log)"))) 
purrr::map(traits, ~lm_bay(cwm = cwm_veg, selectedtrait = .x, 
                                                selectedhabitat = "forests", selecteddata = "Vegetation"))
purrr::map(traits, ~lm_bay(cwm = cwm_veg, selectedtrait = .x, 
                                      selectedhabitat = "grasslands", selecteddata = "Vegetation"))

cwm_pol <- dfCWM_pol %>% 
  filter(treatment == "Correction factor") %>% 
  drop_na(elevation) %>% 
  filter(correction == "correction") %>% 
  mutate(trait = factor(trait, level = c("LA", "PlantHeight", "SLA"),
                        labels = c(LA = "Leaf~area~(mm^{2})(log)",
                                   PlantHeight = "Height~(cm)(log)",
                                   SLA = "SLA~(mm^{2}/mg)(log)"))) 
purrr::map(traits, ~lm_bay(cwm = cwm_pol, selectedtrait = .x, 
                                                selectedhabitat = "forests", selecteddata = "Pollen"))
purrr::map(traits, ~lm_bay(cwm = cwm_pol, selectedtrait = .x, 
                                      selectedhabitat = "grasslands", selecteddata = "Pollen"))


lm_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "05_lm")

lm_values <- lm_files %>% 
  folderpath.fun(.) %>% 
  purrr::map_df(~readRDS(.)) 

# plot model results
cwm_range <- cwm_veg %>%
  left_join(cwm_pol, by = c("sitename", "trait")) %>% 
  group_by(trait) %>%
  summarise(across(c("cwm_mean.x","cwm_mean.y"), ~range(., na.rm = T))) %>%
  transmute(max_lim = max(cwm_mean.x, cwm_mean.y, na.rm = T),
            min_lim = min(cwm_mean.x, cwm_mean.y, na.rm = T)) %>%
  pivot_longer(-trait, names_to = "limits", values_to = "cwm_mean") %>%
  mutate(elevation = rep(c(400,2600,400,2600)), habitat01 = NA) 
 

(p <- lm_values %>% 
    filter(data == "Vegetation" ) %>% 
    mutate(trait = factor(trait, level = c("LA", "PlantHeight", "SLA"),
                          labels = c(LA = "Leaf~area~(mm^{2})(log)",
                                            PlantHeight = "Height~(cm)(log)",
                                            SLA = "SLA~(mm^{2}/mg)(log)"))) %>% 
    ggplot(aes(x = elevation, y = cwm_mean, colour = habitat01)) +
    geom_ribbon(aes(ymin = uncertainty_lower, ymax = uncertainty_upper, fill = habitat01), 
                alpha = 0.2, colour = NA) +
    geom_ribbon(aes(ymin = credible_lower, ymax = credible_upper, fill = habitat01), 
                alpha = 0.2, colour = NA) +
    geom_line() +
    geom_point(data = cwm_veg[cwm_veg$standardized == "no",], aes(x = elevation, y = cwm_mean), size = 0.75) +
    geom_blank(data = cwm_range) +
    scale_color_manual(values = c(forests = "darkorange", grasslands = "purple"),
                       label = c(forests = "Forests", grasslands = "Grassland")) + 
    scale_fill_manual(values = c(forests = "darkorange", grasslands = "purple"),
                       label = c(forests = "Forest", grasslands = "Grassland")) + 
    scale_y_continuous("") + 
    scale_x_continuous("Elevation (m.s.l)", limits = c(400,2600)) +
    facet_wrap(~trait, scales = "free",
               labeller = label_parsed) + 
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()) +
   guides(fill = "none"))
ggsave("Figures/CWM_elevation_vegetaion.png", p, width = 7, height = 4)

(p <- lm_values %>% 
    filter(data == "Pollen" ) %>% 
    mutate(trait = factor(trait, level = c("LA", "PlantHeight", "SLA"),
                          labels = c(LA = "Leaf~area~(mm^{2})(log)",
                                     PlantHeight = "Height~(cm)(log)",
                                     SLA = "SLA~(mm^{2}/mg)(log)"))) %>% 
    ggplot(aes(x = elevation, y = cwm_mean, colour = habitat01)) +
    geom_ribbon(aes(ymin = uncertainty_lower, ymax = uncertainty_upper, fill = habitat01), 
                alpha = 0.2, colour = NA) +
    geom_ribbon(aes(ymin = credible_lower, ymax = credible_upper, fill = habitat01), 
                alpha = 0.2, colour = NA) +
    geom_line() +
    geom_point(data = cwm_pol[cwm_pol$standardized == "no",], aes(x = elevation, y = cwm_mean), size = 0.75) +
    geom_blank(data = cwm_range) +
    scale_color_manual(values = c(forests = "darkorange", grasslands = "purple"),
                       label = c(forests = "Forests", grasslands = "Grassland")) + 
    scale_fill_manual(values = c(forests = "darkorange", grasslands = "purple"),
                      label = c(forests = "Forest", grasslands = "Grassland")) + 
    scale_y_continuous("") + 
    scale_x_continuous("Elevation (m.s.l)", limits = c(400,2600)) +
    facet_wrap(~trait, scales = "free",
               labeller = label_parsed) + 
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()) +
    guides(fill = "none"))
ggsave("Figures/CWM_elevation_pollen.png", p, width = 7, height = 4)


