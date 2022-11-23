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

saveRDS(dfCWM_pol, "RDS_files/04_zCWM_estimates_pollen_Switzerland.rds")
saveRDS(dfCWM_veg, "RDS_files/04_zCWM_estimates_vegetation_Switzerland.rds")

selectedtrait <- "PlantHeight"

# Linear model ----
lm_bay <- function(cwm, selectedtrait, selecteddata){
  Mean <- cwm %>% 
    filter(trait == selectedtrait) %>% 
    filter(standardized == "no") %>% 
    pull(cwm_mean)
  zMean <- cwm %>% 
    filter(trait == selectedtrait) %>% 
    filter(standardized == "yes") %>% 
    pull(cwm_mean) 
  zSD <- cwm %>% 
    filter(trait == selectedtrait) %>% 
    filter(standardized == "yes") %>% 
    pull(cwm_sd) 
  Elev <- cwm %>% 
    filter(trait == selectedtrait) %>% 
    filter(standardized == "yes") %>% 
    pull(elevation)
  zElev <- (Elev - mean(Elev)) / sd(Elev) # standardize elevation
  N <- length(zMean)
  Habitat <- cwm %>% 
    filter(trait == selectedtrait) %>% 
    filter(standardized == "no") %>% 
    pull(habitat01) %>% 
    as.factor()
  # values to predict on
  Elev_pred <- seq(400, 2600, length.out = 250)
  zElev_pred <- rep((Elev_pred - mean(Elev_pred))/sd(Elev_pred), 2)
  Habitat_pred <- c(rep("grasslands", 250), rep("forests", 250)) %>% as.factor()
  Npred <- 500
  
  model_string <- "
      model {
      for(i in 1:N) {
      # regression on standardized data
        zMean[i] ~ dt(regression_fitted[i], 1/zsd^2, df)
        regression_fitted[i] <- zintercept +  zslope_elevation * zElev[i] + zhabitat_effect[Habitat[i]]
        regression_residual[i] <- zMean[i] - regression_fitted[i]
        zSD[i] ~ dnorm(zsd, tau_sd)
        }
    
      # priors
      zsd ~ dlnorm(0, 10^-4)
      tau_sd ~ dgamma(0.001, 0.001)
      zintercept ~ dnorm(0, 10^-4)
      zslope_elevation ~ dnorm(0,10^-4)
      df ~ dexp(1/30) # df of t distribution, Values of d) larger than about 30.0 make the t distribution roughly normal
      zhabitat_effect[1] <- 0
      zhabitat_effect[2] ~ dnorm(0, 10^-4)
  
      # predict
        for(i in 1:Npred){
        zmean_pred[i] ~ dt(zintercept +  zslope_elevation * zElev_pred[i] + zhabitat_effect[Habitat_pred[i]], 1/zsd^2, df)
        mean_pred[i] <- zmean_pred[i] * sd(Mean) + mean(Mean) # transform to original scale
        }
      #data# N, zMean, zSD, zElev, Mean, Habitat, zElev_pred, Habitat_pred, Npred
      #monitor# zintercept, zslope_elevation, zhabitat_effect[2], mean_pred, zsd
      #residual# regression_residual
      #fitted# regression_fitted
      #inits# zintercept, zslope_elevation, zhabitat_effect
    }"
  
  zintercept <- list(chain1 = min(zMean), chain2 = mean(zMean), chain3 = max(zMean))
  zslope_elevation <- list(chain1 = c(-10), chain2 = c(10), chain3 = 0)
  zhabitat_effect <- list(chain1 = c(NA, min(zMean)), chain2 = c(NA, mean(zMean)), chain3 = c(NA, max(zMean)))
  (results <- run.jags(model_string, n.chains = 3, 
                      modules = "glm on", sample = 10000))
  
  results_summary <- as.data.frame(summary(results)) %>% 
    mutate(trait = selectedtrait,
           data = selecteddata)
  if (!all(results_summary[,"SSeff"] >= 10000)) 
    warn("Effective sample size too low", x = results_summary)
 
  # save convergence diagnostics
  plot(results, vars = c("zintercept", "zslope_elevation", "zhabitat_effect[2]", "zsd"), 
       file = paste0("Convergence_diagnostics/04_Covergence_correlation_regression_", selectedtrait,
                              "_", selecteddata, ".pdf" ))
  # gelman plots
  pdf(paste0("Convergence_diagnostics/04_Gelman_correlation_regression_", selectedtrait,
             "_", selecteddata, ".pdf"))
  gelman.plot(results, vars = c("zintercept", "zslope_elevation", "zhabitat_effect[2]", "zsd"),  ask = FALSE)
  dev.off()
  
  # Create points to predict on
  lm_mcmc <- as.mcmc(results)
  pred_cwm_mean <- apply(lm_mcmc[, grep("mean_pred", colnames(lm_mcmc), value = FALSE)], MARGIN = 2, mean)
  # Calculate confidence interval
  credible_lower <- apply(lm_mcmc[, grep("mean_pred", colnames(lm_mcmc), value = FALSE)], MARGIN = 2, quantile, prob = 0.05)
  credible_upper <- apply(lm_mcmc[, grep("mean_pred", colnames(lm_mcmc), value = FALSE)], MARGIN = 2, quantile, prob = 0.95)
  df_pred <- tibble(elevation = rep(Elev_pred, 2),
                    habitat = Habitat_pred,
                    cwm_mean = pred_cwm_mean,
                    lower = credible_lower,
                    upper = credible_upper)
  
  # plot model results
  (p <- ggplot(df_pred, aes(x = elevation, y = cwm_mean, colour = habitat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "grey") +
    geom_line(color = "darkgreen") +
    # add data original data points
    geom_point(data = cwm[cwm$trait == selectedtrait & cwm$standardized == "no",], aes(x = elevation, y = cwm_mean)) +
    scale_y_continuous("") + 
    scale_x_continuous("Elevation (m)") +
    ggtitle( paste(labs_trait[selectedtrait], selecteddata)) +
    theme_bw())
  ggsave(paste0("Figures/CWM_elevation_", selectedtrait,"_", selecteddata, ".png"),p)
  
  # check model assumptions
  # Extract residuals and fitted estimates:
  residuals <- resid(results)
  fitted <- fitted(results)
  
  png(paste0("Figures/CWM_elevation_model_check_", selectedtrait,"_", selecteddata, ".png"))
  par(mfrow = c(2,1))
  # Plot of residuals indicates a potential problem:
  qqnorm(residuals); qqline(residuals)
  
  # Looking at residuals vs fitted indicates increasing variance:
  plot(fitted, residuals); abline(h = 0)
  dev.off()
}

labs_trait <- c("Height (log)(cm)",
                "Leaf area (log)(mm2)",
                "SLA (log) (mm2/mg)")
names(labs_trait) <- c("PlantHeight", "LA", "SLA")

cwm_veg <- dfCWM_veg %>% 
  filter(treatment == "Correction factor") 
lm_veg <- purrr::map(names(labs_trait), ~lm_bay(cwm = cwm_veg, selectedtrait = .x, selecteddata = "Vegetation"))

cwm_pol <- dfCWM_pol %>% 
  filter(treatment == "Correction factor") %>% 
  drop_na(elevation) %>% 
  filter(correction == "correction")
lm_pol <- purrr::map(names(labs_trait), ~lm_bay(cwm = cwm_pol, selectedtrait = .x, selecteddata = "Pollen"))


# Quick plots  ----
## Vegetation
### Elevation
dfCWM_veg %>% 
  filter(treatment == "Correction factor")  %>% 
  filter(standardized == "no") %>% 
ggplot(.,aes(x = elevation, y = cwm_mean, colour = habitat01)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  facet_wrap(.~trait, scales = "free") +
  theme_bw()
### MAT
dfCWM_veg %>% 
  filter(treatment == "Correction factor")  %>% 
  filter(standardized == "no") %>% 
ggplot(., aes(x = MAT, y = cwm_mean, colour = habitat01)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  facet_wrap(.~trait, scales = "free") +
  theme_bw()
### TAP
dfCWM_veg %>% 
  filter(treatment == "Correction factor")  %>% 
  filter(standardized == "no") %>% 
  ggplot(., aes(x = TAP, y = cwm_mean, colour = habitat01)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  facet_wrap(.~trait, scales = "free") +
  theme_bw()

## Pollen
### Elevation
dfCWM_pol %>% 
  filter(treatment == "Correction factor") %>% 
  filter(correction == "correction") %>% 
  filter(standardized == "no") %>% 
ggplot(.,aes(x = elevation, y = cwm_mean, color = habitat01)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  facet_wrap(.~trait, scales = "free") +
  theme_bw()

### MAT
dfCWM_pol %>% 
  filter(treatment == "Correction factor") %>% 
  filter(correction == "correction") %>% 
  filter(standardized == "no") %>% 
ggplot(.,aes(x = MAT, y = cwm_mean, color = habitat01)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  facet_wrap(.~trait, scales = "free") +
  theme_bw()

### TAP
dfCWM_pol %>% 
  filter(treatment == "Correction factor") %>% 
  filter(correction == "correction") %>% 
  filter(standardized == "no") %>% 
ggplot(.,aes(x = TAP, y = cwm_mean, color = habitat01)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  facet_wrap(.~trait, scales = "free") +
  theme_bw()


