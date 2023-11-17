## ---------------------------
##
## Script name: SI_reduced_splist
##
## Purpose of script: Testing the effect of reducing the species list (using the data collected in the field) on the CWM reconstruction
##
## Author: Annegreet Veeken
##
## Date Created: 2022-11-03
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   
## References:
##
## ---------------------------

set.seed(1)
options(scipen = 9999)

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rjags))  install.packages("rjags")
if (!require(runjags)) install.packages("runjags")
if (!require(mcmc)) install.packages("mcmc")
if (!require(coda)) install.packages("coda")
if (!require(furrr))  install.packages("furrr")

source("Cwm_functions.R")

## Load and prepare data -----
# Trait data ----
dfTRAIT <- readRDS("RDS_files/02_Gapfilled_traits_pollen.rds")

# Pollen data ----
dfPOL_Scot <- readRDS("RDS_files/01_pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_pollen_data_Swiss.rds")
dfPOL <- bind_rows("Scotland" = dfPOL_Scot, "Switzerland" =
                     dfPOL_Swiss, .id = "country") %>% 
  ungroup()

# Species list vegetation 
sp_a <- readRDS("RDS_files/01_Species_abundance_a.rds") %>% pull(stand.spec) %>% unique()
sp_b <- readRDS("RDS_files/01_Species_abundance_bc.rds") %>% pull(stand.spec) %>% unique()
sp_scot <- c(sp_a, sp_b) %>% unique()
sp_swiss <- readRDS("RDS_files/01_Species_abundance_Swiss.rds") %>% pull(stand.spec) %>% unique()

selectedcountry <- "Scotland"
selectedtrait <- "LA"
selectedabun <- str_subset(colnames(dfPOL), "percent")[1]
normaltraits <- c("PlantHeight","SLA","LA")

# Pollen  ----
labs_trait <- as_labeller(c(PlantHeight = "Height~(cm)(log)",
                            LA = "Leaf~area~(cm^{2})(log)",
                            SLA ="SLA~(mm^{2}/mg)(log)"),
                          default = label_parsed)

# Run model
plan(multisession(workers = 2))
furrr::future_map(normaltraits,
                    ~cwm_pol_field(selectedtrait = .x,
                                selectedcountry = "Scotland",
                                selectedabun = "percent",
                                taxtrait = TRUE,
                                fieldspec = sp_scot),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(normaltraits,
                  ~cwm_pol_field(selectedtrait = .x,
                                 selectedcountry = "Switzerland",
                                 selectedabun = "percent",
                                 taxtrait = TRUE,
                                 fieldspec = sp_swiss),
                  .options = furrr_options(seed = TRUE))
# Scotland ----
# Compile estimates
pollen_files <- 
  list.files("RDS_files/",full.names = TRUE) %>% 
  str_subset(pattern = "pollen") %>% 
  str_subset(pattern = "field") %>% 
  str_subset(pattern = "Scotland")

pollen_files_all <-  c("RDS_files/03_CWM_estimates_pollen_Scotland_LA_percent.rds", 
                       "RDS_files/03_CWM_estimates_pollen_Scotland_SLA_percent.rds",
                       "RDS_files/03_CWM_estimates_pollen_Scotland_PlantHeight_percent.rds")

dfCWM_pol_field <-  purrr::map(pollen_files, ~readRDS(.x) %>% 
                as.data.frame())  %>% 
  bind_rows()  %>% 
  rownames_to_column(var = "standardized") %>% 
  filter(str_detect(standardized, pattern = "zcwm"))  %>% 
  dplyr::select(country, sitename, pollen_mean = Mean, 
                pollen_sd = SD, trait)

labs_trait <- as_labeller(c(PlantHeight = "Plant~height~(cm)(log)",
                            LA = "Leaf~area~(cm^{2})(log)",
                            SLA ="SLA~(mm^{2}/mg)(log)"),
                          default = label_parsed)
  
cor_bay <- function(cwm, selectedzone, selectedtrait, selectedtreatment){
  CWM <- cwm %>% 
    filter(zone == selectedzone,
           trait == selectedtrait)
  pollen_cwm <- CWM %>% 
    pull(pollen_mean) 
  veg_cwm <- CWM %>% 
    pull(veg_mean) 
  pollen_sd <- CWM %>% 
    pull(pollen_sd) 
  veg_sd <- CWM %>% 
    pull(veg_sd) 
  N <- length(pollen_cwm)
  
  model <- "model {
    for(i in 1:N) {
      pollen_cwm[i] ~ dt(regression_fitted[i], 1/pollen_sd[i]^2, df)
      regression_fitted[i] <- 0 + slope * veg_cwm[i] 
    }
    # Priors
    slope ~ dunif(-1,1)
    df ~ dexp(1/30) # df of t distribution, Values of d) larger than about 30.0 make the t distribution roughly normal
  }
  #monitor# slope, df
  #data# pollen_cwm, veg_cwm, pollen_sd, N 
  #inits# slope, df
  "
  slope <- list(chain1 = -1, chain2 = 0, chain3 = 1)
  df <- list(chain1 = 1, chain2 = 15, chain3 = 30)
  
  # run model 
  (results <- run.jags(model, n.chains = 3, sample = 30000,
                       modules = "glm on"))
  
  # save convergence diagnostics
  plot(results, file = paste0("Convergence_diagnostics/04_Covergence_regression_",selectedtrait,
                              "_", selectedzone,"_",selectedtreatment, "field.pdf" ))
  # gelman plots
  pdf(paste0("Convergence_diagnostics/04_Gelman_regression_", selectedtrait,
             "_", selectedzone,"_",selectedtreatment, "field.pdf"))
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  results_summary <- as.data.frame(summary(results)) %>% 
    rownames_to_column("parameter") %>% 
    mutate(trait = selectedtrait,
           zone = selectedzone,
           treatment = selectedtreatment,
           cor_nonbay = cor(pollen_cwm, veg_cwm))
  
  saveRDS(results_summary, paste0("RDS_files/04_slope_regression_", selectedtrait,
                                  "_", selectedzone, "_", selectedtreatment, "field.rds"))
  return(results_summary)
  
}

dfCWM_veg <- readRDS("RDS_files/04_zCWM_estimates_vegetation_Scotland.rds") %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         taxres == "stand.spec")

cwm <- dfCWM_pol_field %>%
  left_join(dfCWM_veg, by = c("country", "sitename", "trait"))

combo <- cwm %>% 
  dplyr::select(trait, zone) %>% 
  distinct() 

slopes <- purrr::map2_dfr(combo$zone, combo$trait, 
                          ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y,
                                   selectedtreatment =  "field_sp"))

# Switzerland ----
pollen_files <- 
  list.files("RDS_files/",full.names = TRUE) %>% 
  str_subset(pattern = "pollen") %>% 
  str_subset(pattern = "field") %>% 
  str_subset(pattern = "Switzerland")
bdm_meta <- read.csv("Data/bdm_metadata.csv") %>% 
  mutate(sitename = as.character(aID_STAO))
dfCWM_pol_field <-  purrr::map(pollen_files, ~readRDS(.x) %>% 
                                 as.data.frame())  %>% 
  bind_rows()  %>% 
  rownames_to_column(var = "standardized") %>% 
  # filter(str_detect(standardized, pattern = "zcwm"))  %>% 
  dplyr::select(country, sitename, pollen_mean = Mean, 
                pollen_sd = SD, trait, standardized) %>% 
  mutate(sitename = str_remove(sitename, pattern = "X"),
         standardized = ifelse(str_detect(standardized, "zCWM"),"yes","no")) %>% 
  left_join(bdm_meta, by = c("sitename"))

lm_bay <- function(cwm, selectedhabitat, selectedpft, selectedtrait, selecteddata){
  cwm <- cwm %>% 
    filter(habitat01 %in% selectedhabitat) %>% 
    filter(growthform == selectedpft) %>% 
    filter(trait == selectedtrait) %>% 
    left_join(sum_pft, by = c("sitename")) %>% 
    filter((growthform == "allpft") |
             (growthform == "trsh" & Woody >= 5) |
             (growthform == "herb" & `Non-Woody` >= 5)) 
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
             "_",selectedpft, "field.png"))
  par(mfrow = c(2,1))
  # Plot of residuals indicates a potential problem:
  qqnorm(residuals); qqline(residuals)
  # Looking at residuals vs fitted indicates increasing variance:
  plot(fitted, residuals); abline(h = 0)
  dev.off()
  
  # save convergence diagnostics
  plot(results, vars = c("zintercept", "zslope_elevation", "zsd", "df"), 
       file = paste0("Convergence_diagnostics/04_Covergence_lm_", selectedtrait,
                     "_", selecteddata,"_",paste(selectedhabitat, collapse = "_"),"_",selectedpft, "field.pdf" ))
  
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
                    data = selecteddata,
                    cwm_mean = pred_cwm_mean,
                    credible_lower = credible_lower,
                    credible_upper = credible_upper,
                    uncertainty_lower = uncertainty_lower,
                    uncertainty_upper = uncertainty_upper,
                    trait = selectedtrait)
  saveRDS(df_pred, paste0("RDS_files/05_lm_switzerland_", selectedtrait, "_", selecteddata,"_",paste(selectedhabitat, collapse = "_"),
                          "_", selectedpft,"field.rds"))

}
traits <- c("LA","SLA","PlantHeight")



purrr::map(traits, ~lm_bay(cwm = dfCWM_pol_field, 
                           selectedtrait = .x, 
                           selectedpft = "herb",
                           selectedhabitat = c("forests", "grasslands"), 
                           selecteddata = "Pollen"))
purrr::map(traits, ~lm_bay(cwm = dfCWM_pol_field, 
                           selectedtrait = .x, 
                           selectedpft = "trsh",
                           selectedhabitat = c("forests", "grasslands"), 
                           selecteddata = "Pollen"))
