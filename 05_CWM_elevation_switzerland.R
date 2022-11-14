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
##
## ---------------------------

set.seed(1)

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rjags)) install.packages("rjags")
if (!require(runjags)) install.packages("runjags")

## Load data
bdm_meta <- read_csv("Data/bdm_metadata.csv") %>% 
  mutate(sitename = as.character(aID_STAO))
# collate data
traits <- c("LA","SLA","PlantHeight")

pollen_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "03_CWM_estimates_") %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(pattern = "pollen") %>% 
  str_subset(pattern = "gapfilled", negate = TRUE) %>% 
  str_subset(pattern = "Switzerland") %>% 
  str_subset(pattern = "draw", negate = TRUE)

pollen_lab <- pollen_files %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_remove(., "03_CWM_estimates_") %>% 
  str_remove(., ".rds")

veg_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "03_CWM_estimates_") %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(pattern = "Switzerland") 

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
  mutate(correction = case_when(str_detect(label, pattern = "adjustedpercent_draw") ~ str_extract(label, pattern = "draw[0-9]*$"),
                                str_detect(label, pattern = "adjustedpercent_mean") ~ "correction",
                                str_detect(label, pattern = "percent") ~ "no correction",
                                TRUE ~ "no correction"),
         pollination = replace_na(pollination, "allpol"),
         growthform = replace_na(growthform, "allpft"),
         sitename = str_remove(sitename, "X")
  ) %>% 
  dplyr::select(country, sitename, cwm_mean = Mean, 
                cwm_sd = SD, trait, correction, pollination, growthform) %>% 
  # create treatment column
  mutate(treatment = case_when((pollination == "allpol" & growthform == "allpft") ~ "Correction factor",
                               (pollination != "allpol" & growthform == "allpft") ~ "Pollination mode",
                               (growthform != "allpft" & pollination == "allpol") ~ "Growth form")) %>% 
  left_join(bdm_meta, by = "sitename")

dfCWM_veg <- veg_files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(veg_lab, ~readRDS(.) %>% 
                as.data.frame() %>% 
                mutate(label = .y) %>% 
                mutate(across(one_of("sitename"),  as.character))) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(country = str_extract(label, pattern = "Scotland|Switzerland"),
         pollination = dplyr::recode(pollination, `no wind` = "not wind")) %>% 
  dplyr::select(country, sitename, trait, cwm_mean = Mean, cwm_sd = SD,
                growthform, pollination, taxres, zone) %>% 
  # create treatment column
  mutate(treatment = case_when((pollination == "allpol" & growthform == "allpft" &
                                  taxres == "stand.spec") ~ "Correction factor",
                               taxres %in% c("family","genus") ~ "Taxonomic resolution",
                               (pollination != "allpol" & growthform == "allpft") ~ "Pollination mode",
                               (growthform != "allpft" & pollination == "allpol") ~ "Growth form"),
         zone = str_remove(zone, pattern = "Scotland |Switzerland ")) %>% 
  left_join(bdm_meta, by = "sitename")

selectedtrait <- "PlantHeight"

saveRDS(dfCWM_pol, "RDS_files/04_CWM_estimates_pollen_Switzerland.rds") 
saveRDS(dfCWM_veg, "RDS_files/04_CWM_estimates_vegetation_Switzerland.rds") 

# Linear model ----
lm_bay <- function(cwm, selectedtrait, selecteddata){
  CWM <- cwm %>% 
    filter(trait == selectedtrait)
  Mean <- CWM %>% 
    pull(cwm_mean) 
  SD <- CWM %>%
    pull(cwm_sd) 
  Elev <- CWM %>% 
    pull(elevation)
  MeanElev <- mean(Elev)
  zElev <- Elev - MeanElev # mean center elevation 
  N <- length(Mean)
  
  model_string <- "
      model {
      for(i in 1:N) {
        Mean[i] ~ dt(regression_fitted[i], 1/SD[i]^2,df)
        # Mean[i] ~ dnorm(regression_fitted[i], 1/SD[i]^2)
        regression_fitted[i] <- intercept + slope_elevation * Elev[i]
        regression_residual[i] <- Mean[i] - regression_fitted[i]
      }
    
      # Priors
      intercept ~ dnorm(0, 10^-4)
      slope_elevation ~ dnorm(0, 10^-4)
      df ~ dexp(1/30) # df of t distribution, Values of d) larger than about 30.0 make the t distribution roughly normal
      
      #data# N, Mean, SD, Elev
      #monitor# intercept, slope_elevation
      #residual# regression_residual
      #fitted# regression_fitted
      #inits# intercept, slope_elevation
    }"
  
  intercept <- list(chain1 = min(Mean), chain2 = max(Mean))
  slope_elevation <- list(chain1 = c(-10), chain2 = c(10))
  
  results <- run.jags(model_string, n.chains = 2, 
                      modules = "glm on")
  plot(results)
  
  results_summary <- as.data.frame(summary(results)) 
  list(model = results,
         summary = results_summary)
  
  # plot model results
  elev_new <- seq(500, 2500, length.out = 100) # elevations for model predictions
  lm_mcmc <- as.mcmc(results)
  pred_cwm_mean <- mean(lm_mcmc[,"intercept"]) + elev_new * mean(lm_mcmc[,"slope_elevation"])
  pred_mean_dist <- matrix(NA, nrow = nrow(lm_mcmc), ncol = length(elev_new))
  for (i in 1:nrow(pred_mean_dist)) {
    pred_mean_dist[i,] <- lm_mcmc[i,"intercept"] + elev_new * lm_mcmc[i,"slope_elevation"]
  }
  credible_lower <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.05)
  credible_upper <- apply(pred_mean_dist, MARGIN = 2, quantile, prob = 0.95)
  df_pred <- tibble(elevation = elev_new,
               cwm_mean = pred_cwm_mean,
               lower = credible_lower,
               upper = credible_upper)
  p <- ggplot(df_pred, aes(x = elevation, y = cwm_mean)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = "grey") +
    geom_line(color = "darkgreen") +
    geom_point(data = CWM, aes(x = elevation, y = cwm_mean)) +
    scale_y_continuous("") + 
    scale_x_continuous("Elevation (m)") +
    ggtitle( paste(labs_trait[selectedtrait], selecteddata)) +
    theme_bw()
  ggsave(paste0("Figures/CWM_elevation_", selectedtrait,"_", selecteddata, ".png"),p)
}


labs_trait <- c("Height (log)(cm)",
                "Leaf area (log)(mm2)",
                "SLA (mm2/mg)")
names(labs_trait) <- c("PlantHeight", "LA", "SLA")

cwm_veg <- dfCWM_veg %>% 
  filter(treatment == "Correction factor") 
lm_veg <- purrr::map(names(labs_trait), ~lm_bay(cwm = cwm_veg, selectedtrait = .x, selecteddata = "Vegetation"))

cwm_pol <- dfCWM_pol %>% 
  filter(treatment == "Correction factor") %>% 
  drop_na(elevation) %>% 
  filter(correction == "no correction")
lm_pol <- purrr::map(names(labs_trait), ~lm_bay(cwm = cwm_pol, selectedtrait = .x, selecteddata = "Pollen"))


# Check residuals ----

# Extract residuals and fitted estimates:
residuals <- resid(lm_veg[[1]]$model)
fitted <- fitted(results)

# Plot of residuals indicates a potential problem:
qqnorm(residuals); qqline(residuals)

# Looking at residuals vs fitted indicates increasing variance:
plot(fitted, residuals); abline(h = 0)

cwm_veg %>%
  filter(trait == "PlantHeight") %>% 
  add_epred_draws(lm_veg[[1]]) %>%
  ggplot(aes(x = elevation, y = cwm_mean)) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = cwm_veg) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")
plot(lm_veg[[1]]$model)

dfCWM_veg %>% 
  filter(treatment == "Correction factor")  %>% 
ggplot(.,aes(x = elevation, y = cwm_mean)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  facet_wrap(.~trait, scales = "free") +
  theme_bw()
dfCWM_pol %>% 
  filter(treatment == "Correction factor") %>% 
  ggplot(.,aes(x = elevation, y = cwm_mean, color = correction)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  facet_wrap(.~trait, scales = "free") +
  theme_bw()

dfCWM %>% 
  filter(treatment == "Growth form") %>% 
  ggplot(.,aes(x = elevation, y = veg_mean, color = growthform)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  facet_wrap(.~trait, scales = "free") +
  theme_bw()
dfCWM %>% 
  filter(treatment == "Growth form") %>% 
  ggplot(.,aes(x = elevation, y = pollen_mean, color = growthform)) +
  geom_point()  +
  geom_smooth(method = "lm") +
  facet_wrap(.~trait, scales = "free") +
  theme_bw()
