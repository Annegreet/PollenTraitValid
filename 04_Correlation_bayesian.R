## ---------------------------
##
## Script name: 04_Correlation_bayesian
##
## Purpose of script: Calculate correlation between vegetation and pollen CWM for different subsets
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
## References:
## http://doingbayesiandataanalysis.blogspot.com/2015/12/prior-on-df-normality-parameter-in-t.html
## https://solomonkurz.netlify.app/post/2019-02-16-bayesian-correlations-let-s-talk-options/
## ---------------------------

set.seed(1)
options(scipen = 9999)

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rjags)) install.packages("rjags")
if (!require(runjags)) install.packages("runjags")
if (!require(furrr)) install.packages("furrr")

## Load data
# load PFT summary table
# to filter sites with to little vegetation cover of particular type
sum_pft <- readRDS("RDS_files/01_Percentage_cover_pft.rds") %>% 
  replace_na(list(`Non-Woody` = 0,
                  Woody = 0,
                  `not wind` = 0,
                  wind = 0)) 

# collate data
traits <- c("LA","SLA","PlantHeight")

pollen_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "03_zCWM_estimates_") %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(pattern = "pollen") %>% 
  str_subset(pattern = "notinpollen", negate = TRUE) 

pollen_lab <- pollen_files %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_remove(., "03_zCWM_estimates_") %>% 
  str_remove(., ".rds")

veg_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "03_zCWM_estimates_") %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(pattern = "pollen", negate = TRUE) 

veg_lab <- veg_files %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_remove(., "03_zCWM_estimates_") %>% 
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
                                str_detect(label, pattern = "adjusted_helinger") ~ "helinger correction",
                                TRUE ~ "no correction"),
         pollination = replace_na(pollination, "allpol"),
         growthform = replace_na(growthform, "allpft"),
         sitename = str_remove(sitename, "X")
  ) %>%
  # don't include cwm values with removed species (supplementary analysis)
  # filter(is.na(removed_spec)) %>% 
  dplyr::select(country, sitename, pollen_mean = Mean, 
                pollen_sd = SD, trait, correction, pollination, growthform
  ) 


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
         notinpollen = if_else(str_detect(label, "notinpollen"), "notinpollen", "inpollen")) %>%
  # only Scottish data
  filter(country == "Scotland") %>% 
  # remove sites with to little cover of pft (5% or less) 
  left_join(sum_pft, by = c("sitename", "zone")) %>% 
  filter((growthform == "allpft" & pollination == "allpol") |
           (growthform == "trsh" & Woody >= 5) |
           (growthform == "herb" & `Non-Woody` >= 5) |
           (pollination == "wind" & wind >= 5) |
           (pollination == "not wind" & `not wind` >= 5)) %>% 
  dplyr::select(country, sitename, trait, veg_mean = Mean, veg_sd = SD,
                growthform, pollination, taxres, zone, notinpollen) 

saveRDS(dfCWM_pol,"RDS_files/04_zCWM_estimates_pollen_Scotland.rds")
saveRDS(dfCWM_veg,"RDS_files/04_zCWM_estimates_vegetation_Scotland.rds")


selectedtrait <- "PlantHeight"
selectedzone <- "zoneB"

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
                              "_", selectedzone,"_",selectedtreatment, ".pdf" ))
  # gelman plots
  pdf(paste0("Convergence_diagnostics/04_Gelman_regression_", selectedtrait,
             "_", selectedzone,"_",selectedtreatment, ".pdf"))
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  results_summary <- as.data.frame(summary(results)) %>% 
    rownames_to_column("parameter") %>% 
    mutate(trait = selectedtrait,
           zone = selectedzone,
           treatment = selectedtreatment,
           cor_nonbay = cor(pollen_cwm, veg_cwm))
  
  saveRDS(results_summary, paste0("RDS_files/04_slope_regression_", selectedtrait,
                                  "_", selectedzone, "_", selectedtreatment, ".rds"))
  return(results_summary)
  
}

trait_lab <-  c("Leaf area", "Plant height", 
                "Specific leaf area")
names(trait_lab) <- c("LA", "PlantHeight","SLA")
zone_lab <- c("Inner ring", "Middle ring", "Outer ring")
names(zone_lab) <- c("zoneA","zoneB","zoneC")

## Plotting function
corrplot <- function(slopes, nobs, labels){
  slopes %>% 
    # gfilter for the slope estimates
    filter(parameter == "slope") %>% 
    # make sure the order of zones is right
    mutate(zone = fct_relevel(zone, "zoneC", "zoneB","zoneA")) %>% 
    ggplot(aes(x = zone, y = Mean, fill = treatment)) + 
    # geoms
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_errorbar(aes(ymin = Lower95, ymax = Upper95), linewidth = 0.5, width = 0.1,
                  position = position_dodge(0.5)) +
    geom_point(shape = 21, size = 2, position = position_dodge(0.5)) +
    # labels and scales
    scale_x_discrete("", labels = c("zoneA" = "Inner ring (0-10m)", 
                                    "zoneB" = "Middle ring (0-100m)", 
                                    "zoneC"  = "Outer ring (0-1000m)")) +
    scale_y_continuous("β", limits = c(-1,1)) +
    scale_fill_manual("", values = c("darkorange", "purple", "cyan4"),
                      labels = labels) +
    coord_flip() +
    # Faceting
    facet_grid(~trait, labeller = labeller(trait = trait_lab)) +
    # Number of observations
    geom_text(data = nobs, aes(x = x, y = Mean,label = label),
              nudge_x = 0.25, size = 2, inherit.aes = FALSE) +
    # Theme
    theme_bw() +
    theme(text = element_text(size = 12),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom",
          legend.margin = margin(-15,0,0,0)
          )
}

if(0){
## No correction + correction----
# No correction
cwm <- dfCWM_pol %>%
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "no correction") %>% 
  left_join(dfCWM_veg, by = c("country", "sitename", "trait", "growthform", "pollination")) %>% 
  # plant height in zone B has bimodal distribution, remove from dataset
  filter(!(trait == "PlantHeight" & zone == "zoneB")) %>% 
  filter(taxres == "stand.spec") %>% 
  filter(notinpollen == "inpollen")

combo <- cwm %>% 
  dplyr::select(trait, zone) %>% 
  distinct() 

slopes <- purrr::map2_dfr(combo$zone, combo$trait, 
                        ~cor_bay(cwm, selectedzone = .x, selectedtrait = .y,
                                 selectedtreatment =  "no correction"))

# Pollen correction 
cwm2 <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "correction") %>% 
  left_join(dfCWM_veg, by = c("country","sitename", "trait", "growthform", "pollination")) %>% 
  # plant height in zone B has bimodal distribution, remove from dataset
  filter(!(trait == "PlantHeight" & zone == "zoneB")) %>% 
  filter(taxres == "stand.spec") %>% 
  filter(notinpollen == "inpollen")

combo <- cwm %>% 
  dplyr::select(trait, zone) %>% 
  distinct() 

slopes2 <- purrr::map2_dfr(combo$zone, combo$trait, 
                         ~cor_bay(cwm2, selectedzone = .x, selectedtrait = .y,
                                  selectedtreatment =  "correction")) 

# Pollen correction 
cwm3 <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "helinger correction") %>% 
  left_join(dfCWM_veg, by = c("country","sitename", "trait", "growthform", "pollination")) %>% 
  # plant height in zone B has bimodal distribution, remove from dataset
  filter(!(trait == "PlantHeight" & zone == "zoneB")) %>% 
  filter(taxres == "stand.spec") %>% 
  filter(notinpollen == "inpollen")

combo <- cwm3 %>% 
  dplyr::select(trait, zone) %>% 
  distinct() 

slopes3 <- purrr::map2_dfr(combo$zone, combo$trait, 
                           ~cor_bay(cwm3, selectedzone = .x, selectedtrait = .y,
                                    selectedtreatment =  "helinger correction")) %>% 
  bind_rows(slopes) %>% 
  bind_rows(slopes2)

# df with number of observations 
nobs <- cwm2 %>% 
  bind_rows(cwm) %>% 
  bind_rows(cwm3) %>% 
  # plant height in zone B has bimodal distribution, remove from dataset
  filter(!(trait == "PlantHeight" & zone == "zoneB")) %>% 
  group_by(trait, zone, correction) %>% 
  summarise(label = n(), .groups = "keep") %>% 
  left_join(slopes3[slopes3$parameter == "slope", ], by = c("trait",  "zone", "correction" = "treatment")) %>% 
  dplyr::select(label,  zone, trait, Mean, treatment = "correction") %>% 
  ungroup() %>% 
  arrange(trait) %>% 
  mutate(x = c(2.66,2.84,3,
               1.66,1.84,2,
               0.66,0.84,1,
               2.66,2.84,3,
               0.66,0.84,1,
               2.66,2.84,3,
               1.66,1.84,2,
               0.66,0.84,1))

# plot
p_slope <- corrplot(slopes = slopes3, nobs, labels = c("correction" = "Pollen productivity estimates",
                                                       "no correction" = "No correction",
                                                       "helinger correction" = "Helinger transformation")) 

ggsave("Figures/Pollen_correction_factors.png", p_slope, width = 7, height = 4.5)

## Pollination mode ----
# Wind
cwm_wind <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "wind",
         correction == "no correction") %>% 
  left_join(dfCWM_veg, by = c("country", "sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 

combo <- cwm_wind %>% 
  dplyr::select(trait, zone) %>% 
  distinct() 

slopes_wind <- purrr::map2_dfr(combo$zone, combo$trait, 
                                    ~cor_bay(cwm_wind, selectedzone = .x, selectedtrait = .y,
                                             selectedtreatment = "wind"))
# No wind
cwm_nowind <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "not wind",
         correction == "no correction") %>% 
  left_join(dfCWM_veg, by = c("country", "sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec")

slopes_pol <- purrr::map2_dfr(combo$zone, combo$trait, 
                                   ~cor_bay(cwm_nowind, selectedzone = .x, selectedtrait = .y,
                                            selectedtreatment =  "not wind")) %>% 
  bind_rows(slopes_wind) 

# df with number of observations 
nobs_pol <- cwm_nowind %>% 
  bind_rows(cwm_wind) %>% 
  group_by(trait, zone, pollination) %>% 
  summarise(label = n(), .groups = "keep") %>% 
  left_join(slopes_pol[slopes_pol$parameter == "slope", ],
            by = c("trait",  "zone", "pollination" = "treatment")) %>% 
  dplyr::select(label,  zone, trait, Mean, treatment = "pollination") %>% 
  ungroup() %>% 
  arrange(trait) %>% 
  mutate(x = rep(c(2.75,3,1.75,2,0.75,1), 3))

p_slope <- corrplot(slopes_pol, nobs_pol, c( "wind" = "Wind pollinated",
                                         "not wind" = "Not wind pollinated"))
ggsave("Figures/Pollination_mode.png", p_slope, width = 7, height = 3.5)

## Plant functional type ----
# Trees and Shrubs
cwm_trsh <- dfCWM_pol %>% 
  filter(growthform == "trsh",
         pollination == "allpol",
         correction == "no correction") %>% 
  left_join(dfCWM_veg, by = c("country", "sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 

combo_trsh <- cwm_trsh %>% 
  dplyr::select(trait, zone) %>% 
  distinct()  # no zone A in the tree data

slopes_trsh <- map2_dfr(combo_trsh$zone, combo_trsh$trait, 
                      ~cor_bay(cwm_trsh, selectedzone = .x, selectedtrait = .y,
                               selectedtreatment = "trsh"))
# Herbs
cwm_herb <- dfCWM_pol %>% 
  filter(growthform == "herb",
         pollination == "allpol",
         correction == "no correction") %>% 
  left_join(dfCWM_veg, by = c("country","sitename", "trait", "growthform", "pollination")) %>% 
  filter(taxres == "stand.spec") 

combo <- cwm_herb %>% 
  dplyr::select(trait, zone) %>% 
  distinct()
slopes_pft <- map2_dfr(combo$zone, combo$trait, 
                     ~cor_bay(cwm_herb, selectedzone = .x, selectedtrait = .y,
                              selectedtreatment =  "herb")) %>% 
  bind_rows(slopes_trsh) 

# df with number of observations 
nobs_pft <- cwm_herb %>% 
  bind_rows(cwm_trsh) %>% 
  group_by(trait, zone, growthform) %>% 
  summarise(label = n(), .groups = "keep") %>% 
  left_join(slopes_pft[slopes_pft$parameter == "slope", ],
            by = c("trait",  "zone", "growthform" = "treatment")) %>% 
  dplyr::select(label,  zone, trait, Mean, treatment = "growthform") %>% 
  ungroup() %>% 
  arrange(trait) %>% 
  mutate(x = rep(c(2.86,1.75,2,0.75,1), 3))

p_slope <- corrplot(slopes_pft, nobs_pft, c("trsh" = "Woody", "herb" = "Non-Woody"))
ggsave("Figures/PFT.png", p_slope, width = 7, height = 3.5)

## Taxonomic resolution ----
# Species
cwm_sp <- cwm
slopes_sp <- slopes %>% 
  mutate(treatment = "stand.spec")

# Genus
cwm_gen <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "no correction") %>% 
  left_join(dfCWM_veg, by = c("country", "sitename", "trait", "growthform", "pollination")) %>% 
  # plant height in zone B has bimodal distribution, remove from dataset
  filter(!(trait == "PlantHeight" & zone == "zoneB")) %>% 
  filter(taxres == "genus") 

combo <- cwm_gen %>% 
  dplyr::select(trait, zone) %>% 
  distinct()

slopes_gen <- purrr::map2_dfr(combo$zone, combo$trait, 
                                   ~cor_bay(cwm_gen, selectedzone = .x, selectedtrait = .y,
                                            selectedtreatment =  "genus"))
# family
cwm_fam <- dfCWM_pol %>% 
  filter(growthform == "allpft",
         pollination == "allpol",
         correction == "no correction") %>% 
  left_join(dfCWM_veg, by = c("country","sitename", "trait", "growthform", "pollination")) %>% 
  # plant height in zone B has bimodal distribution, remove from dataset
  filter(!(trait == "PlantHeight" & zone == "zoneB")) %>% 
  filter(taxres == "family") 

slopes_taxres <- purrr::map2_dfr(combo$zone, combo$trait, 
                                      ~cor_bay(cwm_fam, selectedzone = .x, selectedtrait = .y,
                                               selectedtreatment =  "family")) %>% 
  bind_rows(slopes_gen) %>%
  bind_rows(slopes_sp) 

# df with number of observations 
nobs_taxres <- cwm_fam %>%
  bind_rows(cwm_gen) %>%
  bind_rows(cwm_sp) %>%
  group_by(trait, zone, taxres) %>%
  summarise(label = n(), .groups = "keep") %>%
  # plant height in zone B has bimodal distribution, remove from dataset
  filter(!(trait == "PlantHeight" & zone == "zoneB")) %>% 
  left_join(slopes_taxres[slopes_taxres$parameter == "slope", ], by = c("trait",  "zone", "taxres" = "treatment")) %>%
  dplyr::select(label,  zone, trait, Mean, treatment = "taxres") %>%
  ungroup() %>%
  arrange(trait) %>% 
  mutate(x = c(2.66,2.84,3,
               1.66,1.84,2,
               0.66,0.84,1,
               2.66,2.84,3,
               0.66,0.84,1,
               2.66,2.84,3,
               1.66,1.84,2,
               0.66,0.84,1))

p_slope <- corrplot(slopes_taxres, nobs_taxres, c("family" = "Family",
                                              "genus" = "Genus",
                                              "stand.spec" = "Species"))

ggsave("Figures/Taxonomic_resolution.png", p_slope, width = 7, height = 4.5)
}


## Pollen correction + uncertainty ----
cor_bay_draw <- function(cwm, selectedzone, selectedtrait, 
                         selectedtreatment, selecteddraw){
  CWM <- cwm %>% 
    filter(zone == selectedzone,
           trait == selectedtrait,
           correction == selecteddraw)
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
  
  results_summary <- as.data.frame(summary(results)) %>% 
    rownames_to_column("parameter") %>% 
    mutate(trait = selectedtrait,
           draw = selecteddraw,
           zone = selectedzone,
           treatment = selectedtreatment)
  return(results_summary)
  saveRDS(results, paste0("RDS_files/04_regression_", selectedtrait,
                          "_", selectedzone,"_", selectedtreatment,"_",
                          selecteddraw, ".rds"))
}

cwm_draw <- dfCWM_pol %>%
  filter(growthform == "allpft",
         pollination == "allpol",
         str_detect(correction, pattern = "draw")) %>%
  left_join(dfCWM_veg, by = c("country","sitename", "trait", "growthform", "pollination")) %>% 
  # plant height in zone B has bimodal distribution, remove from dataset
  filter(!(trait == "PlantHeight" & zone == "zoneB")) %>% 
  filter(taxres == "stand.spec")

combo <- expand_grid(zone = unique(cwm_draw$zone), correction = unique(cwm_draw$correction))

slopes_LA <- purrr::map2_dfr(combo$zone, combo$correction, 
                                  ~cor_bay_draw(cwm_draw, selectedzone = .x, selecteddraw = .y,
                                                selectedtrait = "LA",
                                                selectedtreatment = "draw"))
slopes_SLA <- purrr:::map2_dfr(combo$zone, combo$correction, 
                                   ~cor_bay_draw(cwm_draw, selectedzone = .x, selecteddraw = .y,
                                                 selectedtrait = "SLA",
                                                 selectedtreatment = "draw"))
slopes_PH <- purrr::map2_dfr(combo$zone, combo$correction, 
                                  ~cor_bay_draw(cwm_draw, selectedzone = .x, selecteddraw = .y,
                                                selectedtrait = "PlantHeight",
                                                selectedtreatment = "draw"))
slopes_draw <- slopes_LA %>% 
  bind_rows(slopes_SLA) %>%
  bind_rows(slopes_PH) %>% 
  filter(parameter == "slope")

(p_slope <- 
    slopes_draw %>% 
    # make sure the order of zones is right
    mutate(zone = fct_relevel(zone, "zoneC", "zoneB","zoneA")) %>% 
    # plant height in zone B has bimodal distribution, remove from dataset
    filter(!(trait == "PlantHeight" & zone == "zoneB")) %>% 
    ggplot(aes(x = zone, y = Mean)) + 
    # geoms
    geom_hline(yintercept = 0, linetype = "dashed") +
    # add points
    geom_point(shape = 21, size = 2, alpha = 0.5, fill = "darkorange") +
    # labels and scales
    scale_x_discrete("", labels = c("zoneA" = "Inner ring (0-10m)", 
                                    "zoneB" = "Middle ring (0-100m)", 
                                    "zoneC"  = "Outer ring (0-1000m)")) +
    scale_y_continuous("β", limits = c(-1,1)) +
    coord_flip() +
    # Faceting
    facet_grid(~trait, labeller = labeller(trait = trait_lab)) +
    # Theme
    theme_bw() +
    theme(text = element_text(size = 12),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom",
          legend.margin = margin(-15,0,0,0))
) 
ggsave("Figures/Correction_draw.png", p_slope, width = 7, height = 3.5)
