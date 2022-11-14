## ---------------------------
##
## Script name: 03_CWM_estimates_pollen.R
##
## Purpose of script: Calculate CWM values of pollen data
##
## Author: Annegreet Veeken
##
## Date Created: 2022-07-18
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

set.seed(1)

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rjags))  install.packages("rjags")
if (!require(runjags)) install.packages("runjags")
if (!require(mcmc)) install.packages("mcmc")
if (!require(coda)) install.packages("coda")
if (!require(furrr)) install.packages("furrr")

source("cwm_functions.R")

## Load and prepare data -----
# Trait data ---- 
dfTRAIT <- readRDS("RDS_files/02_Gapfilled_traits_pollen.rds")

# Pollen data ----
dfPOL_Scot <- readRDS("RDS_files/01_pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_pollen_data_Swiss.rds")
dfPOL <- bind_rows("Scotland" = dfPOL_Scot, "Switzerland" =
                     dfPOL_Swiss, .id = "country") %>% 
  ungroup()

selectedcountry <- "Scotland"
selectedtrait <- "PlantHeight"
selectedabun <- str_subset(colnames(dfPOL), "percent")[1]

normaltraits <- c("PlantHeight","SLA","LA")

# set number of cores to use
if (Sys.info()[4] == "LUIN67311") {
  plan(multisession, workers = 1)
} else {
  plan(multisession, workers = 6)
}

# No subset ----
if (1) {

  furrr::future_map(normaltraits,
                     ~cwm_pol(selectedtrait = .x,
                              selectedcountry = "Scotland",
                              selectedabun = "percent"
                              # taxtrait = TRUE
                              ),
                    .options = furrr_options(seed = TRUE))
  furrr::future_map(normaltraits,
                     ~cwm_pol(selectedtrait = .x,
                              selectedcountry = "Switzerland",
                              selectedabun = "percent",
                              taxtrait = TRUE
                              ),
                    .options = furrr_options(seed = TRUE))
  furrr::future_map(normaltraits,
                    ~cwm_pol(selectedtrait = .x,
                                selectedcountry = "Scotland",
                                selectedabun = "adjustedpercent_mean"),
                    .options = furrr_options(seed = TRUE))
  furrr::future_map(normaltraits,
                    ~cwm_pol(selectedtrait = .x,
                                selectedcountry = "Switzerland",
                                selectedabun = "adjustedpercent_mean"),
                    .options = furrr_options(seed = TRUE))
  
}


# Growth form ----
if (1) {

trsh <- c("tree", "shrub")
herb <- c("grass", "herb")

furrr::future_map(normaltraits, 
                  ~cwm_pol_pft(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpft = trsh,
                              selectedabun = "percent"),
                  .options = furrr_options(seed = TRUE))

furrr::future_map(normaltraits,
                  ~cwm_pol_pft(selectedtrait = .x,
                              selectedcountry = "Switzerland",
                              selectedpft = trsh,
                              selectedabun = "percent"),
                  .options = furrr_options(seed = TRUE))

furrr::future_map(normaltraits, 
                  ~cwm_pol_pft(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpft = herb,
                              selectedabun = "percent"),
                  .options = furrr_options(seed = TRUE))

furrr::future_map(normaltraits,
                  ~cwm_pol_pft(selectedtrait = .x,
                              selectedcountry = "Switzerland",
                              selectedpft = herb,
                              selectedabun = "percent"),
                  .options = furrr_options(seed = TRUE))

}
# Pollination mode ----
if (1) {

furrr::future_map(normaltraits,
                  ~cwm_pol_pol(selectedtrait = .x,
                              selectedcountry = "Scotland",
                              selectedpolmode = "wind",
                              selectedabun = "percent"),
                  .options = furrr_options(seed = TRUE))

furrr::future_map(normaltraits,
                  ~cwm_pol_pol(selectedtrait = .x,
                              selectedcountry = "Switzerland",
                              selectedpolmode = "wind",
                              selectedabun = "percent"),
                  .options = furrr_options(seed = TRUE))

furrr::future_map(normaltraits,
                  ~cwm_pol_pol(selectedtrait = .x,
                              selectedcountry = "Scotland",
                              selectedpolmode = "not wind",
                              selectedabun = "percent"),
                  .options = furrr_options(seed = TRUE))

furrr::future_map(normaltraits,
                  ~cwm_pol_pol(selectedtrait = .x,
                              selectedcountry = "Switzerland",
                              selectedpolmode = "not wind",
                              selectedabun = "percent"),
                  .options = furrr_options(seed = TRUE))
}

# Pollen correction factors ----
if (1) {
  
  combo <- expand_grid(trait = normaltraits, pollen_abun = str_subset(colnames(dfPOL), "percent"))
  

  furrr::future_map2(combo$trait, combo$pollen_abun,
                    ~cwm_pol(selectedtrait = .x,
                             selectedcountry = "Scotland",
                             selectedabun = .y))

  furrr::future_map2(combo$trait, combo$pollen_abun,
                     ~cwm_pol(selectedtrait = .x,
                                 selectedcountry = "Switzerland",
                                 selectedabun = .y))
}