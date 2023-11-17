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
# Trait data  
dfTRAIT <- readRDS("RDS_files/02_Gapfilled_traits_pollen.rds")

# Pollen data 
dfPOL_Scot <- readRDS("RDS_files/01_pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_pollen_data_Swiss.rds")
dfPOL <- bind_rows("Scotland" = dfPOL_Scot, "Switzerland" =
                     dfPOL_Swiss, .id = "country") %>% 
  ungroup()


# Params 
# selectedcountry <- "Scotland"
# selectedtrait <- "PlantHeight"
# selectedabun <- str_subset(colnames(dfPOL), "percent")[1]

normaltraits <- c("PlantHeight","SLA","LA")

# set number of cores to use
plan(multisession)

# No subset ----
if (0) {

furrr::future_map(normaltraits,
                   ~cwm_pol(selectedtrait = .x,
                            selectedcountry = "Scotland",
                            selectedabun = "percent",
                            taxtrait = TRUE
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
                           selectedcountry = "Scotland",
                           selectedabun = "adjusted_helinger"),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(normaltraits,
                    ~cwm_pol(selectedtrait = .x,
                                selectedcountry = "Switzerland",
                                selectedabun = "adjustedpercent_mean"))
  
}


# Growth form ----
if (0) {

trsh <- c("tree", "shrub")
herb <- c("grass", "herb")

purrr::map(normaltraits,
                  ~cwm_pol_pft(selectedtrait = .x,
                              selectedcountry = "Scotland",
                              selectedpft = trsh,
                              selectedabun = "percent"))
purrr::map(normaltraits, 
                  ~cwm_pol_pft(selectedtrait = .x, 
                               selectedcountry = "Scotland",
                               selectedpft = trsh,
                               selectedabun = "adjustedpercent_mean"))
purrr::map(normaltraits,
                  ~cwm_pol_pft(selectedtrait = .x,
                               selectedcountry = "Scotland",
                               selectedpft = herb,
                               selectedabun = "percent"))
purrr::map(normaltraits, 
                  ~cwm_pol_pft(selectedtrait = .x, 
                               selectedcountry = "Scotland",
                               selectedpft = herb,
                               selectedabun = "adjustedpercent_mean"))
purrr::map(normaltraits,
           ~cwm_pol_pft(selectedtrait = .x,
                        selectedcountry = "Switzerland",
                        selectedpft = herb,
                        selectedabun = "percent"))
purrr::map(normaltraits,
           ~cwm_pol_pft(selectedtrait = .x,
                        selectedcountry = "Switzerland",
                        selectedpft = trsh,
                        selectedabun = "percent"))


}
# Pollination mode ----
if (0) {

purrr::map(normaltraits,
                  ~cwm_pol_pol(selectedtrait = .x,
                              selectedcountry = "Scotland",
                              selectedpolmode = "wind",
                              selectedabun = "percent"))

purrr::map(normaltraits,
                    ~cwm_pol_pol(selectedtrait = .x,
                                 selectedcountry = "Scotland",
                                 selectedpolmode = "wind",
                                 selectedabun = "adjustedpercent_mean"))
purrr::map(normaltraits,
                  ~cwm_pol_pol(selectedtrait = .x,
                              selectedcountry = "Scotland",
                              selectedpolmode = "not wind",
                              selectedabun = "percent"),
                  .options = furrr_options(seed = TRUE))
purrr::map(normaltraits,
                  ~cwm_pol_pol(selectedtrait = .x,
                               selectedcountry = "Scotland",
                               selectedpolmode = "not wind",
                               selectedabun = "adjustedpercent_mean"))

purrr::map(normaltraits,
           ~cwm_pol_pol(selectedtrait = .x,
                        selectedcountry = "Switzerland",
                        selectedpolmode = "wind",
                        selectedabun = "percent"))
purrr::map(normaltraits,
           ~cwm_pol_pol(selectedtrait = .x,
                        selectedcountry = "Switzerland",
                        selectedpolmode = "not wind",
                        selectedabun = "percent"),
           .options = furrr_options(seed = TRUE))

}

# Pollen correction factors ----
if (1) {
  # create data frame with combinations that have runned already
  file_names <- list.files("RDS_files") %>%
    str_subset("03_CWM_estimates") %>%
    str_subset("draw") %>%
    str_subset("gapfilled", negate = T) %>%
    str_remove("03_CWM_estimates_pollen_Scotland_") %>%
    str_remove(".rds") %>%
    data.frame(label = .) %>%
    separate(., label, into = c("trait", "pollen_abun"), sep = "_", extra = "merge") %>%
    mutate(run = "yes")
  combo <- expand_grid(trait = normaltraits, pollen_abun = 
                         str_subset(colnames(dfPOL), "draw")) %>% 
    left_join(file_names, by = c("trait", "pollen_abun")) %>% 
    filter(is.na(run))


  future_map2(combo$trait, combo$pollen_abun,
                    ~cwm_pol(selectedtrait = .x,
                             selectedcountry = "Scotland",
                             selectedabun = .y),
              .options = furrr_options(seed = TRUE))
}

