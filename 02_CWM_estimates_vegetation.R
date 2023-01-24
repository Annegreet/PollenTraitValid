## ---------------------------
##
## Script name:02_CWM_estimates_vegetation.R
##
## Purpose of script: Calculating CWM values of the vegetation
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

## Load functions
source("cwm_functions.R")

## Load and prepare data -----
# Trait data ----
# Scotland ----
plantheight <- readRDS("RDS_files/01_PlantHeight.rds") %>% 
  ungroup() 
sla <- readRDS("RDS_files/01_SLA.rds") %>% 
  ungroup()
la <- readRDS("RDS_files/01_LA.rds") %>% 
  ungroup()


lTrait <- list(PlantHeight = plantheight,
               SLA = sla,
               LA = la)

dfTRY_raw <- readRDS("RDS_files/02_Gapfilled_traits.rds")

# Pollination and PFT ----
dfPFT <- readRDS("RDS_files/01_Polmode_pft_vegetation.rds")

# Species to pollen table
polspec <- readRDS("RDS_files/02_PollenType_species.rds") %>% 
  dplyr::select(pollentaxon, stand.spec) %>% distinct() %>% as_tibble() %>% 
  drop_na()

# Species abundance data ----
dfABUN_a <- readRDS("RDS_files/01_Species_abundance_a.rds") 
dfABUN_bc <- readRDS("RDS_files/01_Species_abundance_bc.rds")

zoneA <- dfABUN_a %>%
  arrange(sitename, stand.spec) %>% 
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  rename(family = fam) %>% 
  ungroup %>% 
  # join with pollen to species table
  left_join(polspec, by = "stand.spec") %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) 
zoneB <- dfABUN_bc %>% 
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  ungroup() %>% 
  dplyr::select(sitename, stand.spec, genus, family = fam, abun = spec_abun_b, growthform, polmode) %>% 
  # join with pollen to species table
  left_join(polspec, by = "stand.spec") %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) 
zoneC <- dfABUN_bc %>% 
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  ungroup %>% 
  dplyr::select(sitename, stand.spec, genus, family = fam, abun = spec_abun_c, growthform, polmode) %>% 
  # join with pollen to species table
  left_join(polspec, by = "stand.spec") %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) 
lABUN <- list("zoneA" = zoneA, "zoneB" = zoneB, "zoneC" = zoneC)

trsh <- c("tree", "shrub")
herb <- c("herb", "grass", "fern")
allpft  <- c("tree","herb","grass","fern","shrub")
wind <- "wind"
nowind <- "not wind"
allpol <- c("wind", "not wind")

# Create joining column for trait data (to allow match for taxa 
# that are not identified to species level)
veg_sp <- lABUN %>% 
  purrr::map(., ~pull(.,stand.spec)) %>% 
  unlist %>% unique
lTrait <-
  lTrait %>% 
  purrr::map(., ~mutate(., stand.spec = 
                          case_when(stand.spec %in% veg_sp ~ stand.spec,
                                    genus %in% veg_sp ~ genus,
                                    family %in% veg_sp ~ family)
           )
  ) 

dfTRY <- 
  dfTRY_raw %>% 
  mutate(., stand.spec = case_when(stand.spec %in% veg_sp ~ stand.spec,
                                   genus %in% veg_sp ~ genus,
                                   family %in% veg_sp ~ family)
                        ) %>% 
  drop_na()
             
selectedtrait <- "LA"
selectedabun <- "zoneB"
selectedtaxres <- "stand.spec"
selectedpolmode <- allpol
selectedpft <- allpft
notinpollen <- FALSE

scottraits <- c("PlantHeight", "LA", "SLA")

plan(multisession, workers = 3)
# Scotland----
if (1) { 
# Zone A ----
# all
 # furrr::future_map(scottraits,
 #                   ~cwm_veg_scot(selectedabun = "zoneA",
 #                             selectedtrait = .,
 #                             selectedtaxres = "stand.spec",
 #                             selectedpolmode = allpol,
 #                             selectedpft = allpft),
 #                   .options = furrr_options(seed = TRUE))
 # # only taxa in pollen data
 # furrr::future_map(scottraits,
 #                   ~cwm_veg_scot(selectedabun = "zoneA",
 #                             selectedtrait = .,
 #                             selectedtaxres = "stand.spec",
 #                             selectedpolmode = allpol,
 #                             selectedpft = allpft,
 #                             notinpollen = TRUE),
 #                   .options = furrr_options(seed = TRUE))
 # # pollination
 # furrr::future_map(scottraits,
 #                   ~cwm_veg_scot(selectedabun = "zoneA",
 #                             selectedtrait = .,
 #                             selectedtaxres = "stand.spec",
 #                             selectedpolmode = nowind,
 #                             selectedpft = allpft),
 #                   .options = furrr_options(seed = TRUE))
 # furrr::future_map(scottraits,
 #                   ~cwm_veg_scot(selectedabun = "zoneA",
 #                             selectedtrait = .,
 #                             selectedtaxres = "stand.spec",
 #                             selectedpolmode = wind,
 #                             selectedpft = allpft),
 #                   .options = furrr_options(seed = TRUE))
 # 
 # # growth form
 # # barely trees in zone A -> only do herbs
 # furrr::future_map(scottraits,
 #                   ~cwm_veg_scot(selectedabun = "zoneA",
 #                             selectedtrait = .,
 #                             selectedtaxres = "stand.spec",
 #                             selectedpolmode = allpol,
 #                             selectedpft = herb
 #                             ),
 #                   .options = furrr_options(seed = TRUE))
 # # Taxonomic resolution
 # furrr::future_map(scottraits,
 #                   ~cwm_veg_scot(selectedabun = "zoneA",
 #                             selectedtrait = .,
 #                             selectedtaxres = "genus",
 #                             selectedpolmode = allpol,
 #                             selectedpft = allpft),
 #                   .options = furrr_options(seed = TRUE))
 # furrr::future_map(scottraits,
 #                   ~cwm_veg_scot(selectedabun = "zoneA",
 #                             selectedtrait = .,
 #                             selectedtaxres = "family",
 #                             selectedpolmode = allpol,
 #                             selectedpft = allpft),
 #                   .options = furrr_options(seed = TRUE))

# Zone B ----
 # all
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneB",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = allpol,
                             selectedpft = allpft,
                             taxtrait = TRUE # save taxon level trait estimates
                   ),
                   .options = furrr_options(seed = TRUE))
 # only taxa in pollen data
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneB",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = allpol,
                             selectedpft = allpft,
                             notinpollen = TRUE),
                   .options = furrr_options(seed = TRUE))
 # pollination
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneB",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = nowind,
                             selectedpft = allpft),
                   .options = furrr_options(seed = TRUE))
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneB",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = wind,
                             selectedpft = allpft),
                   .options = furrr_options(seed = TRUE))
 # growth form
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneB",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = allpol,
                             selectedpft = trsh),
                   .options = furrr_options(seed = TRUE))
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneB",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = allpol,
                             selectedpft = herb),
                   .options = furrr_options(seed = TRUE))
 # Taxonomic resolution
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneB",
                             selectedtrait = .,
                             selectedtaxres = "genus",
                             selectedpolmode = allpol,
                             selectedpft = allpft),
                   .options = furrr_options(seed = TRUE))
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneB",
                             selectedtrait = .,
                             selectedtaxres = "family",
                             selectedpolmode = allpol,
                             selectedpft = allpft),
                   .options = furrr_options(seed = TRUE))
# Zone C ----
 # all
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneC",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = allpol,
                             selectedpft = allpft),
                   .options = furrr_options(seed = TRUE))
 # only taxa in pollen data
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneC",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = allpol,
                             selectedpft = allpft,
                             notinpollen = TRUE),
                   .options = furrr_options(seed = TRUE))
 # pollination
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneC",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = nowind,
                             selectedpft = allpft),
                   .options = furrr_options(seed = TRUE))
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneC",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = wind,
                             selectedpft = allpft),
                   .options = furrr_options(seed = TRUE))
 # growth form
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneC",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = allpol,
                             selectedpft = trsh),
                   .options = furrr_options(seed = TRUE))
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneC",
                             selectedtrait = .,
                             selectedtaxres = "stand.spec",
                             selectedpolmode = allpol,
                             selectedpft = herb
                   ),
                   .options = furrr_options(seed = TRUE))
 # Taxonomic resolution
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneC",
                             selectedtrait = .,
                             selectedtaxres = "genus",
                             selectedpolmode = allpol,
                             selectedpft = allpft),
                   .options = furrr_options(seed = TRUE))
 furrr::future_map(scottraits,
                   ~cwm_veg_scot(selectedabun = "zoneC",
                             selectedtrait = .,
                             selectedtaxres = "family",
                             selectedpolmode = allpol,
                             selectedpft = allpft),
                   .options = furrr_options(seed = TRUE))
}


# Switzerland ----
if (1) { 
  # Load data
  bdm_traits <- readRDS("RDS_files/01_Traits_Swiss.rds") %>% 
    # trees over 200 cm not sampled for traits, use TRY data instead
    filter(!PFT == "TR") %>% 
    filter(!stand.spec == "Unindentified") %>% 
    mutate(stand.spec = str_remove(stand.spec, pattern = " sp\\."))
  bdm_abun <- readRDS("RDS_files/01_Species_abundance_Swiss.rds") %>% 
    mutate(stand.spec = str_remove(stand.spec, pattern = " sp\\."))

  zoneA <- bdm_abun %>%
    arrange(sitename, stand.spec) %>% 
    filter(!stand.spec == "Unindentified") %>% 
    ungroup 
 
  lABUN <- list("zoneA" = zoneA)

  # match taxon names abundance data with trait data
  veg_sp <- lABUN %>% 
    purrr::map(., ~pull(.,stand.spec)) %>% 
    unlist %>% unique 
  bdm_traits <- bdm_traits %>% 
    mutate(., stand.spec =
             case_when(stand.spec %in% veg_sp ~ stand.spec,
                                      genus %in% veg_sp ~ genus,
                                      family %in% veg_sp ~ family)) 
  dfTRY <- readRDS("RDS_files/02_Gapfilled_traits.rds") %>% 
    mutate(., stand.spec = case_when(stand.spec %in% veg_sp ~ stand.spec,
                                     genus %in% veg_sp ~ genus,
                                     family %in% veg_sp ~ family)
    ) %>% 
    drop_na()
  


## Zone A
# all
plan(multisession(workers = 1))
traits <- c("PlantHeight", "LA", "SLA")
furrr::future_map(traits,
                  ~cwm_veg_swiz(selectedabun = "zoneA",
                    selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))

}

