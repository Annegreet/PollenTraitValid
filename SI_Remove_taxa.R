## ---------------------------
##
## Script name: SI_Remove_taxa
##
## Purpose of script: Testing the effect of a single species or taxon on the CMW
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

# Species data ----
dfABUN_a <- readRDS("RDS_files/01_Species_abundance_a.rds") 
dfABUN_bc <- readRDS("RDS_files/01_Species_abundance_bc.rds")

spec <- c(dfABUN_a$stand.spec, dfABUN_bc$stand.spec) %>% unique()

# use only species that are observed in the field?
# dfTRAIT <- dfTRAIT %>% 
#   filter(stand.spec %in% spec)

selectedcountry <- "Scotland"
selectedtrait <- "LA"
selectedabun <- str_subset(colnames(dfPOL), "percent")[1]
removedspec <- unique(dfPOL$pollentaxon)[3]
normaltraits <- c("PlantHeight","SLA","LA")

# Pollen  ----
labs_trait <- c("Height (cm)(log)",
                "Leaf area (mm2)(log)",
                "SLA (mm2/mg)(log)")
names(labs_trait) <- c("PlantHeight", "LA", "SLA")

# Run model
plan(multisession(workers = 2))
furrr::future_map(unique(dfPOL_Scot$pollentaxon),
                    ~cwm_pol_remove(selectedtrait = "LA",
                                selectedcountry = "Scotland",
                                selectedabun = "percent",
                                removedspec = .x),
                  .options = furrr_options(seed = TRUE))
# furrr::future_map(unique(dfPOL_Scot$pollentaxon),
#                   ~cwm_pol_remove(selectedtrait = "SLA",
#                               selectedcountry = "Scotland",
#                               selectedabun = "percent",
#                               removedspec = .x),
#                   .options = furrr_options(seed = TRUE))   
furrr::future_map(unique(dfPOL_Scot$pollentaxon),
                  ~cwm_pol_remove(selectedtrait = "PlantHeight",
                              selectedcountry = "Scotland",
                              selectedabun = "percent",
                              removedspec = .x),
                  .options = furrr_options(seed = TRUE))  


# Compile estimates
pollen_files <- 
  list.files("RDS_files/",full.names = TRUE) %>% 
  str_subset(pattern = "pollen") %>% 
  str_subset(pattern = paste(unique(dfPOL_Scot$pollentaxon), collapse = "|")) %>% 
  str_subset(pattern = "reduced_sp", negate = T)

pollen_files_all <-  list.files("RDS_files/",full.names = TRUE) %>% 
  str_subset(pattern = "03_CWM_estimates_pollen_Scotland_LA_percent.rds|03_CWM_estimates_pollen_Scotland_SLA_percent.rds|03_CWM_estimates_pollen_Scotland_PlantHeight_percent.rds")

dfCWM_pol <- c(pollen_files, pollen_files_all) %>% 
  purrr::map2(., str_remove(c(pollen_files, pollen_files_all), "RDS_files/"), ~readRDS(.x) %>% 
                as.data.frame())  %>% 
  bind_rows()  %>% 
  mutate(removed_spec = ifelse(is.na(removed_spec), "All species", removed_spec))
p <- dfCWM_pol %>% 
  filter(removed_spec %in% c("Pinus", "Betula", "Poaceae","Asteraceae",
                             "Ericales (tetrad)", "Pteridophyte",
                             "All species")) %>%
  mutate(trait = factor(trait, level = c("LA", "PlantHeight", "SLA"),
                      labels = c(LA = "Leaf~area~(mm^{2})(log)",
                                 PlantHeight = "Height~(cm)(log)",
                                 SLA = "SLA~(mm^{2}/mg)(log)"))) %>% 
  ggplot(aes(x = Mean, y = removed_spec)) +
  geom_boxplot(fill = "grey") +
  ylab("") + 
  xlab("") +
  facet_grid(~trait, labeller = label_parsed, scales = "free") +
  theme_bw()  +
  theme(legend.position = "none")
ggsave("Figures/CWM_pollen_removed_taxa.png", height = 4, width = 7)

# Vegetation ----
plantheight <- readRDS("RDS_files/01_PlantHeight.rds") %>% 
  ungroup() 
ldmc <- readRDS("RDS_files/01_LDMC.rds") %>% 
  ungroup()
sla <- readRDS("RDS_files/01_SLA.rds") %>% 
  ungroup()
la <- readRDS("RDS_files/01_LA.rds") %>% 
  ungroup()


lTrait <- list(PlantHeight = plantheight,
               LDMC = ldmc,
               SLA = sla,
               LA = la)

dfTRY_raw <- readRDS("RDS_files/02_Gapfilled_traits.rds")

# Pollination and PFT 
dfPFT <- readRDS("RDS_files/Polmode_pft_vegetation.rds")

# Species to pollen table
polspec <- readRDS("RDS_files/02_PollenType_species.rds") %>% 
  dplyr::select(pollentaxon, stand.spec) %>% distinct() %>% as_tibble()

# Species data ----
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
allpft  <- c("tree","herb","grass","fern","shrub", NA)
wind <- "wind"
nowind <- "not wind"
allpol <- c("wind", "not wind", NA)

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
  )

selectedtrait <- "LA"
selectedabun <- "zoneB"
selectedtaxres <- "stand.spec"
selectedpolmode <- allpol
selectedpft <- allpft
notinpollen <- TRUE

scottraits <- c("PlantHeight", "LA", "SLA")

# function

plan(multisession(workers = 4))
zoneB %>%
  group_by(sitename, stand.spec) %>% 
  summarise(abun = sum(abun)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  filter(!is.nan(abun)) %>% 
  pull(stand.spec) %>% 
  unique() %>% 
  map(.,
      ~cwm_veg_remove(selectedabun = "zoneB",
                selectedtrait = "LA",
                removedspec = .,
                selectedtaxres = "stand.spec",
                selectedpolmode = allpol,
                selectedpft = allpft))
zoneB %>%
  group_by(sitename, stand.spec) %>% 
  summarise(abun = sum(abun)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  filter(!is.nan(abun)) %>% 
  filter(!stand.spec %in% spec_done) %>% 
  pull(stand.spec) %>% 
  unique() %>% 
  map(.,
      ~cwm_veg_remove(selectedabun = "zoneB",
                selectedtrait = "PlantHeight",
                removedspec = .,
                selectedtaxres = "stand.spec",
                selectedpolmode = allpol,
                selectedpft = allpft))
zoneB %>%
  group_by(sitename, stand.spec) %>% 
  summarise(abun = sum(abun)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  filter(!is.nan(abun)) %>% 
  filter(!stand.spec %in% spec_done) %>% 
  pull(stand.spec) %>% 
  unique() %>% 
  map(.,
      ~cwm_veg_remove(selectedabun = "zoneB",
                selectedtrait = "SLA",
                removedspec = .,
                selectedtaxres = "stand.spec",
                selectedpolmode = allpol,
                selectedpft = allpft))
veg_files <- 
  list.files("RDS_files/",full.names = TRUE) %>% 
  str_subset(pattern = paste(unique(zoneB$stand.spec), collapse = "|")) %>% 
  str_subset(pattern = "zoneB")

dfCWM_veg <- veg_files %>% 
  purrr::map2(., str_remove(veg_files, "RDS_files/"), ~readRDS(.x) %>% 
                as.data.frame()) %>% 
  bind_rows() 

ggplot(dfCWM_veg, aes(x = Mean, y = removed_spec)) +
  geom_boxplot(fill = "grey") +
  ylab("") + 
  xlab("") +
  facet_grid(~trait, labeller = labeller(trait = labs_trait)) +
  theme_bw()  +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"))
