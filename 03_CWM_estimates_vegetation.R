## ---------------------------
##
## Script name:03_CWM_estimates_vegetation.R
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

## Load packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(rjags))  install.packages("rjags")
if(!require(runjags)) install.packages("runjags")
if(!require(mcmc)) install.packages("mcmc")
if(!require(coda)) install.packages("coda")
if(!require(Hmisc)  )install.packages("Hmisc") 
if(!require(furrr))  install.packages("furrr")
if(!require(LCVP)) install.packages("LCVP")
if(!require(lcvplants)) install.packages("lcvplants")

set.seed(1)

## Load and prepare data -----
# Trait data ----
# Scotland ----
plantheight <- readRDS("RDS_files/01_PlantHeight.rds") %>% 
  mutate(sitename = as.factor(sitename), species = as.factor(species),
         PlantHeight = PlantHeight) %>% # plantheight in centimeter
  ungroup()
ldmc <- readRDS("RDS_files/01_LDMC.rds") %>% 
  mutate(sitename = as.factor(sitename), species = as.factor(species)) %>% 
  ungroup()
sla <- readRDS("RDS_files/01_SLA.rds") %>% 
  mutate(sitename = as.factor(sitename), species = as.factor(species)) %>% 
  ungroup()
la <- readRDS("RDS_files/01_LA.rds") %>% 
  mutate(sitename = as.factor(sitename), species = as.factor(species)) %>% 
  ungroup()
c <- readRDS("RDS_files/01_CN.rds") %>% 
  mutate(sitename = as.factor(sitename), species = as.factor(species)) %>% 
  ungroup() %>% 
  rename(LeafC = C)
n <- readRDS("RDS_files/01_CN.rds") %>% 
  mutate(sitename = as.factor(sitename), species = as.factor(species)) %>% 
  ungroup() %>% 
  rename(LeafN = N)

lTrait <- list(PlantHeight = plantheight,
               LDMC = ldmc,
               SLA = sla,
               LA = la,
               LeafC = c,
               LeafN = n)
# Switzerland ----
bdm_traits <- readRDS("RDS_files/01_Traits_Swiss.rds") %>% 
  rename(stand.spec = binom)

# Gapdata ---- 
dfTRAIT <- readRDS("RDS_files/02_Gapfilled_traits.rds")

# Pollination and PFT ----
dfPFT <- readRDS("RDS_files/Polmode_pft_vegetation.rds")

# Species data ----
dfABUN_a <- readRDS("RDS_files/01_Species_abundance_a.rds") 
dfABUN_bc <- readRDS("RDS_files/01_Species_abundance_bc.rds")
bdm_abun <- readRDS("RDS_files/01_Species_abundance_Swiss.rds") %>% 
  # join with pft and polmode data
  left_join(dfPFT[,c("AccSpeciesName", "growthform", "polmode")], by = c("stand.spec" = "AccSpeciesName")) %>% 
  rename(genus = Genus, family = Family)

zoneA <- dfABUN_a %>%
  arrange(sitename, stand.spec) %>% 
  # select plots that are comparable with Switzerland
  filter(distance %in% c("0 meter", "1.5-3 meter")) %>%
  # join with pft and polmode data
  left_join(dfPFT[,c("AccSpeciesName", "growthform", "polmode")], 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  rename(family = fam)
zoneB <- dfABUN_bc %>% 
  # join with pft and polmode data
  left_join(dfPFT[,c("AccSpeciesName", "growthform", "polmode")], 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  ungroup() %>% 
  dplyr::select(sitename, stand.spec, genus, family = fam, abun = spec_abun_b, growthform, polmode)
zoneC <- dfABUN_bc %>% 
  # join with pft and polmode data
  left_join(dfPFT[,c("AccSpeciesName", "growthform", "polmode")], 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  ungroup %>% 
  dplyr::select(sitename, stand.spec, genus, family = fam, abun = spec_abun_c, growthform, polmode)

selectedtrait <- "PlantHeight"
abun <- zoneB
selectedpolmode <- unique(dfPFT$polmode)
selectedpft <- unique(dfPFT$growthform)
selectedtaxres <- "stand.spec"

# Calculate CWM Scotland----
cwm_func <- function(abun, selectedtrait,selectedtaxres, selectedpolmode, selectedpft){
  taxa <- abun %>% 
    arrange(sitename,selectedtaxres) %>% 
    # optional filter for polmode and pft
    filter(polmode %in% selectedpolmode) %>% 
    filter(growthform %in% selectedpft) %>% 
    drop_na() %>% 
    pull(selectedtaxres) %>% 
    unique()
  
  Ab <- abun %>%
    arrange(sitename,selectedtaxres) %>% 
    # optional filter for polmode and pft
    filter(polmode %in% selectedpolmode) %>% 
    filter(growthform %in% selectedpft) %>% 
    rename(tax = all_of(selectedtaxres)) %>% 
    drop_na() %>% 
    # calculate species abundance on the plot level
    group_by(sitename, tax) %>% 
    summarise(abun = sum(abun)) %>% 
    mutate(abun = abun/sum(abun)) %>% 
    filter(!is.nan(abun)) %>%  # filters sites without species of particular selection
    # convert to wide format
    pivot_wider(names_from = tax, values_from = abun,
                values_fill = 0) %>% 
    ungroup() 
  
  # save site names with data
  sitenames <- Ab %>% pull(sitename) 
  # get order of taxa
  taxorder <- colnames(Ab)[-1]
  
  # missing trait data
  traittax <- lTrait %>%
    pluck(selectedtrait) %>%
    arrange(selectedtaxres) %>% 
    pull(selectedtaxres) %>% 
    unique
  missingtax <- taxa[!taxa %in% traittax] 
  
  # filter species with too little observations (2 or less)
  nobs <- lTrait %>%
    pluck(selectedtrait) %>%
    arrange(selectedtaxres) %>% 
    group_by(across(all_of(selectedtaxres))) %>%
    summarise(n = n()) %>%
    filter(n <= 2) %>%
    pull(selectedtaxres)
  
  # find in gapdata
  TRY <- dfTRAIT %>% 
    filter(country == "Scotland") %>%
    dplyr::select(tax = all_of(selectedtaxres), all_of(selectedtrait)) %>% 
    arrange(tax) %>% 
    filter(tax %in% c(missingtax, nobs)) 
    
  Trait <- lTrait %>%
    pluck(selectedtrait) %>%
    dplyr::select(tax = all_of(selectedtaxres), all_of(selectedtrait)) %>% 
    # join with try data for the missing traits
    bind_rows(TRY) %>% 
    filter(tax %in% taxa) %>%
    arrange(match(tax, taxorder)) %>% # order same way as abundance data
    pull(selectedtrait) 
  
  Tax <- lTrait %>%
    pluck(selectedtrait) %>%
    arrange(selectedtaxres) %>% 
    dplyr::select(tax = all_of(selectedtaxres), all_of(selectedtrait)) %>% 
    # join with try data for the missing traits
    bind_rows(TRY) %>% 
    filter(tax %in% taxa) %>%
    arrange(match(tax, taxorder)) %>% # order same way as abundance data
    pull(tax) %>% 
    as.factor()

  # still missing trait data? drop from abundance data
  dropspecies <- taxorder[!taxorder %in% unique(Tax)]
  Ab <- Ab %>% 
    dplyr::select(-sitename, -all_of(dropspecies)) %>% 
    as.matrix()
  
  N <- length(Trait)
  
  Ntax <- ncol(Ab)
  Nsites <- nrow(Ab) 
  
  # for setting priors taxon level
  Mean <- mean(Trait)
  SD <- sd(Trait)
  
  # initial values
  mean.tax <- list(chain1 = rep(min(Trait), Ntax), 
                   chain2 = rep(max(Trait), Ntax))
  tau.tax <- list(chain1 = rep(1/(SD*0.01)^2,Ntax), 
                  chain2 = rep(1/(SD*100)^2, Ntax))
  
  results <- run.jags("CWM_log.txt", n.chains = 2,monitor = "cwm")
  
  # trace plots
  plot(results, file =                          
         paste0("Convergence_diagnostics/Conv_", deparse(substitute(abun)),"_Scotland_", selectedtrait,
                "_", selectedtaxres, "_", selectedpolmode,"_",selectedpft, ".pdf")
  )
  
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_", deparse(substitute(abun)),"_Scotland_", selectedtrait,
             "_", selectedtaxres, "_", selectedpolmode,"_",selectedpft, ".pdf")
  )
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  res <- as.data.frame(summary(results)) %>% 
    mutate(sitename = sitenames)
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_", deparse(substitute(abun)),"_Scotland_", selectedtrait,
                      "_", selectedtaxres, "_", 
                      paste(selectedpolmode,collapse = ""),"_",paste(selectedpft,collapse = ""),".rds"))
}
## Zone A
# all 
plan(multisession(workers = 2))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneA, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(zoneA$polmode),
                            selectedpft = unique(zoneA$growthform)),
                  .options = furrr_options(seed = TRUE))
# pollination
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneA, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = "not wind",
                            selectedpft = unique(zoneA$growthform)),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneA, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = "wind",
                            selectedpft = unique(zoneA$growthform)),
                  .options = furrr_options(seed = TRUE))
# growth form
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneA, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(zoneA$polmode),
                            selectedpft = c("tree", "shrub")),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneA, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(zoneA$polmode),
                            selectedpft = c("grass", "herb", "fern")
                            ),
                  .options = furrr_options(seed = TRUE))
# Taxonomic resolution
plan(multisession(workers = 2))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneA, 
                            selectedtrait = .,
                            selectedtaxres = "genus",
                            selectedpolmode = unique(zoneA$polmode),
                            selectedpft = unique(zoneA$growthform)),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneA, 
                            selectedtrait = .,
                            selectedtaxres = "family",
                            selectedpolmode = unique(zoneA$polmode),
                            selectedpft = unique(zoneA$growthform)),
                  .options = furrr_options(seed = TRUE))

## Zone B
# all 
plan(multisession(workers = 2))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneB, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(zoneB$polmode),
                            selectedpft = unique(zoneB$growthform)),
                  .options = furrr_options(seed = TRUE))
# pollination
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneB, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = "not wind",
                            selectedpft = unique(zoneB$growthform)),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneB, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = "wind",
                            selectedpft = unique(zoneB$growthform)),
                  .options = furrr_options(seed = TRUE))
# growth form
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneB, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(zoneB$polmode),
                            selectedpft = c("tree", "shrub")),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneB, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(zoneB$polmode),
                            selectedpft = c("grass", "herb", "fern")
                  ),
                  .options = furrr_options(seed = TRUE))
# Taxonomic resolution
plan(multisession(workers = 2))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneB, 
                            selectedtrait = .,
                            selectedtaxres = "genus",
                            selectedpolmode = unique(zoneB$polmode),
                            selectedpft = unique(zoneB$growthform)),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneB, 
                            selectedtrait = .,
                            selectedtaxres = "family",
                            selectedpolmode = unique(zoneB$polmode),
                            selectedpft = unique(zoneB$growthform)),
                  .options = furrr_options(seed = TRUE))

## Zone C
# all 
plan(multisession(workers = 2))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneC, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(zoneC$polmode),
                            selectedpft = unique(zoneC$growthform)),
                  .options = furrr_options(seed = TRUE))
# pollination
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneC, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = "not wind",
                            selectedpft = unique(zoneC$growthform)),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneC, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = "wind",
                            selectedpft = unique(zoneC$growthform)),
                  .options = furrr_options(seed = TRUE))
# growth form
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneC, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(zoneC$polmode),
                            selectedpft = c("tree", "shrub")),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneC, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(zoneC$polmode),
                            selectedpft = c("grass", "herb", "fern")
                  ),
                  .options = furrr_options(seed = TRUE))
# Taxonomic resolution
plan(multisession(workers = 2))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneC, 
                            selectedtrait = .,
                            selectedtaxres = "genus",
                            selectedpolmode = unique(zoneC$polmode),
                            selectedpft = unique(zoneC$growthform)),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func(abun = zoneC, 
                            selectedtrait = .,
                            selectedtaxres = "family",
                            selectedpolmode = unique(zoneC$polmode),
                            selectedpft = unique(zoneC$growthform)),
                  .options = furrr_options(seed = TRUE))



# Switzerland ----
cwm_func_swiz <- function(abun, selectedtrait,selectedtaxres, selectedpolmode, selectedpft){
  taxa <- abun %>% 
    arrange(sitename,selectedtaxres) %>% 
    # optional filter for polmode and pft
    filter(polmode %in% selectedpolmode) %>% 
    filter(growthform %in% selectedpft) %>% 
    drop_na(-Notes) %>% 
    pull(selectedtaxres) %>% 
    unique()
  
  Ab <- abun %>%
    arrange(sitename,selectedtaxres) %>% 
    # optional filter for polmode and pft
    filter(polmode %in% selectedpolmode) %>% 
    filter(growthform %in% selectedpft) %>% 
    rename(tax = all_of(selectedtaxres)) %>%
    drop_na(-Notes) %>% 
    # calculate species abundance on the plot level
    group_by(sitename, tax) %>% 
    summarise(abun = sum(abun)) %>% 
    mutate(abun = abun/sum(abun)) %>% 
    filter(!is.nan(abun)) %>%  # filters sites without species of particular selection
    # convert to wide format
    pivot_wider(names_from = tax, values_from = abun,
                values_fill = 0) %>% 
    ungroup() 
  
  # save site names with data
  sitenames <- Ab %>% pull(sitename) 
  # get order of taxa
  taxorder <- colnames(Ab)[-1]
 
  # missing trait data
  traittax <- bdm_traits %>%
    arrange(selectedtaxres) %>% 
    pull(selectedtaxres) %>% 
    unique
  missingtax <- taxa[!taxa %in% traittax] 
  
  # filter species with too little observations (2 or less)
  nobs <- bdm_traits %>%
    arrange(selectedtaxres) %>%
    group_by(across(all_of(selectedtaxres))) %>%
    summarise(n = n()) %>%
    filter(n <= 2) %>%
    pull(selectedtaxres)

  # find in gapdata
  TRY <- dfTRAIT %>% 
    filter(country == "Switzerland") %>%
    dplyr::select(tax = all_of(selectedtaxres), all_of(selectedtrait)) %>% 
    arrange(tax) %>% 
    filter(tax %in% c(missingtax,nobs)) 
  
  Trait <- bdm_traits %>%
    drop_na(all_of(selectedtrait)) %>% 
    dplyr::select(tax = all_of(selectedtaxres), all_of(selectedtrait)) %>% 
    # join with try data for the missing traits
    bind_rows(TRY) %>% 
    filter(tax %in% taxa) %>%
    arrange(match(tax, taxorder)) %>% # order same way as abundance data
    pull(selectedtrait) 
  
  Tax <- bdm_traits %>%
    drop_na(all_of(selectedtrait)) %>% 
    dplyr::select(tax = all_of(selectedtaxres), all_of(selectedtrait)) %>% 
    # join with try data for the missing traits
    bind_rows(TRY) %>% 
    filter(tax %in% taxa) %>%
    arrange(match(tax, taxorder)) %>% # order same way as abundance data
    pull(tax) %>% 
    as.factor()
  
  # still missing trait data? drop from abundance data
  dropspecies <- taxorder[!taxorder %in% unique(Tax)]
  Ab <- Ab %>% 
    dplyr::select(-sitename, -all_of(dropspecies)) %>% 
    as.matrix()
  
  N <- length(Trait)
  
  Ntax <- ncol(Ab)
  Nsites <- nrow(Ab) 
  
  # for setting priors taxon level
  Mean <- mean(Trait)
  SD <- sd(Trait)
  
  # initial values
  mean.tax <- list(chain1 = rep(min(Trait), Ntax), 
                   chain2 = rep(max(Trait), Ntax))
  tau.tax <- list(chain1 = rep(1/(SD*0.01)^2,Ntax), 
                  chain2 = rep(1/(SD*100)^2, Ntax))
  
  results <- run.jags("CWM_log.txt", n.chains = 2,monitor = "cwm")
  
  # trace plots
  plot(results, file =                          
         paste0("Convergence_diagnostics/Conv_", deparse(substitute(abun)),"_Switzerland_", selectedtrait,
                "_", selectedtaxres, "_", selectedpolmode,"_",selectedpft, ".pdf")
  )
  
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_", deparse(substitute(abun)),"_Switzerland_", selectedtrait,
             "_", selectedtaxres, "_", selectedpolmode,"_",selectedpft, ".pdf")
  )
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  res <- as.data.frame(summary(results)) %>% 
    mutate(sitename = sitenames)
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_", deparse(substitute(abun)),"_Switzerland_", selectedtrait,
                      "_", selectedtaxres, "_", 
                      paste(selectedpolmode,collapse = ""),"_",paste(selectedpft,collapse = ""),".rds"))
}
## Zone A
# all 
plan(multisession(workers = 2))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func_swiz(abun = bdm_abun, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(bdm_abun$polmode),
                            selectedpft = unique(bdm_abun$growthform)),
                  .options = furrr_options(seed = TRUE))
# pollination
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func_swiz(abun = bdm_abun, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = "not wind",
                            selectedpft = unique(bdm_abun$growthform)),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func_swiz(abun = bdm_abun, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = "wind",
                            selectedpft = unique(bdm_abun$growthform)),
                  .options = furrr_options(seed = TRUE))
# growth form
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func_swiz(abun = bdm_abun, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(bdm_abun$polmode),
                            selectedpft = c("tree", "shrub")),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func_swiz(abun = bdm_abun, 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = unique(bdm_abun$polmode),
                            selectedpft = c("grass", "herb", "fern")
                  ),
                  .options = furrr_options(seed = TRUE))
# Taxonomic resolution
plan(multisession(workers = 2))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func_swiz(abun = bdm_abun, 
                            selectedtrait = .,
                            selectedtaxres = "genus",
                            selectedpolmode = unique(bdm_abun$polmode),
                            selectedpft = unique(bdm_abun$growthform)),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(c("PlantHeight", "LA", "SLA"),
                  ~cwm_func_swiz(abun = bdm_abun, 
                            selectedtrait = .,
                            selectedtaxres = "family",
                            selectedpolmode = unique(bdm_abun$polmode),
                            selectedpft = unique(bdm_abun$growthform)),
                  .options = furrr_options(seed = TRUE))



