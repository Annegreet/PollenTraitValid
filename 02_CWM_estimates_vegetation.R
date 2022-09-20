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

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rjags))  install.packages("rjags")
if (!require(runjags)) install.packages("runjags")
if (!require(mcmc)) install.packages("mcmc")
if (!require(coda)) install.packages("coda")
if (!require(furrr))  install.packages("furrr")

set.seed(1)

## Load and prepare data -----
# Trait data ----
# Scotland ----
plantheight <- readRDS("RDS_files/01_PlantHeight.rds") %>% 
  ungroup() 
ldmc <- readRDS("RDS_files/01_LDMC.rds") %>% 
  ungroup()
sla <- readRDS("RDS_files/01_SLA.rds") %>% 
  ungroup()
la <- readRDS("RDS_files/01_LA.rds") %>% 
  ungroup()
c <- readRDS("RDS_files/01_CN.rds") %>% 
  ungroup() %>% 
  rename(LeafC = C)
n <- readRDS("RDS_files/01_CN.rds") %>% 
  ungroup() %>% 
  rename(LeafN = N)

lTrait <- list(PlantHeight = plantheight,
               LDMC = ldmc,
               SLA = sla,
               LA = la,
               LeafC = c,
               LeafN = n)

dfTRY <- readRDS("RDS_files/01_Clean_TRY_data.rds")

# Pollination and PFT ----
dfPFT <- readRDS("RDS_files/Polmode_pft_vegetation.rds")

# Species data ----
dfABUN_a <- readRDS("RDS_files/01_Species_abundance_a.rds") 
dfABUN_bc <- readRDS("RDS_files/01_Species_abundance_bc.rds")

zoneA <- dfABUN_a %>%
  arrange(sitename, stand.spec) %>% 
  # select plots that are comparable with Switzerland
  filter(distance %in% c("0 meter", "1.5-3 meter")) %>%
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  rename(family = fam) %>% 
  ungroup
zoneB <- dfABUN_bc %>% 
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  ungroup() %>% 
  dplyr::select(sitename, stand.spec, genus, family = fam, abun = spec_abun_b, growthform, polmode)
zoneC <- dfABUN_bc %>% 
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  ungroup %>% 
  dplyr::select(sitename, stand.spec, genus, family = fam, abun = spec_abun_c, growthform, polmode)
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
  dfTRY %>% 
  mutate(., stand.spec = case_when(stand.spec %in% veg_sp ~ stand.spec,
                                   genus %in% veg_sp ~ genus,
                                   family %in% veg_sp ~ family)
                        )
             
# Scotland----
if (0) { 
cwm_func <- function(selectedabun, selectedtrait, selectedtaxres, selectedpolmode, 
                     selectedpft){
  abun <- lABUN %>% 
    pluck(selectedabun) %>% 
    # optional filter for polmode and pft
    filter(polmode %in% selectedpolmode) %>% 
    filter(growthform %in% selectedpft) %>% 
    # rename column to selected taxonomic resolution
    rename(tax = all_of(selectedtaxres)) %>% 
    filter(!tax == "Unindentified") %>% 
    dplyr::select(sitename, tax, abun) %>% 
    drop_na() 
  
  ## prepare trait data
  # taxa in the abundance data
  taxa <- abun %>% 
    pull(tax) %>% 
    unique() %>% 
    sort
  
  # missing field collected trait data
  traittax <- lTrait %>%
    pluck(selectedtrait) %>% 
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
  TRY <- dfTRY %>% 
    dplyr::select(tax = all_of(selectedtaxres), 
                  all_of(selectedtrait)) %>% 
    filter(tax %in% c(missingtax,nobs)) 
  
  trait <- lTrait %>%
    pluck(selectedtrait) %>% 
    dplyr::select(tax = all_of(selectedtaxres), 
                  all_of(selectedtrait)) %>% 
    # join with try data for the missing traits
    bind_rows(TRY) %>% 
    filter(tax %in% taxa) %>%
    # sort alphabetically
    arrange(tax) %>%
    # drop NA's and 0's
    mutate(across(where(is.double), ~na_if(.,0))) %>%
    drop_na() %>% 
    # check number of observations
    group_by(tax) %>% 
    mutate(ntrait = n()) %>% 
    filter(ntrait > 3)
  
  ## Model data (capitalized)
  Trait <- trait %>% 
    pull(selectedtrait)
  
  Tax <- trait  %>% 
    pull(tax) %>% 
    as.factor()
  
    # Prepare abundance matrix
  Ab <- abun %>%
    # filter species without enough data (more than 2)
    filter(tax %in% unique(Tax)) %>% 
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
  
  Ab <- Ab %>% 
    dplyr::select(-sitename) %>% 
    # sort alphabetically
    select(order(colnames(.))) %>% 
    as.matrix() 
  
  N <- length(Trait)
  
  Ntax <- ncol(Ab)
  Nsites <- nrow(Ab) 
  
  # for setting priors taxon level
  MeanLog <- aggregate(Trait, list(Tax), FUN = mean)[,"x"] %>% log()
  SDLog <- aggregate(Trait, list(Tax), FUN = sd)[,"x"] %>%
    log() %>% 
    na_if(., 0) # prior cant take 0 or negative
  SDLog[SDLog < 0] <- NA  
  SDLog <- replace_na(SDLog, mean(SDLog, na.rm = TRUE))  # in case of only 1 obs set sd to mean sd
  
  # initial values
  mean.tax <- list(chain1 = MeanLog*0.01, chain2 = MeanLog*100)
  tau.tax <- list(chain1 = 1/(SDLog*0.01)^2, chain2 = 1/(SDLog*100)^2)
  cwm <- list(chain1 = rep(min(Trait), Nsites), 
              chain2 = rep(max(Trait), Nsites))
  
  # run model
  results <- run.jags("CWM_log.txt", n.chains = 2, monitor = "cwm",
                      burnin = 10000, sample = 15000)
  
  # trace plots
  plot(results, file =                          
         paste0("Convergence_diagnostics/Conv_", selectedabun,"_Scotland_", selectedtrait,
                "_", selectedtaxres, "_", paste(selectedpolmode, collapse = ""),"_",paste(selectedpft, collapse = ""), ".pdf")
  )
  
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_", selectedabun,"_Scotland_", selectedtrait,
             "_", selectedtaxres, "_", paste(selectedpolmode, collapse = ""),"_",paste(selectedpft, collapse = ""), ".pdf")
  )
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  res <- as.data.frame(summary(results)) %>% 
    mutate(sitename = sitenames,
    trait = selectedtrait,
    growthform = case_when(all(selectedpft == trsh) ~ "trsh",
                          all(selectedpft == herb) ~ "herb",
                          TRUE ~ "allpft"),
    pollination = case_when(all(selectedpolmode == wind) ~ "wind",
                            all(selectedpolmode == nowind) ~ "no wind",
                            TRUE ~ "allpol"),
    taxres = selectedtaxres,
    zone = selectedabun)
  
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_",selectedabun,"_Scotland_", selectedtrait,
                      "_", selectedtaxres, "_", 
                      paste(selectedpolmode, collapse = ""),"_",paste(selectedpft,collapse = ""),".rds"))
}

## Zone A
# all 
scottraits <- c("PlantHeight", "LA", "SLA")
plan(multisession(workers = 4))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneA",
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# pollination
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneA",
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = nowind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneA",
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = wind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))

# growth form
# barely trees in zone A -> only do herbs
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneA",
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = herb
                            ),
                  .options = furrr_options(seed = TRUE))
# Taxonomic resolution
plan(multisession(workers = 4))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneA", 
                            selectedtrait = .,
                            selectedtaxres = "genus",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneA", 
                            selectedtrait = .,
                            selectedtaxres = "family",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))

## Zone B
# all 
plan(multisession(workers = 4))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# pollination
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = nowind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = wind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# growth form
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = trsh),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = herb
                  ),
                  .options = furrr_options(seed = TRUE))
# Taxonomic resolution
plan(multisession(workers = 4))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "genus",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "family",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))

## Zone C
# all 
plan(multisession(workers = 4))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# pollination
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = nowind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = wind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# growth form
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft =trsh),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = herb
                  ),
                  .options = furrr_options(seed = TRUE))
# Taxonomic resolution
plan(multisession(workers = 4))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "genus",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(scottraits,
                  ~cwm_func(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "family",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
}


# Switzerland ----
if (1) { 
  bdm_traits <- readRDS("RDS_files/01_Traits_Swiss.rds")
  bdm_abun <- readRDS("RDS_files/01_Species_abundance_Swiss.rds") 
  bdm_abun_bc <- readRDS("RDS_files/01_Species_abundance_zoneBC_swiss.rds")
  bdm_pollen <- readRDS("RDS_files/01_Pollen_data_Swiss.rds")
  pollen_ids <- bdm_pollen$sitename %>% unique %>% str_remove("X")
  
  dfPFT <- readRDS("RDS_files/Polmode_pft_vegetation.rds")
  
  zoneA <- bdm_abun %>%
    arrange(sitename, stand.spec) %>% 
    # join with pft and polmode data
    left_join(dfPFT, 
              by = c("stand.spec" = "AccSpeciesName")) %>% 
    ungroup %>% 
    filter(sitename %in% pollen_ids)
  zoneB <- bdm_abun_bc %>% 
    # join with pft and polmode data
    left_join(dfPFT, 
              by = c("stand.spec" = "AccSpeciesName")) %>% 
    ungroup() %>% 
    dplyr::select(sitename, stand.spec, genus, family, abun = spec_abun_b, growthform, polmode) %>% 
    filter(sitename %in% pollen_ids)
  zoneC <- bdm_abun_bc %>% 
    # join with pft and polmode data
    left_join(dfPFT, 
              by = c("stand.spec" = "AccSpeciesName")) %>% 
    ungroup %>% 
    dplyr::select(sitename, stand.spec, genus, family, abun = spec_abun_c, growthform, polmode) %>% 
    filter(sitename %in% pollen_ids)
  lABUN <- list("zoneA" = zoneA, "zoneB" = zoneB, "zoneC" = zoneC)

  # match taxon names abundance data with trait data
  veg_sp <- lABUN %>% 
    purrr::map(., ~pull(.,stand.spec)) %>% 
    unlist %>% unique
  bdm_traits <- bdm_traits %>% 
    mutate(., stand.spec =
             case_when(stand.spec %in% veg_sp ~ stand.spec,
                                      genus %in% veg_sp ~ genus,
                                      family %in% veg_sp ~ family)
    ) 
  dfTRY <- readRDS("RDS_files/01_Clean_TRY_data.rds") %>% 
    mutate(., stand.spec = case_when(stand.spec %in% veg_sp ~ stand.spec,
                                     genus %in% veg_sp ~ genus,
                                     family %in% veg_sp ~ family)
    )
  
cwm_func_swiz <- function(selectedabun, selectedtrait, selectedtaxres, selectedpolmode, selectedpft){
  ## Select abundance data
    abun <- lABUN %>% 
      pluck(selectedabun) %>% 
      # optional filter for polmode and pft
      filter(polmode %in% selectedpolmode) %>% 
      filter(growthform %in% selectedpft) %>% 
      # rename column to selected taxonomic resolution
      rename(tax = all_of(selectedtaxres)) %>% 
      filter(!tax == "Unindentified") %>% 
      dplyr::select(sitename, tax, abun) %>% 
      drop_na() 
    
    ## prepare trait data
    # taxa in the abundance data
    taxa <- abun %>% 
      pull(tax) %>% 
      unique() %>% 
      sort
    
    # missing field collected trait data
    traittax <- bdm_traits %>%
      pull(selectedtaxres) %>% 
      unique
    missingtax <- taxa[!taxa %in% traittax] 
    
    # find species with too little observations (2 or less)
    nobs <- bdm_traits %>%
      arrange(selectedtaxres) %>%
      group_by(across(all_of(selectedtaxres))) %>%
      summarise(n = n()) %>%
      filter(n <= 2) %>%
      pull(selectedtaxres)
     
    # find in missing species in  try data
    TRY <- dfTRY %>% 
      dplyr::select(tax = all_of(selectedtaxres), all_of(selectedtrait)) %>% 
      arrange(tax) %>% 
      filter(tax %in% c(missingtax,nobs)) 
    
    trait <- bdm_traits %>%
      dplyr::select(tax = all_of(selectedtaxres), all_of(selectedtrait)) %>% 
      # join with try data for the missing traits
      bind_rows(TRY) %>% 
      filter(tax %in% taxa) %>%
      # sort alphabetically
      arrange(tax) %>%
      # drop NA's and 0's
      mutate(across(where(is.double), ~na_if(.,0))) %>%
      drop_na() %>% 
      # check number of observations
      group_by(tax) %>% 
      mutate(ntrait = n()) %>% 
      filter(ntrait > 3)
    
    ## Model data (capitalized)
    Trait <- trait %>% 
      pull(selectedtrait)
    
    Tax <- trait  %>% 
      pull(tax) %>% 
      as.factor()
    
    # Prepare abundance matrix
    Ab <- abun %>%
      # filter species without enough data (more than 2)
      filter(tax %in% unique(Tax)) %>% 
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
    
    Ab <- Ab %>% 
      dplyr::select(-sitename) %>% 
      # sort alphabetically
      select(order(colnames(.))) %>% 
      as.matrix() 
    
    N <- length(Trait)
    
    Ntax <- ncol(Ab)
    Nsites <- nrow(Ab) 
    
    # for setting priors taxon level
    MeanLog <- aggregate(Trait, list(Tax), FUN = mean)[,"x"] %>% log()
    SDLog <- aggregate(Trait, list(Tax), FUN = sd)[,"x"] %>%
      log() %>% 
      na_if(., 0) # prior cant take 0 or negative
    SDLog[SDLog < 0] <- NA  
    SDLog <- replace_na(SDLog, mean(SDLog, na.rm = TRUE))  # in case of only 1 obs set sd to mean sd
  
    # initial values
    mean.tax <- list(chain1 = MeanLog*0.01, chain2 = MeanLog*100)
    tau.tax <- list(chain1 = 1/(SDLog*0.01)^2, chain2 = 1/(SDLog*100)^2)
    cwm <- list(chain1 = rep(min(Trait), Nsites), 
                chain2 = rep(max(Trait), Nsites))
    
    # run model
    results <- run.jags("CWM_log.txt", n.chains = 2, monitor = "cwm",
                        burnin = 10000, sample = 15000)
    res <- summary(results)
    
  # trace plots
  plot(results, file =                          
         paste0("Convergence_diagnostics/Conv_Switzerland_",selectedabun,"_", selectedtrait,
                "_", selectedtaxres, "_", paste(selectedpolmode, collapse = ""),"_",
                paste(selectedpft, collapse = ""), ".pdf")
  )
  
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_Switzerland_", selectedtrait,
             "_",selectedabun,"_", selectedtaxres, "_", paste(selectedpolmode, collapse = ""),"_",
             paste(selectedpft, collapse = ""), ".pdf")
  )
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  res <- as.data.frame(summary(results)) %>% 
    mutate(sitename = sitenames,
           trait = selectedtrait,
           growthform = case_when(all(selectedpft == trsh) ~ "trsh",
                                  all(selectedpft == herb) ~ "herb",
                                  TRUE ~ "allpft"),
           pollination = case_when(all(selectedpolmode == wind) ~ "wind",
                                   all(selectedpolmode == nowind) ~ "no wind",
                                   TRUE ~ "allpol"),
           taxres = selectedtaxres,
           zone = selectedabun)
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_Switzerland_", selectedtrait,
                      "_",selectedabun,"_", selectedtaxres, "_", 
                      paste(selectedpolmode,collapse = ""),"_",paste(selectedpft,collapse = ""),".rds"))
}

## Zone A
# all
traits <- c("PlantHeight", "LA", "SLA")[1]
plan(multisession(workers = 4))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneA",
                    selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# pollination
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneA",
                                 selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = nowind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneA",
                                 selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = wind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# growth form
# barely trees in zone A -> only run herbs
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneA",
                                 selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = herb
                  ),
                  .options = furrr_options(seed = TRUE))
# Taxonomic resolution
plan(multisession(workers = 4))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneA",
                                 selectedtrait = .,
                            selectedtaxres = "genus",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneA",
                                 selectedtrait = .,
                            selectedtaxres = "family",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# Zone B
# all 
plan(multisession(workers = 4))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# pollination
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = nowind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = wind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# growth form
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = trsh),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = herb
                  ),
                  .options = furrr_options(seed = TRUE))
# Taxonomic resolution
plan(multisession(workers = 4))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "genus",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneB", 
                            selectedtrait = .,
                            selectedtaxres = "family",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))

## Zone C
# all 
plan(multisession(workers = 4))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# pollination
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = nowind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = wind,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
# growth form
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft =trsh),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "stand.spec",
                            selectedpolmode = allpol,
                            selectedpft = herb
                  ),
                  .options = furrr_options(seed = TRUE))
# Taxonomic resolution
plan(multisession(workers = 4))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "genus",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(traits,
                  ~cwm_func_swiz(selectedabun = "zoneC", 
                            selectedtrait = .,
                            selectedtaxres = "family",
                            selectedpolmode = allpol,
                            selectedpft = allpft),
                  .options = furrr_options(seed = TRUE))

}


