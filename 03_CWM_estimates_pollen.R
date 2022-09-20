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
dfTRAIT <- readRDS("RDS_files/01_TRY_raw_traits.rds")

# Pollen data ----
dfPOL_Scot <- readRDS("RDS_files/01_pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_pollen_data_Swiss.rds")
dfPOL <- bind_rows("Scotland" = dfPOL_Scot, "Switzerland" =
                     dfPOL_Swiss, .id = "country") %>% 
  ungroup()

selectedcountry <- "Scotland"
selectedtrait <- "LA"
selectedabun <- str_subset(colnames(dfPOL), "percent")[1]

normaltraits <- c("PlantHeight","SLA","LA")

# No subset ----
if (0) {
  cwmlogfunc <- function(selectedtrait, selectedcountry, selectedabun){
    trait <- dfTRAIT %>% 
      filter(country == selectedcountry) %>% 
      # drop NA's and 0's
      dplyr::select(pollentaxon, all_of(selectedtrait)) %>% 
      mutate(across(where(is.double), ~na_if(.,0))) %>% 
      drop_na() %>% 
      # Sort alphabetically
      arrange(pollentaxon) %>% 
      # check number of observations
      group_by(pollentaxon) %>% 
      mutate(ntrait = n()) %>% 
      filter(ntrait > 3)
    
    pol <- dfPOL %>% 
      # select country
      filter(country == selectedcountry) %>%
      dplyr::select(sitename, pollentaxon, pollen_abun = all_of(selectedabun)) %>% 
      drop_na()
    
    # save pollen taxa
    pollentaxa <- pol %>% pull(pollentaxon) %>% unique()
    
    ## Model data (capitalized)
    Trait <- trait %>% 
      filter(pollentaxon %in% pollentaxa) %>% 
      pull(selectedtrait) 
    
    Tax <- trait %>% 
      filter(pollentaxon %in% pollentaxa) %>% 
      pull(pollentaxon) %>% 
      as.factor() 
    
    N <- length(Trait)
    
    Ab <- pol %>%
      # this filters taxa with too little trait data
      filter(pollentaxon %in% unique(Tax)) %>% 
      # calculate percentage after removing taxa
      group_by(sitename, pollentaxon) %>% 
      summarise(abun = sum(pollen_abun)) %>% 
      mutate(abun = abun/sum(abun)) %>% 
      # convert to wide format
      pivot_wider(names_from = pollentaxon, values_from = abun) %>% 
      # filter empty columns (taxon not in that country)
      discard(~all(is.na(.))) %>% 
      # replace remaining nas by 0
      replace(is.na(.), 0) %>% 
      # sort alphabetically
      ungroup()
    # save site names with data
    sitenames <- Ab %>% pull(sitename)
    Ab <- Ab %>% 
      dplyr::select(-sitename, sort(colnames(.[-1]))) %>% 
      as.matrix() 
    
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
    plot(results, file = paste0("Convergence_diagnostics/Conv_pollen_", selectedcountry, "_",
                                selectedtrait, "_", selectedabun, ".pdf"))
    # gelman plots
    pdf(paste0("Convergence_diagnostics/Gelman_pollen_", selectedcountry,
               "_", selectedtrait,"_", selectedabun, ".pdf"))
    
    gelman.plot(results, ask = FALSE)
    dev.off()
    
    # summarize results
    res <- data.frame(summary(results),
                      sitename = sitenames,
                      trait = selectedtrait,
                      country = selectedcountry,
                      correction = selectedabun)
    saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", selectedcountry,"_",
                        selectedtrait, "_", selectedabun, ".rds"))
  }
  if (Sys.info()[4] == "LUIN67311") {
  plan(multisession, workers = 1)
    } else {
  plan(multisession, workers = 6)
  }
  furrr::future_map(normaltraits,
                     ~cwmlogfunc(selectedtrait = .x,
                                 selectedcountry = "Scotland",
                                 selectedabun = "percent"))
  furrr::future_map(normaltraits,
                     ~cwmlogfunc(selectedtrait = .x,
                                 selectedcountry = "Switzerland",
                                 selectedabun = "percent"))
  furrr::future_map(normaltraits,
                    ~cwmlogfunc(selectedtrait = .x,
                                selectedcountry = "Scotland",
                                selectedabun = "adjustedpercent_mean"))
  furrr::future_map(normaltraits,
                    ~cwmlogfunc(selectedtrait = .x, 
                                selectedcountry = "Switzerland",
                                selectedabun = "adjustedpercent_mean"))
  
}


# Growth form ----
if (0) {
cwmlogfunc <- function(selectedtrait, selectedcountry, selectedabun,
                       selectedpft){
  trait <- dfTRAIT %>% 
    filter(country == selectedcountry) %>% 
    # drop NA's and 0's
    dplyr::select(pollentaxon, all_of(selectedtrait)) %>% 
    mutate(across(where(is.double), ~na_if(.,0))) %>% 
    drop_na() %>% 
    # Sort alphabetically
    arrange(pollentaxon) %>% 
    # check number of observations
    group_by(pollentaxon) %>% 
    mutate(ntrait = n()) %>% 
    filter(ntrait > 3)
  
  pol <- dfPOL %>% 
    # select country
    filter(country == selectedcountry) %>%
    filter(growthform %in% selectedpft) %>%
    dplyr::select(sitename, pollentaxon, 
                  pollen_abun = percent) %>% 
    drop_na()
  
   # save pollen taxa
  pollentaxa <- pol %>% pull(pollentaxon) %>% unique()
  
  ## Model data (capitalized)
  Trait <- trait %>% 
    filter(pollentaxon %in% pollentaxa) %>% 
    pull(selectedtrait) 
  
  Tax <- trait %>% 
    filter(pollentaxon %in% pollentaxa) %>% 
    pull(pollentaxon) %>% 
    as.factor() 
  
  N <- length(Trait)
  
  Ab <- pol %>%
    # this filters taxa with too little trait data
    filter(pollentaxon %in% unique(Tax)) %>% 
    # calculate percentage after removing taxa
    group_by(sitename, pollentaxon) %>% 
    summarise(abun = sum(pollen_abun)) %>% 
    mutate(abun = abun/sum(abun)) %>% 
    # convert to wide format
    pivot_wider(names_from = pollentaxon, values_from = abun) %>% 
    # filter empty columns (taxon not in that country)
    discard(~all(is.na(.))) %>% 
    # replace remaining nas by 0
    replace(is.na(.), 0) %>% 
    # sort alphabetically
    ungroup()
  # save site names with data
  sitenames <- Ab %>% pull(sitename)
  Ab <- Ab %>% 
    dplyr::select(-sitename, sort(colnames(.[-1]))) %>% 
    as.matrix() 
  
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
  plot(results, file = paste0("Convergence_diagnostics/Conv_pollen_", 
                              selectedcountry,"_", selectedtrait, "_", 
                              deparse(substitute(selectedpft)),  ".pdf"))
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_pollen_",
             selectedcountry,"_", selectedtrait, "_", 
             deparse(substitute(selectedpft)),  ".pdf"))
  
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  # summarize results
  res <- data.frame(summary(results),
                    sitename = sitenames, 
                    trait = selectedtrait,
                    country = selectedcountry,
                    growthform = case_when(all(selectedpft == trsh) ~ "trsh",
                                           all(selectedpft == herb) ~ "herb",
                                           TRUE ~ "allpft"),
                    pollination = "allpol")
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", 
                      selectedcountry,"_", selectedtrait, "_", 
                      deparse(substitute(selectedpft)), ".rds"))
}

trsh <- c("tree", "shrub")
herb <- c("grass", "herb")

if (Sys.info()[4] == "LUIN67311") {  
  plan(multisession, workers = 1)
  } else {
  plan(multisession, workers = 6)
  }

furrr::future_map(normaltraits, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpft = trsh))

furrr::future_map(normaltraits[1:3], 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedpft = trsh))

furrr::future_map(normaltraits, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpft = herb))

furrr::future_map(normaltraits[1:3], 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedpft = herb))

}
# Pollination mode ----
if (0) {
cwmlogfunc <- function(selectedtrait, selectedcountry, selectedabun,
                       selectedpolmode){
  trait <- dfTRAIT %>% 
    filter(country == selectedcountry) %>% 
    # drop NA's and 0's
    dplyr::select(pollentaxon, all_of(selectedtrait)) %>% 
    mutate(across(where(is.double), ~na_if(.,0))) %>% 
    drop_na() %>% 
    # Sort alphabetically
    arrange(pollentaxon) %>% 
    # check number of observations
    group_by(pollentaxon) %>% 
    mutate(ntrait = n()) %>% 
    filter(ntrait > 3)
  
  pol <- dfPOL %>% 
    # select country
    filter(country == selectedcountry) %>%
    filter(pollination %in% selectedpolmode) %>%
    dplyr::select(sitename, pollentaxon, 
                  pollen_abun = percent) %>% 
    drop_na()
  
  # save pollen taxa
  pollentaxa <- pol %>% pull(pollentaxon) %>% unique()
  
  ## Model data (capitalized)
  Trait <- trait %>% 
    filter(pollentaxon %in% pollentaxa) %>% 
    pull(selectedtrait) 
  
  Tax <- trait %>% 
    filter(pollentaxon %in% pollentaxa) %>% 
    pull(pollentaxon) %>% 
    as.factor() 
  
  N <- length(Trait)
  
  Ab <- pol %>%
    # this filters taxa with too little trait data
    filter(pollentaxon %in% unique(Tax)) %>% 
    # calculate percentage after removing taxa
    group_by(sitename, pollentaxon) %>% 
    summarise(abun = sum(pollen_abun)) %>% 
    mutate(abun = abun/sum(abun)) %>% 
    # convert to wide format
    pivot_wider(names_from = pollentaxon, values_from = abun) %>% 
    # filter empty columns (taxon not in that country)
    discard(~all(is.na(.))) %>% 
    # replace remaining nas by 0
    replace(is.na(.), 0) %>% 
    # sort alphabetically
    ungroup()
  # save site names with data
  sitenames <- Ab %>% pull(sitename)
  Ab <- Ab %>% 
    dplyr::select(-sitename, sort(colnames(.[-1]))) %>% 
    as.matrix() 
  
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
  plot(results, file = paste0("Convergence_diagnostics/Conv_pollen_",  
                              selectedcountry,"_", selectedtrait, "_", 
                              selectedpolmode,  ".pdf"))
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_pollen_", 
             selectedcountry,"_", selectedtrait, "_", 
             selectedpolmode,  ".pdf"))
  
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  # summarize results
  res <- data.frame(summary(results),
                    sitename = sitenames, 
                    trait = selectedtrait,
                    country = selectedcountry,
                    pollination = selectedpolmode,
                    growthform = "allpft")
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", 
                      selectedcountry,"_", selectedtrait, "_", 
                      selectedpolmode, ".rds"))
}

if (Sys.info()[4] == "LUIN67311") {
plan(multisession, workers = 1)
} else {
  plan(multisession, workers = 6)
}
furrr::future_map(normaltraits,
                  ~cwmlogfunc(selectedtrait = .x,
                              selectedcountry = "Scotland",
                              selectedpolmode = "wind"))

furrr::future_map(normaltraits[1:3],
                  ~cwmlogfunc(selectedtrait = .x,
                              selectedcountry = "Switzerland",
                              selectedpolmode = "wind"))

furrr::future_map(normaltraits,
                  ~cwmlogfunc(selectedtrait = .x,
                              selectedcountry = "Scotland",
                              selectedpolmode = "not wind"))

furrr::future_map(normaltraits[1:3], 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedpolmode = "not wind"))
}

# Pollen correction factors ----
if (1) {
  cwmlogfunc <- function(selectedtrait, selectedcountry, selectedabun){
    trait <- dfTRAIT %>% 
      filter(country == selectedcountry) %>% 
      # drop NA's and 0's
      dplyr::select(pollentaxon, all_of(selectedtrait)) %>% 
      mutate(across(where(is.double), ~na_if(.,0))) %>% 
      drop_na() %>% 
      # Sort alphabetically
      arrange(pollentaxon) %>% 
      # check number of observations
      group_by(pollentaxon) %>% 
      mutate(ntrait = n()) %>% 
      filter(ntrait > 3)
    
    pol <- dfPOL %>% 
      # select country
      filter(country == selectedcountry) %>%
      dplyr::select(sitename, pollentaxon, pollen_abun = all_of(selectedabun)) %>% 
      drop_na()
    
    # save pollen taxa
    pollentaxa <- pol %>% pull(pollentaxon) %>% unique()
    
    ## Model data (capitalized)
    Trait <- trait %>% 
      filter(pollentaxon %in% pollentaxa) %>% 
      pull(selectedtrait) 
    
    Tax <- trait %>% 
      filter(pollentaxon %in% pollentaxa) %>% 
      pull(pollentaxon) %>% 
      as.factor() 
    
    N <- length(Trait)
    
    Ab <- pol %>%
      # this filters taxa with too little trait data
      filter(pollentaxon %in% unique(Tax)) %>% 
      # calculate percentage after removing taxa
      group_by(sitename, pollentaxon) %>% 
      summarise(abun = sum(pollen_abun)) %>% 
      mutate(abun = abun/sum(abun)) %>% 
      # convert to wide format
      pivot_wider(names_from = pollentaxon, values_from = abun) %>% 
      # filter empty columns (taxon not in that country)
      discard(~all(is.na(.))) %>% 
      # replace remaining nas by 0
      replace(is.na(.), 0) %>% 
      # sort alphabetically
      ungroup()
    # save site names with data
    sitenames <- Ab %>% pull(sitename)
    Ab <- Ab %>% 
      dplyr::select(-sitename, sort(colnames(.[-1]))) %>% 
      as.matrix() 
    
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
    plot(results, file = paste0("Convergence_diagnostics/Conv_pollen_", selectedcountry, "_",
                                selectedtrait, "_", selectedabun, ".pdf"))
    # gelman plots
    pdf(paste0("Convergence_diagnostics/Gelman_pollen_", selectedcountry,
               "_", selectedtrait,"_", selectedabun, ".pdf"))
    
    gelman.plot(results, ask = FALSE)
    dev.off()
    
    # summarize results
    res <- data.frame(summary(results),
                      sitename = sitenames,
                      trait = selectedtrait,
                      country = selectedcountry,
                      correction = selectedabun)
    saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", selectedcountry,"_",
                        selectedtrait, "_", selectedabun, ".rds"))
  }
  
  combo <- expand_grid(trait = normaltraits, pollen_abun = str_subset(colnames(dfPOL), "percent"))
  
  if (Sys.info()[4] == "LUIN67311") {
    plan(multisession, workers = 1)
  } else {
    plan(multisession, workers = 6)
  }
  # furrr::future_map2(combo$trait, combo$pollen_abun,
  #                   ~cwmlogfunc(selectedtrait = .x,
  #                               selectedcountry = "Scotland",
  #                               selectedabun = .y))
  
  combo <- expand_grid(trait = normaltraits[1:3], 
                       pollen_abun = str_subset(colnames(dfPOL), "percent"))
  furrr::future_map2(combo$trait, combo$pollen_abun,
                     ~cwmlogfunc(selectedtrait = .x, 
                                 selectedcountry = "Switzerland",
                                 selectedabun = .y))
}