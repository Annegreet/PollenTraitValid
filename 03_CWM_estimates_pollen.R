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
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(rjags))  install.packages("rjags")
if(!require(runjags)) install.packages("runjags")
if(!require(mcmc)) install.packages("mcmc")
if(!require(coda)) install.packages("coda")
if(!require(Hmisc)  )install.packages("Hmisc") 
if(!require(furrr))  install.packages("furrr")

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

# selectedcountry <- "Switzerland"
# selectedtrait <- "SLA"
# selectedabun <- str_subset(colnames(dfPOL), "percent")[52]

normaltraits <- c("PlantHeight","SLA","LA", "LeafN", "LDMC")

# Growth form ----
if(1){
cwmlogfunc <- function(selectedtrait, selectedcountry, selectedabun,
                       selectedpft){
  pol <- dfPOL %>% 
    # select country
    filter(country == selectedcountry) %>%
    # filter out specific pollination mode
    filter(growthform %in% selectedpft) %>%
    dplyr::select(sitename, pollentaxon, 
                  abun = adjustedpercent_mean) %>% 
    # filter out unidentified
    filter(!pollentaxon == "unindentified",
           !is.na(abun)) %>% 
    # recalculate pollen abundance on the plot level
    group_by(sitename, pollentaxon) %>% 
    summarise(abun = sum(abun)) %>% 
    mutate(abun = abun/sum(abun)) %>% 
    filter(!is.nan(abun)) %>%  # filters sites without species of particular selection
    # convert to wide format
    pivot_wider(names_from = pollentaxon, values_from = abun,
                values_fill = 0) %>% 
    ungroup() %>% 
    discard(~all(is.na(.))) 
  
  # save sitenames
  sitenames <- pol %>% 
    pull(sitename)
  
  # save pollen taxa
  pollentaxa <- colnames(pol[, 2:ncol(pol)])
  
  Trait <- dfTRAIT %>% 
    filter(country == selectedcountry) %>% 
    filter(pollentaxon %in% pollentaxa) %>%
    # drop NA's and 0's
    dplyr::select(pollentaxon, all_of(selectedtrait)) %>% 
    mutate(across(where(is.double), ~na_if(.,0))) %>% 
    drop_na() %>% 
    arrange(match(pollentaxon, pollentaxa)) %>% # order in the same way as pollen data
    pull(selectedtrait) 
  
  Tax <- dfTRAIT %>%
    filter(country == selectedcountry) %>% 
    filter(pollentaxon %in% pollentaxa) %>%
    # drop NA's and 0's
    dplyr::select(pollentaxon, all_of(selectedtrait)) %>% 
    mutate(across(where(is.double), ~na_if(.,0))) %>% 
    drop_na() %>% 
    arrange(match(pollentaxon, pollentaxa)) %>% # order in the same way as pollen data
    pull(pollentaxon) %>% 
    as.factor() 
  
  N <- length(Trait)
  
  Ab <- pol %>%
    dplyr::select(-sitename) %>%
    as.matrix() 
    
  Ntax <- ncol(Ab)
  Nsites <- nrow(Ab) 
  
  # for setting priors taxon level
  Mean <- aggregate(Trait, list(Tax), FUN=mean)[,"x"] 
  SD <- aggregate(Trait, list(Tax), FUN=sd)[,"x"] %>% 
    replace_na(., mean(., na.rm=TRUE)) # in case of only 1 obs set sd to mean sd
  
  # initial values
  mean.tax <- list(chain1 = rep(min(Trait), Ntax), chain2 = rep(max(Trait),Ntax))
  tau.tax <- list(chain1 = 1/SD^2*0.001, chain2 = 1/SD^2*1000)
  
  # run model
  results <- run.jags("CWM_log.txt", n.chains = 2, monitor = "cwm")
  
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

plan(multisession, workers = 6)
furrr::future_map(normaltraits, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpft = trsh))

plan(multisession, workers = 6)
furrr::future_map(normaltraits[1:3], 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedpft = trsh))

plan(multisession, workers = 6)
furrr::future_map(normaltraits, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpft = herb))

plan(multisession, workers = 6)
furrr::future_map(normaltraits[1:3], 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedpft = herb))

}
# Pollination mode ----
if(1){
cwmlogfunc <- function(selectedtrait, selectedcountry, selectedabun,
                       selectedpolmode){
  pol <- dfPOL %>% 
    # select country
    filter(country == selectedcountry) %>%
    # filter out specific pollination mode
    filter(pollination %in% selectedpolmode) %>%
    dplyr::select(sitename, pollentaxon, abun = adjustedpercent_mean) %>% 
    # filter out unidentified
    filter(!pollentaxon == "unindentified",
           !is.na(abun)) %>% 
    # recalculate pollen abundance on the plot level
    group_by(sitename, pollentaxon) %>% 
    summarise(abun = sum(abun)) %>% 
    mutate(abun = abun/sum(abun)) %>% 
    filter(!is.nan(abun)) %>%  # filters sites without species of particular selection
    # convert to wide format
    pivot_wider(names_from = pollentaxon, values_from = abun,
                values_fill = 0) %>% 
    ungroup() %>% 
    discard(~all(is.na(.))) 
  
  # save sitenames
  sitenames <- pol %>% 
    pull(sitename)
  # save pollen taxa
  pollentaxa <- colnames(pol[, 2:ncol(pol)])
  
  Trait <- dfTRAIT %>% 
    filter(country == selectedcountry) %>% 
    filter(pollentaxon %in% pollentaxa) %>%
    # drop NA's and 0's
    dplyr::select(pollentaxon, all_of(selectedtrait)) %>% 
    mutate(across(where(is.double), ~na_if(.,0))) %>% 
    drop_na() %>% 
    arrange(match(pollentaxon, pollentaxa)) %>% # order in the same way as pollen data
    pull(selectedtrait) 
  
  Tax <- dfTRAIT %>%
    filter(country == selectedcountry) %>% 
    filter(pollentaxon %in% pollentaxa) %>%
    # drop NA's and 0's
    dplyr::select(pollentaxon, all_of(selectedtrait)) %>% 
    mutate(across(where(is.double), ~na_if(.,0))) %>% 
    drop_na() %>% 
    arrange(match(pollentaxon, pollentaxa)) %>% # order in the same way as pollen data
    pull(pollentaxon) %>% 
    as.factor() 
  
  N <- length(Trait)
  
  Ab <- pol %>%
    dplyr::select(-sitename) %>%
    as.matrix() 
  
  Ntax <- ncol(Ab)
  Nsites <- nrow(Ab) 
  
  # for setting priors taxon level
  Mean <- aggregate(Trait, list(Tax), FUN=mean)[,"x"] 
  SD <- aggregate(Trait, list(Tax), FUN=sd)[,"x"] %>% 
    replace_na(., mean(., na.rm=TRUE)) # in case of only 1 obs set sd to mean sd
  
  # initial values
  mean.tax <- list(chain1 = rep(min(Trait), Ntax), chain2 = rep(max(Trait),Ntax))
  tau.tax <- list(chain1 = 1/SD^2*0.001, chain2 = 1/SD^2*1000)
  
  # run model
  results <- run.jags("CWM_log.txt", n.chains = 2, monitor = "cwm")
  
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

plan(multisession, workers = 6)
furrr::future_map(normaltraits,
                  ~cwmlogfunc(selectedtrait = .x,
                              selectedcountry = "Scotland",
                              selectedpolmode = "wind"))

plan(multisession, workers = 6)
furrr::future_map(normaltraits[1:3],
                  ~cwmlogfunc(selectedtrait = .x,
                              selectedcountry = "Switzerland",
                              selectedpolmode = "wind"))

plan(multisession, workers = 6)
furrr::future_map(normaltraits, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpolmode = "not wind"))

plan(multisession, workers = 6)
furrr::future_map(normaltraits[1:3], 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedpolmode = "not wind"))
}

# Pollen correction factors ----
if(1){
  cwmlogfunc <- function(selectedtrait, selectedcountry, selectedabun){
    pol <- dfPOL %>% 
      # select country
      filter(country == selectedcountry) %>%
      dplyr::select(sitename, pollentaxon, pollen_abun = all_of(selectedabun)) %>% 
      # filter out unidentified
      filter(!pollentaxon %in% c("unindentified", "Liliaceae", "Pedicularis palustris"),
             !is.na(pollen_abun)) %>% 
      pivot_wider(names_from = pollentaxon, values_from = pollen_abun) %>% 
      # filter empty columns (taxon not in that country)
      discard(~all(is.na(.))) 
    
    # save sitenames
    sitenames <- pol %>% 
      pull(sitename)
    # save pollen taxa
    pollentaxa <- colnames(pol[, 2:ncol(pol)])
    
    Trait <- dfTRAIT %>% 
      filter(country == selectedcountry) %>% 
      filter(pollentaxon %in% pollentaxa) %>%
      # drop NA's and 0's
      dplyr::select(pollentaxon, all_of(selectedtrait)) %>% 
      mutate(across(where(is.double), ~na_if(.,0))) %>% 
      drop_na() %>% 
      arrange(match(pollentaxon, pollentaxa)) %>% # order in the same way as pollen data
      pull(selectedtrait) 
    
    Tax <- dfTRAIT %>%
      filter(country == selectedcountry) %>% 
      filter(pollentaxon %in% pollentaxa) %>%
      # drop NA's and 0's
      dplyr::select(pollentaxon, all_of(selectedtrait)) %>% 
      mutate(across(where(is.double), ~na_if(.,0))) %>% 
      drop_na() %>% 
      pull(pollentaxon) %>% 
      as.factor() 
    
    N <- length(Trait)
    
    Ab <- pol %>%
      dplyr::select(-sitename) %>%
      as.matrix() 
    
    Ntax <- ncol(Ab)
    Nsites <- nrow(Ab) 
    
    # for setting priors taxon level
    Mean <- aggregate(Trait, list(Tax), FUN=mean)[,"x"] 
    SD <- aggregate(Trait, list(Tax), FUN=sd)[,"x"] %>% 
      replace_na(., mean(., na.rm=TRUE)) # in case of only 1 obs set sd to mean sd
    
    # initial values
    mean.tax <- list(chain1 = rep(min(Trait), Ntax), chain2 = rep(max(Trait),Ntax))
    tau.tax <- list(chain1 = 1/SD^2*0.001, chain2 = 1/SD^2*1000)
    
    # run model
    results <- run.jags("CWM_log.txt", n.chains = 2, monitor = "cwm")
    
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
  
  plan(multisession, workers = 6)
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