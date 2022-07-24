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
if(!require(LCVP)) install.packages("LCVP")
if(!require(lcvplants)) install.packages("lcvplants")

## Load and prepare data -----
# Trait data ----
dfTRAIT <- readRDS("RDS_files/02_Gapfilled_traits.rds")

# Pollen data ----
dfPOL_Scot <- readRDS("RDS_files/01_pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_pollen_data_Swiss.rds")
dfPOL <- bind_rows("Scotland" = dfPOL_Scot, "Switzerland" =
                     dfPOL_Swiss, .id = "country") %>% 
  ungroup()

selectedcountry <- "Scotland"
selectedtrait <- "LA"
selectedabun <- str_subset(colnames(dfPOL), "percent")[57]
# Calculate CWM ----
cwmlogfunc <- function(selectedtrait, selectedcountry, selectedabun){
  pol <- dfPOL %>% 
    # select country
    filter(country == selectedcountry) %>%
    dplyr::select(sitename, pollentaxon, pollen_abun = all_of(selectedabun)) %>% 
    # filter out unidentified
    filter(!pollentaxon == "unindentified",
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
    arrange(match(pollentaxon, pollentaxa)) %>% # order in the same way as pollen data
    pull(selectedtrait) 
  
  Tax <- dfTRAIT %>%
    filter(country == selectedcountry) %>% 
    filter(pollentaxon %in% pollentaxa) %>%
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
  Mean <- mean(Trait)
  SD <- sd(Trait)
  
  cat(
    "model {
      for (i in 1:Nsites){
    cwm[i] ~ dnorm(mean.tax[r[i]], tau.tax[r[i]])
    r[i] ~ dcat(Ab[i, 1:Ntax])
    }
    for (i in 1:N) {
      Trait[i] ~ dlnorm(mean.tax[Tax[i]], tau.tax[Tax[i]])
    }
    for (l in 1:Ntax){
    # vague priors taxon level
      mean.tax[l] ~ dnorm(0, 10^-4)
      tau.tax[l] ~ dgamma(0.001, 0.001)
      sd.tax[l] <- sqrt(1/tau.tax[l])
    }
    #data# N, Nsites, Trait, Tax, Ntax, Ab
    #monitor# mean.tax, sd.tax, cwm
    #inits# mean.tax, tau.tax
  }", file = "CWM_log.txt"
  )
  
  # initial values
  mean.tax <- list(chain1 = rep(min(Trait), Ntax), chain2 = rep(max(Trait),Ntax))
  tau.tax <- list(chain1 = rep(1/(SD*0.1)^2, Ntax), chain2 = rep(1/(SD*10)^2,
                                                                 Ntax))
  
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
                    sitename = sitenames)
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", selectedcountry,"_",
                      selectedtrait, "_", selectedabun, ".rds"))
}

normaltraits <- c("PlantHeight","SLA","LA")
combo <- expand_grid(trait = normaltraits, pollen_abun = str_subset(colnames(dfPOL), "percent"))

plan(multisession, workers = 4)
furrr::future_map2(combo$trait, combo$pollen_abun,
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedabun = .y))

plan(multisession, workers = 4)
furrr::future_map2(combo$trait, combo$pollen_abun,
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedabun = .y))

# Growth form ----
cwmlogfunc <- function(selectedtrait, selectedcountry, selectedabun,
                       selectedpft){
  pol <- dfPOL %>% 
    # select country
    filter(country == selectedcountry) %>%
    # filter out specific pollinationmode
    filter(growthform %in% selectedpft) %>%
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
    arrange(match(pollentaxon, pollentaxa)) %>% # order in the same way as pollen data
    pull(selectedtrait) 
  
  Tax <- dfTRAIT %>%
    filter(country == selectedcountry) %>% 
    filter(pollentaxon %in% pollentaxa) %>%
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
  Mean <- mean(Trait)
  SD <- sd(Trait)
  
  # initial values
  mean.tax <- list(chain1 = rep(min(Trait), Ntax), chain2 = rep(max(Trait),Ntax))
  tau.tax <- list(chain1 = rep(1/(SD*0.1)^2, Ntax), chain2 = rep(1/(SD*10)^2,
                                                                 Ntax))
  
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
                   sitename = sitenames)
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", 
                      selectedcountry,"_", selectedtrait, "_", 
                      deparse(substitute(selectedpft)), ".rds"))
}

trsh <- c("tree", "shrub")
herb <- c("grass", "herb")

plan(multisession, workers = 4)
furrr::future_map(combo$trait, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpft = trsh))

plan(multisession, workers = 4)
furrr::future_map(combo$trait, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedpft = trsh))

plan(multisession, workers = 4)
furrr::future_map(combo$trait, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpft = herb))

plan(multisession, workers = 4)
furrr::future_map(combo$trait, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedpft = herb))


# Pollination mode ----
cwmlogfunc <- function(selectedtrait, selectedcountry, selectedabun,
                       selectedpolmode){
  pol <- dfPOL %>% 
    # select country
    filter(country == selectedcountry) %>%
    # filter out specific pollinationmode
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
    arrange(match(pollentaxon, pollentaxa)) %>% # order in the same way as pollen data
    pull(selectedtrait) 
  
  Tax <- dfTRAIT %>%
    filter(country == selectedcountry) %>% 
    filter(pollentaxon %in% pollentaxa) %>%
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
  Mean <- mean(Trait)
  SD <- sd(Trait)
  
  # initial values
  mean.tax <- list(chain1 = rep(min(Trait), Ntax), chain2 = rep(max(Trait),Ntax))
  tau.tax <- list(chain1 = rep(1/(SD*0.1)^2, Ntax), chain2 = rep(1/(SD*10)^2,
                                                                 Ntax))
  
  # run model
  results <- run.jags("CWM_log.txt", n.chains = 2, monitor = "cwm")
  
  # trace plots
  plot(results, file = paste0("Convergence_diagnostics/Conv_pollen_",  
                              selectedcountry,"_", selectedtrait, "_", 
                              deparse(substitute(selectedpolmode)),  ".pdf"))
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_pollen_", 
             selectedcountry,"_", selectedtrait, "_", 
             deparse(substitute(selectedpolmode)),  ".pdf"))
  
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  # summarize results
  res <- data.frame(summary(results),
                    sitename = sitenames)
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", 
                      selectedcountry,"_", selectedtrait, "_", 
                      deparse(substitute(selectedpolmode)), ".rds"))
}

plan(multisession, workers = 4)
furrr::future_map(combo$trait, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpolmode = "wind"))

plan(multisession, workers = 4)
furrr::future_map(combo$trait, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedpolmode = "wind"))

plan(multisession, workers = 4)
furrr::future_map(combo$trait, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Scotland",
                              selectedpolmode = "no wind"))

plan(multisession, workers = 4)
furrr::future_map(combo$trait, 
                  ~cwmlogfunc(selectedtrait = .x, 
                              selectedcountry = "Switzerland",
                              selectedpolmode = "no wind"))