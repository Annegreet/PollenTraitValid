# CWM vegetation Scotland ----
cwm_veg_scot <- function(selectedabun, selectedtrait, selectedtaxres, selectedpolmode, 
                         selectedpft, notinpollen = FALSE, taxtrait = FALSE){
  ## prepare abundance data
  abun <- lABUN %>% 
    # select zone 
    pluck(selectedabun) %>% 
    # optional filter for polmode and pft
    filter(polmode %in% selectedpolmode) %>% 
    filter(growthform %in% selectedpft) %>% 
    # optional to filter out taxa that are not in pollen data
    {if (notinpollen) filter(., !pollentaxon == "Not in pollen data") else . } %>% 
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
  
  # filter species with too little observations (10 or less)
  nobs <- lTrait %>%
    pluck(selectedtrait) %>%
    arrange(selectedtaxres) %>% 
    group_by(across(all_of(selectedtaxres))) %>%
    summarise(n = n()) %>%
    filter(n <= 10) %>%
    pull(selectedtaxres)
  
  # find taxa with not or too little data in TRY gapfilled data
  TRY <- dfTRY %>% 
    dplyr::select(tax = all_of(selectedtaxres), 
                  all_of(selectedtrait)) %>% 
    filter(tax %in% c(missingtax, nobs)) 
  
  # Create trait data df with field an try data
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
    filter(ntrait >= 10)
  
  ## Model data (capitalized)
  Trait <- trait %>% 
    pull(selectedtrait)
  
  Tax <- trait  %>% 
    pull(tax) %>% 
    as.factor()
  
  # Save missing species that are not in TRY neither
  taxa[!taxa %in% unique(Tax)] %>% saveRDS(file = paste0("RDS_files/Missing_taxa_", selectedabun,"_Scotland_", selectedtrait,
                                                         "_", selectedtaxres, "_", paste(selectedpolmode, collapse = ""),"_",
                                                         paste(selectedpft, collapse = ""),ifelse(notinpollen == T, "_notinpollen",""), ".rds")
  )
  
  # Prepare abundance matrix
  Ab <- abun %>%
    # filter species without enough data 
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
  
  # for setting inits taxon level
  MeanLog <- aggregate(Trait, list(Tax), FUN = mean)[,"x"] %>% log()
  SDLog <- aggregate(Trait, list(Tax), FUN = sd)[,"x"] %>%
    log() %>% 
    na_if(., 0) # prior cant take 0 or negative
  SDLog[SDLog < 0] <- NA  
  SDLog <- replace_na(SDLog, mean(SDLog, na.rm = TRUE))  # in case of only 1 obs set sd to mean sd
  
  # initial values
  mean.tax.log <- list(chain1 = MeanLog*0.01, chain2 = MeanLog*100)
  sd.tax.log <- list(chain1 = SDLog*0.01, chain2 = SDLog*100)
  cwm <- list(chain1 = rep(min(log(Trait)), Nsites), 
              chain2 = rep(max(log(Trait)), Nsites))
  
  # run model
  results <- run.jags("CWM_log_prior_data.txt", n.chains = 2)
  
  # trace plots
  plot(results, file =
         paste0("Convergence_diagnostics/Conv_", selectedabun,"_Scotland_", selectedtrait,
                "_", selectedtaxres, "_", paste(selectedpolmode, collapse = ""),"_",
                paste(selectedpft, collapse = ""),ifelse(notinpollen == T, "_notinpollen",""), ".pdf")
  )
  
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_", selectedabun,"_Scotland_", selectedtrait,
             "_", selectedtaxres, "_", paste(selectedpolmode, collapse = ""),"_",
             paste(selectedpft, collapse = ""),ifelse(notinpollen == T, "_notinpollen",""), ".pdf")
  )
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  # save CWM values
  res <- as.data.frame(summary(results, vars = "cwm")) %>%
    # filter cwm values
    rownames_to_column("parameter") %>%
    filter(!str_detect(parameter, "zcwm")) %>% 
    mutate(sitename = sitenames,
           trait = selectedtrait,
           growthform = case_when(all(selectedpft %in% trsh) ~ "trsh",
                                  all(selectedpft %in% herb) ~ "herb",
                                  TRUE ~ "allpft"),
           pollination = case_when(all(selectedpolmode %in% wind) ~ "wind",
                                   all(selectedpolmode %in% nowind) ~ "no wind",
                                   TRUE ~ "allpol"),
           taxres = selectedtaxres,
           zone = selectedabun,
           notinpollendata = notinpollen)
  
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_",selectedabun,"_Scotland_", selectedtrait,
                      "_", selectedtaxres, "_",
                      paste(selectedpolmode, collapse = ""),"_", 
                      paste(selectedpft, collapse = ""), ifelse(notinpollen == T, "_notinpollen",""),".rds"))
  
  # save standardized CWM values
  zres <- as.data.frame(summary(results, vars = "zcwm")) %>%
    # filter cwm values
    rownames_to_column("parameter") %>%
    mutate(sitename = sitenames,
           trait = selectedtrait,
           growthform = case_when(all(selectedpft %in% trsh) ~ "trsh",
                                  all(selectedpft %in% herb) ~ "herb",
                                  TRUE ~ "allpft"),
           pollination = case_when(all(selectedpolmode %in% wind) ~ "wind",
                                   all(selectedpolmode %in% nowind) ~ "no wind",
                                   TRUE ~ "allpol"),
           taxres = selectedtaxres,
           zone = selectedabun,
           notinpollendata = notinpollen)
  
  saveRDS(zres, paste0("RDS_files/03_zCWM_estimates_",selectedabun,"_Scotland_", selectedtrait,
                      "_", selectedtaxres, "_",
                      paste(selectedpolmode, collapse = ""),"_", 
                      paste(selectedpft,collapse = ""), ifelse(notinpollen == T, "_notinpollen",""),".rds"))
  
  # Save taxon level trait values
  if (taxtrait) {
    res_tax <- as.data.frame(summary(results, vars = c("mean.tax.log", "sd.tax.log"))) %>%
      # filter taxon trait estimates
      rownames_to_column("parameter") %>%
      dplyr::select(parameter, Mean) %>%
      mutate(parameter = str_extract(parameter, pattern = "mean.tax|sd.tax"),
             tax = c(unique(Tax), unique(Tax))) %>% 
      pivot_wider(names_from = parameter, values_from = Mean) 
    saveRDS(res_tax, paste0("RDS_files/03_taxon_trait_estimates_", selectedabun,
                            "_Scotland_", selectedtrait, ".rds"))
    
  }
  
}

# CWM vegetation Switzerland ----
cwm_veg_swiz <- function(selectedabun, selectedtrait, selectedtaxres, selectedpolmode, selectedpft,
                         notinpollen = FALSE, taxtrait = FALSE){
  ## Select abundance data
  abun <- lABUN %>% 
    pluck(selectedabun) %>% 
    # optional filter for polmode and pft
    filter(polmode %in% selectedpolmode) %>% 
    filter(growthform %in% selectedpft) %>% 
    {if (notinpollen) filter(., !pollentaxon == "Not in pollen data") else . } %>% 
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
  
  # find species with too little observations (10 or less)
  nobs <- bdm_traits %>%
    arrange(selectedtaxres) %>%
    group_by(across(all_of(selectedtaxres))) %>%
    summarise(n = n()) %>%
    filter(n <= 10) %>%
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
    filter(ntrait >= 10)
  
   ## Model data (capitalized)
  Trait <- trait %>% 
    pull(selectedtrait)
  
  Tax <- trait  %>% 
    pull(tax) %>% 
    as.factor()
  # Save missing species that are not in TRY neither
  taxa[!taxa %in% unique(Tax)] %>% saveRDS(file = paste0("RDS_files/Missing_taxa_", selectedabun,"_Switzerland_", selectedtrait,
                                                         "_", selectedtaxres, "_", paste(selectedpolmode, collapse = ""),"_",
                                                         paste(selectedpft, collapse = ""), ifelse(notinpollen == T, "_notinpollen",""), ".rds"))
  
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
  
  # for setting inits taxon level
  MeanLog <- aggregate(Trait, list(Tax), FUN = mean)[,"x"] %>% log()
  SDLog <- aggregate(Trait, list(Tax), FUN = sd)[,"x"] %>%
    log() %>% 
    na_if(., 0) # prior cant take 0 or negative
  SDLog[SDLog < 0] <- NA  
  SDLog <- replace_na(SDLog, mean(SDLog, na.rm = TRUE))  # in case of only 1 obs set sd to mean sd
  
  # initial values
  mean.tax.log <- list(chain1 = MeanLog*0.01, chain2 = MeanLog*100)
  sd.tax.log <- list(chain1 = SDLog*0.01, chain2 = SDLog*100)
  cwm <- list(chain1 = rep(min(log(Trait)), Nsites), 
              chain2 = rep(max(log(Trait)), Nsites))
  
  # run model
  results <- run.jags("CWM_log_prior_data.txt", n.chains = 2)
  
  # trace plots
  plot(results, file =
         paste0("Convergence_diagnostics/Conv_", selectedabun,"_Switzerland_", selectedtrait,
                "_", selectedtaxres, "_", paste(selectedpolmode, collapse = ""),"_",
                paste(selectedpft, collapse = ""),ifelse(notinpollen == T, "_notinpollen",""), ".pdf")
  )
  
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_", selectedabun,"_Switzerland_", selectedtrait,
             "_", selectedtaxres, "_", paste(selectedpolmode, collapse = ""),"_",
             paste(selectedpft, collapse = ""),ifelse(notinpollen == T, "_notinpollen",""), ".pdf")
  )
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  # save CWM values
  res <- as.data.frame(summary(results, vars = "cwm")) %>% 
    # filter cwm values
    rownames_to_column("parameter") %>%
    filter(!str_detect(parameter, "zcwm")) %>% 
    mutate(sitename = sitenames,
           trait = selectedtrait,
           growthform = case_when(all(selectedpft %in% trsh) ~ "trsh",
                                  all(selectedpft %in% herb) ~ "herb",
                                  TRUE ~ "allpft"),
           pollination = case_when(all(selectedpolmode %in% wind) ~ "wind",
                                   all(selectedpolmode %in% nowind) ~ "no wind",
                                   TRUE ~ "allpol"),
           taxres = selectedtaxres,
           zone = selectedabun)
  
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_Switzerland_", selectedtrait,
                      "_",selectedabun,"_", selectedtaxres, "_", 
                      paste(selectedpolmode,collapse = ""),"_", 
                      paste(selectedpft,collapse = ""), ifelse(notinpollen == T, "_notinpollen",""),".rds"))
  
  # save standardized CWM values
  zres <- as.data.frame(summary(results, vars = "zcwm")) %>% 
    mutate(sitename = sitenames,
           trait = selectedtrait,
           growthform = case_when(all(selectedpft %in% trsh) ~ "trsh",
                                  all(selectedpft %in% herb) ~ "herb",
                                  TRUE ~ "allpft"),
           pollination = case_when(all(selectedpolmode %in% wind) ~ "wind",
                                   all(selectedpolmode %in% nowind) ~ "no wind",
                                   TRUE ~ "allpol"),
           taxres = selectedtaxres,
           zone = selectedabun)
  
  saveRDS(zres, paste0("RDS_files/03_zCWM_estimates_Switzerland_", selectedtrait,
                      "_",selectedabun,"_", selectedtaxres, "_", 
                      paste(selectedpolmode,collapse = ""),"_", 
                      paste(selectedpft,collapse = ""), ifelse(notinpollen == T, "_notinpollen",""),".rds"))
  
  # Save taxon level trait values
  if (taxtrait) {
    res_tax <- as.data.frame(summary(results, vars = c("mean.tax.log", "sd.tax.log"))) %>%
      # filter taxon trait estimates
      rownames_to_column("parameter") %>%
      dplyr::select(parameter, Mean) %>%
      mutate(parameter = str_extract(parameter, pattern = "mean.tax|sd.tax"),
             tax = c(unique(Tax), unique(Tax))) %>% 
      pivot_wider(names_from = parameter, values_from = Mean) 
    saveRDS(res_tax, paste0("RDS_files/03_taxon_trait_estimates_", selectedabun,
                            "_Switzerland_", selectedtrait, ".rds"))
  }
  
}

# CWM remove species ----
cwm_veg_remove <- function(selectedabun, selectedtrait, removedspec){
  abun <- lABUN %>% 
    pluck(selectedabun) %>% 
    filter(!stand.spec %in% removedspec) %>% 
    # rename column to selected taxonomic resolution
    rename(tax = stand.spec) %>% 
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
    filter(n <= 10) %>%
    pull(selectedtaxres)
  
  # find in gapdata
  TRY <- dfTRY %>% 
    dplyr::select(tax = all_of(selectedtaxres), 
                  all_of(selectedtrait)) %>% 
    filter(tax %in% c(missingtax, nobs)) 
  
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
    filter(ntrait >= 10)
  
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
  
  # for setting inits taxon level
  MeanLog <- aggregate(Trait, list(Tax), FUN = mean)[,"x"] %>% log()
  SDLog <- aggregate(Trait, list(Tax), FUN = sd)[,"x"] %>%
    log() %>% 
    na_if(., 0) # prior cant take 0 or negative
  SDLog[SDLog < 0] <- NA  
  SDLog <- replace_na(SDLog, mean(SDLog, na.rm = TRUE))  # in case of only 1 obs set sd to mean sd
  
  # initial values
  mean.tax.log <- list(chain1 = MeanLog*0.01, chain2 = MeanLog*100)
  sd.tax.log <- list(chain1 = SDLog*0.01, chain2 = SDLog*100)
  cwm <- list(chain1 = rep(min(log(Trait)), Nsites), 
              chain2 = rep(max(log(Trait)), Nsites))
  
  # run model
  results <- run.jags("CWM_log_prior_data.txt", n.chains = 2)
  
  
  res <- as.data.frame(summary(results, vars = "cwm")) %>%
    # filter cwm values
    rownames_to_column("parameter") %>%
    filter(!str_detect(parameter, "zcwm")) %>% 
    mutate(sitename = sitenames,
           trait = selectedtrait,
           zone = selectedabun,
           removed_spec = removedspec)
  
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_",selectedabun,"_Scotland_", selectedtrait,"_",
                      removedspec,".rds"))
}

# CWM pollen no subset ----
cwm_pol <- function(selectedtrait, selectedcountry, selectedabun, taxtrait = FALSE){
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
    filter(ntrait >= 10)
  
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
  mean.tax.log <- list(chain1 = MeanLog*0.01, chain2 = MeanLog*100)
  sd.tax.log <- list(chain1 = SDLog*0.01, chain2 = SDLog*100)
  cwm <- list(chain1 = rep(min(log(Trait)), Nsites), 
              chain2 = rep(max(log(Trait)), Nsites))
  # run model
  results <- run.jags("CWM_log_prior_data.txt", n.chains = 2)
  
  # trace plots
  plot(results, 
       file = paste0("Convergence_diagnostics/Conv_pollen_", selectedcountry, "_",
                              selectedtrait, "_", selectedabun, ".pdf"))
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_pollen_", selectedcountry,
             "_", selectedtrait,"_", selectedabun, ".pdf"))
  
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  # Save cwm values
  res <-  res <- as.data.frame(summary(results, vars = "cwm")) %>%
    # filter cwm values
    rownames_to_column("parameter") %>%
    filter(!str_detect(parameter, "zcwm")) %>% 
    mutate(sitename = sitenames,
           trait = selectedtrait,
           country = selectedcountry,
           correction = selectedabun)
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", selectedcountry,"_",
                      selectedtrait, "_", selectedabun, ".rds"))
  
  # Save standardized cwm values
  zres <- data.frame(summary(results, vars = "zcwm"),
                    sitename = sitenames,
                    trait = selectedtrait,
                    country = selectedcountry,
                    correction = selectedabun)
  saveRDS(zres, paste0("RDS_files/03_zCWM_estimates_pollen_", selectedcountry,"_",
                      selectedtrait, "_", selectedabun, ".rds"))
  
  # Save taxon level trait values
  if (taxtrait) {
    res_tax <- as.data.frame(summary(results, vars = c("mean.tax.log", "sd.tax.log"))) %>%
      # filter taxon trait estimates
      rownames_to_column("parameter") %>%
      dplyr::select(parameter, Mean) %>%
      mutate(parameter = str_extract(parameter, pattern = "mean.tax|sd.tax"),
             tax = rep(as.character(unique(Tax)), 2)) %>% 
      pivot_wider(names_from = parameter, values_from = Mean) 
    
    saveRDS(res_tax, paste0("RDS_files/03_taxon_trait_estimates_",selectedcountry,"_",
                            selectedabun, "_", selectedtrait, ".rds"))
  }
}

# CWM pollen plant functional type ----
cwm_pol_pft <- function(selectedtrait, selectedcountry, selectedabun,
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
    filter(ntrait >= 10)
  
  pol <- dfPOL %>% 
    # select country
    filter(country == selectedcountry) %>%
    filter(growthform %in% selectedpft) %>%
    dplyr::select(sitename, pollentaxon, 
                  pollen_abun = all_of(selectedabun)) %>% 
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
  mean.tax.log <- list(chain1 = MeanLog*0.01, chain2 = MeanLog*100)
  sd.tax.log <- list(chain1 = SDLog*0.01, chain2 = SDLog*100)
  cwm <- list(chain1 = rep(min(log(Trait)), Nsites), 
              chain2 = rep(max(log(Trait)), Nsites))
  # run model
  results <- run.jags("CWM_log_prior_data.txt", n.chains = 2)
  
  # trace plots
  plot(results, file = paste0("Convergence_diagnostics/Conv_pollen_", 
                              selectedcountry,"_", selectedabun, "_", selectedtrait, "_", 
                              deparse(substitute(selectedpft)),  ".pdf"))
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_pollen_",
             selectedcountry,"_", selectedabun,"_", selectedtrait, "_", 
             deparse(substitute(selectedpft)),  ".pdf"))
  
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  # save cwm values
  res <- as.data.frame(summary(results, vars = "cwm")) %>%
    # filter cwm values
    rownames_to_column("parameter") %>%
    filter(!str_detect(parameter, "zcwm")) %>% 
    mutate( sitename = sitenames, 
           trait = selectedtrait,
           country = selectedcountry,
           correction = selectedabun,
           growthform = case_when(all(selectedpft %in% trsh) ~ "trsh",
                                  all(selectedpft %in% herb) ~ "herb",
                                  TRUE ~ "allpft"),
           pollination = "allpol")
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", 
                      selectedcountry,"_",selectedabun,"_", selectedtrait, "_", 
                      deparse(substitute(selectedpft)), ".rds"))
  # save standardized cwm values
  zres <- data.frame(summary(results, vars = "zcwm"),
                    sitename = sitenames, 
                    trait = selectedtrait,
                    country = selectedcountry,
                    correction = selectedabun,
                    growthform = case_when(all(selectedpft %in% trsh) ~ "trsh",
                                           all(selectedpft %in% herb) ~ "herb",
                                           TRUE ~ "allpft"),
                    pollination = "allpol")
  saveRDS(zres, paste0("RDS_files/03_zCWM_estimates_pollen_", 
                      selectedcountry,"_",selectedabun,"_", selectedtrait, "_", 
                      deparse(substitute(selectedpft)), ".rds"))
  
}

# CWM pollen pollination ----
cwm_pol_pol <- function(selectedtrait, selectedcountry, selectedabun,
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
    filter(ntrait >= 10)
  
  pol <- dfPOL %>% 
    # select country
    filter(country == selectedcountry) %>%
    filter(pollination %in% selectedpolmode) %>%
    dplyr::select(sitename, pollentaxon, 
                  pollen_abun = all_of(selectedabun)) %>% 
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
  mean.tax.log <- list(chain1 = MeanLog*0.01, chain2 = MeanLog*100)
  sd.tax.log <- list(chain1 = SDLog*0.01, chain2 = SDLog*100)
  cwm <- list(chain1 = rep(min(log(Trait)), Nsites), 
              chain2 = rep(max(log(Trait)), Nsites))
  # run model
  results <- run.jags("CWM_log_prior_data.txt", n.chains = 2)
  # trace plots
  plot(results, file = paste0("Convergence_diagnostics/Conv_pollen_",  
                              selectedcountry,"_", selectedabun,"_", selectedtrait, "_", 
                              selectedpolmode,  ".pdf"))
  # gelman plots
  pdf(paste0("Convergence_diagnostics/Gelman_pollen_", 
             selectedcountry,"_", selectedabun,"_", selectedtrait, "_", 
             selectedpolmode,  ".pdf"))
  
  gelman.plot(results, ask = FALSE)
  dev.off()
  
  # save cwm values
  res <- as.data.frame(summary(results, vars = "cwm")) %>%
    # filter cwm values
    rownames_to_column("parameter") %>%
    filter(!str_detect(parameter, "zcwm")) %>% 
    mutate(sitename = sitenames, 
                    trait = selectedtrait,
                    country = selectedcountry,
                    pollination = selectedpolmode,
                    correction = selectedabun,
                    growthform = "allpft")
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", 
                      selectedcountry,"_",selectedabun,"_", selectedtrait, "_", 
                      selectedpolmode, ".rds"))
  
  # save standardized cwm values
  zres <- data.frame(summary(results, vars = "zcwm"),
                    sitename = sitenames, 
                    trait = selectedtrait,
                    country = selectedcountry,
                    pollination = selectedpolmode,
                    correction = selectedabun,
                    growthform = "allpft")
  saveRDS(zres, paste0("RDS_files/03_zCWM_estimates_pollen_", 
                      selectedcountry,"_",selectedabun,"_", selectedtrait, "_", 
                      selectedpolmode, ".rds"))
  }


# CWM pollen remove taxa ----
cwm_pol_remove <- function(selectedtrait, selectedcountry, selectedabun,
                       removedspec){
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
    filter(ntrait >= 10)
  
  pol <- dfPOL %>% 
    # select country
    filter(country == selectedcountry) %>%
    # filter species
    filter(!pollentaxon == removedspec) %>% 
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
  mean.tax.log <- list(chain1 = MeanLog*0.01, chain2 = MeanLog*100)
  sd.tax.log <- list(chain1 = SDLog*0.01, chain2 = SDLog*100)
  cwm <- list(chain1 = rep(min(log(Trait)), Nsites), 
              chain2 = rep(max(log(Trait)), Nsites))
  # run model
  results <- run.jags("CWM_log_prior_data.txt", n.chains = 2, monitor = "cwm")
  
  # summarize results
  res <- data.frame(summary(results),
                    sitename = sitenames,
                    trait = selectedtrait,
                    country = selectedcountry,
                    correction = selectedabun,
                    removed_spec = removedspec)
  saveRDS(res, paste0("RDS_files/03_CWM_estimates_pollen_", selectedcountry,"_",
                      selectedtrait, "_", selectedabun,"_", removedspec, ".rds"))
}