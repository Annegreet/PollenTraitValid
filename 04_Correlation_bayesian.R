## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Annegreet Veeken
##
## Date Created: 2022-07-24
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
if(!require(rjags)) install.packages("rjags")
if(!require(runjags)) install.packages("runjags")
if(!require(mvtnorm)) install.packages("mvtnorm") #to generate correlated data with rmvnorm.
if(!require(car)) install.packages("car") #To plot the estimated bivariate normal distribution.

# collate data
traits <- c("LA","SLA","PlantHeight")
# all taxa
files <- 
  list.files("RDS_files/") %>% 
  str_subset(pattern = "03_CWM_estimates_") %>% 
  str_subset(pattern = "adjusted", negate = TRUE) %>% 
  str_subset(paste(traits, collapse = "|"))

lab <- files %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(., "multivariate", negate = TRUE) %>% 
  str_remove(., "03_CWM_estimates_") %>% 
  str_remove(., ".rds")

folderpath.fun <- function(x)
{paste("RDS_files/", x, sep = "/")}

dfCWM <- files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(lab, ~readRDS(.) %>% 
                as.data.frame() %>% 
                rownames_to_column("parameter") %>% 
                filter(., str_detect(parameter, "cwm\\[")) %>% 
                dplyr::select(Mean, SD) %>%  
                mutate(label = .y,
                       id = seq(1:nrow(.))
                )
  ) %>% 
  bind_rows() %>% 
  mutate(trait = str_extract(label, pattern = paste(traits, collapse = "|")) %>% 
           str_remove_all("_"),
         zone = str_extract(label, pattern = "pollen|zoneA|zoneB|zoneC"),
         growthform = case_when(str_detect(label, pattern = "grassherbfernshrubtreeNA") ~ "all",
                                str_detect(label, pattern = "grassherbfern") ~ "herb",
                                str_detect(label, pattern = "shrubtree") ~ "trsh",
                                TRUE ~ "all"),
         pollination = case_when(str_detect(label, pattern = "windnot wind") ~ "all",
                                 str_detect(label, pattern = "not wind") ~ "not wind",
                                 str_detect(label, pattern = "treeshrub") ~ "wind",
                                 TRUE ~ "all"),
         taxres = case_when(str_detect(label, pattern = "stand.spec") ~ "species",
                            str_detect(label, pattern = "genus") ~ "genus",
                            str_detect(label, pattern = "genus") ~ "family",
                            zone == "pollen" ~ "pollentaxon",
                            TRUE ~ "stand.spec"),
         country = str_extract(label, pattern = "Scotland|Switzerland")
  ) 


Trait <- dfCWM %>%
  filter(label %in% c("pollen_Scotland_PlantHeight_percent",
                      "zoneB_Scotland_PlantHeight_stand.spec_windNAnot wind_treeshrubNAherbgrassfern")) %>%
  dplyr::select(label, id, Mean, SD) %>%
  pivot_wider(names_from = label, values_from = c(Mean,SD)) %>% 
  rename(pollen_mean = 2, zoneA_mean =3, pollen_sd=4, zoneA_sd=5)

CWM_mean <- Trait %>% 
  dplyr::select(pollen_mean, zoneA_mean) %>% 
  as.matrix()
CWM_sd <- Trait %>% 
  dplyr::select(pollen_sd, zoneA_sd) %>% 
  as.matrix()
N <- nrow(CWM_mean)

# https://www.sumsar.net/blog/2013/08/bayesian-estimation-of-correlation/
model_string <- "
  model {
    for(i in 1:N) {
      CWM_mean[i,1:2] ~ dmnorm(mu[], prec[,,i])
    }

    # Constructing the covariance matrix and the corresponding precision matrix.
    for(i in 1:N){
      prec[1:2,1:2,i] <- inverse(cov[,,i])
      cov[1,1,i] <- CWM_sd[i,1] * CWM_sd[i,1]
      cov[1,2,i] <- CWM_sd[i,1] * CWM_sd[i,2] * rho
      cov[2,1,i] <- CWM_sd[i,1] * CWM_sd[i,2] * rho
      cov[2,2,i] <- CWM_sd[i,2] * CWM_sd[i,2]
    }

    # Priors
    rho ~ dunif(-1, 1)
    mu[1] ~ dnorm(0, 0.001)
    mu[2] ~ dnorm(0, 0.001)

    # Generate random draws from the estimated bivariate normal distribution
    for(i in 1:N){
      x_rand[i,1:2] ~ dmnorm(mu[], prec[,,i])
      }
    #data# N, CWM_mean, CWM_sd
    #monitor# rho
  }"

inits_list <-  list(mu = c(0, 10),
                  rho = cor(CWM_mean[, 1], CWM_mean[, 2]))

results <- run.jags(model_string, inits = inits_list)

plot(Trait$pollen_mean, Trait$zoneA_mean)
