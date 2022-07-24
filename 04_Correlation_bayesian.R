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
  filter(trait == "PlantHeight",
         growthform == "all",
         pollination == "all",
         taxres %in% c("pollentaxon", "species")) %>%
  filter(zone %in% c("pollen","zoneA"),
         country == "Scotland") %>%
  dplyr::select(zone, id, Mean) %>%
  pivot_wider(id_cols = c("id", "zone"), names_from = zone, values_from = Mean) 

plot(Trait)
N <- nrow(Trait)

set.seed(31415)

mu <- c(10, 30)
sigma <- c(20, 40)
rho <- -0.7
cov_mat <- rbind(c(sigma[1]^2, sigma[1]*sigma[2]*rho),
                 c(sigma[1]*sigma[2]*rho, sigma[2]^2))
N <- 30
Trait <- rmvnorm(N, mu, cov_mat)
plot(x, xlim=c(-125, 125), ylim=c(-100, 150))


# https://www.sumsar.net/blog/2013/08/bayesian-estimation-of-correlation/
model_string <- "
  model {
    for(i in 1:N) {
      Trait[i,1:2] ~ dmnorm(mu[], prec[ , ])
    }

    # Constructing the covariance matrix and the corresponding precision matrix.
    prec[1:2,1:2] <- inverse(cov[,])
    cov[1,1] <- sigma[1] * sigma[1]
    cov[1,2] <- sigma[1] * sigma[2] * rho
    cov[2,1] <- sigma[1] * sigma[2] * rho
    cov[2,2] <- sigma[2] * sigma[2]

    # Uninformative priors on all parameters which could, of course, be made more informative.
    sigma[1] ~ dunif(0, 1000)
    sigma[2] ~ dunif(0, 1000)
    rho ~ dunif(-1, 1)
    mu[1] ~ dnorm(0, 0.001)
    mu[2] ~ dnorm(0, 0.001)

    # Generate random draws from the estimated bivariate normal distribution
    x_rand ~ dmnorm(mu[], prec[ , ])

    #data# N, Trait
    #monitor# mu, rho, sigma, x_rand
  }"

inits_list <-  list(mu = c(mean(x[, 1]), mean(x[, 2])),
                  rho = cor(x[, 1], x[, 2]),
                  sigma = c(sd(x[, 1]), sd(x[, 1])))

run.jags(model_string, inits = inits_list)
