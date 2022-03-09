## ---------------------------
##
## Script name: Correlation 
##
## Purpose of script: Calculating correlation between CWM species based and pollen based
##
## Author: Annegreet Veeken
##
## Date created: 08-02-2022
## ---------------------------
##
## Notes:
##   - re-run when leafC and leafN are available
##
## ---------------------------

## Load libraries
library(tidyverse)
library(rjags) 
library(runjags)
library(mcmc)
library(coda)

## Load data and prepare
dfCWM <- readRDS("RDS_files/03_Posterior_CWM.rds") %>% 
  as_tibble() %>% 
  filter(!trait %in% c("LeafC", "LeafN")) # not yet available

# Make list per trait, make tidy
lCWM <- dfCWM %>% 
  
  split(f = as.factor(.$trait)) %>% 
  purrr::map(., ~pivot_wider(., id_cols = sitename, names_from = zone ,
                             values_from = Mean)) %>% 
  keep(~ncol(.) == 5) 
traitnames <- names(lCWM)

## Plots ----
plot_func <- function(df,traitname){
p <- ggplot(df, aes(x = pollen, y = speciesA)) +
  geom_point(col = "darkorange") +
  geom_point(data = df, aes(x = pollen, y = speciesB), col = "darkorchid") +
  geom_point(data = df, aes(x = pollen, y = speciesC), col = "cyan4") +
  scale_x_continuous("Pollen-based CWM", limits = c(min(df[2:5]),max(df[2:5]))) +
  scale_y_continuous("Vegetation CWM", limits = c(min(df[2:5]),max(df[2:5]))) +
  ggtitle(traitname) +
  theme_bw()
}

pall <- purrr::map2(lCWM,traitnames, ~plot_func(df = .x, traitname = .y))
#

## Regression bayesian----

