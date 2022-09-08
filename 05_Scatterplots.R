## ---------------------------
##
## Script name: 04_Create_scatterplots
##
## Purpose of script: Create scatterplots to go in to the appendix
##
## Author: Annegreet Veeken
##
## Date Created: 2022-09-08
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

## Load data
dfCWM_veg <- readRDS("RDS_files/04_CWM_estimates_vegetation.rds")
dfCWM_pol <- readRDS("RDS_files/04_CWM_estimates_pollen.rds")

# join data

# filter draws 

# Produce plots
(p_scatter <-
    ggplot(cwm, aes(x = pollen_mean, y = veg_mean, color = zone)) +
    geom_point() +
    scale_fill_manual() +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    facet_wrap(country~trait, scales = "free") +
    ggtitle("Corrected pollen data") +
    theme_bw() +
    theme(legend.position =  "bottom"))

