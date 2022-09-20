## ---------------------------
##
## Script name: 05_Scatterplots
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
if (!require(tidyverse)) install.packages("tidyverse")

## Load data
dfCWM_veg <- readRDS("RDS_files/04_CWM_estimates_vegetation.rds")
dfCWM_pol <- readRDS("RDS_files/04_CWM_estimates_pollen.rds") %>% 
  filter(!str_detect(correction,pattern = "draw"))

# join data
dfCWM <- dfCWM_veg %>% 
  left_join(dfCWM_pol, by = c("country","sitename", "trait", "growthform","pollination")) %>% 
  # create treatment column
  mutate(treatment = case_when((pollination == "allpol" & growthform == "allpft" &
                                  taxres == "stand.spec") ~ "Correction factor",
                               taxres %in% c("family","genus") ~ "Taxonomic resolution",
                               (pollination != "allpol" & growthform == "allpft") ~ "Pollination mode",
                               (growthform != "allpft" & pollination == "allpol") ~ "Growth form"),
         zone = str_remove(zone, pattern = "Scotland |Switzerland "))

# Produce plots
labs_trait <- c(PlantHeight = "Height (cm)",
          LA = "Leaf area (mm2)",
          SLA = "SLA (mg/mm2)")
labs_zone <- c(zoneA = "",
               zoneB = "",
               zoneC = "")

windows()
(p_correction <-
    dfCWM %>% 
    filter(treatment == "Correction factor") %>% 
    ggplot(aes(x = pollen_mean, y = veg_mean, color = correction, shape = country)) +
    geom_point() +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c("darkorange", "purple"),
                       label = c("Correction", "No correction")) +
    facet_wrap(zone~trait, scales = "free",
               labeller = labeller(trait = labs_trait,
                                   zone = labs_zone)) +
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank())) 
ggsave("Figures/Scatter_correction.png", p_correction, height = 7, width = 7)

(p_taxres <-
    dfCWM %>% 
    filter(treatment == "Taxonomic resolution" |
          (treatment == "Correction factor" & 
             correction == "no correction")) %>% 
    ggplot(aes(x = pollen_mean, y = veg_mean, color = taxres, shape = country)) +
    geom_point() +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c("darkorange", "purple","cyan4"),
                       labels = c("Family", "Genus", "Species")) +
    facet_wrap(zone~trait, scales = "free",
               labeller = labeller(trait = labs_trait,
                                   zone = labs_zone)) +
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()))
ggsave("Figures/Scatter_taxonomic_resolution.png", p_taxres, 
       height = 7, width = 7)

(p_pft <-
    dfCWM %>% 
    filter(treatment == "Growth form" ) %>% 
    ggplot(aes(x = pollen_mean, y = veg_mean, color = growthform, shape = country)) +
    geom_point() +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c("darkorange", "purple"),
                       labels = c("Non-woody", "Woody")) +
    facet_wrap(zone~trait, scales = "free",
               labeller = labeller(trait = labs_trait,
                                   zone = labs_zone)) +
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()))
ggsave("Figures/Scatter_growthform.png", p_pft, 
       height = 7, width = 7)

(p_pol <-
    dfCWM %>% 
    filter(treatment == "Pollination mode") %>% 
    ggplot(aes(x = pollen_mean, y = veg_mean, color = pollination, shape = country)) +
    geom_point() +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c("darkorange", "purple")) +
    facet_wrap(zone~trait, scales = "free",
               labeller = labeller(trait = labs_trait,
                                   zone = labs_zone)) +
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()))
ggsave("Figures/Scatter_pollination.png", p_pol, 
       height = 7, width = 7)
