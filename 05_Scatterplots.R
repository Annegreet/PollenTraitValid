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
dfCWM_veg <- readRDS("RDS_files/04_CWM_estimates_vegetation_Scotland.rds") %>%
  mutate(across(c("veg_mean", "veg_sd"),~exp(.)))
dfCWM_pol <- readRDS("RDS_files/04_CWM_estimates_pollen_Scotland.rds") %>%
  filter(!str_detect(correction,pattern = "draw")) %>%
  mutate(across(c("pollen_mean", "pollen_sd"),~exp(.)))

# join data
dfCWM <- dfCWM_veg %>%
  left_join(dfCWM_pol, by = c("country","sitename", "trait", "growthform","pollination")) %>%
  # create treatment column
  mutate(treatment = case_when((pollination == "allpol" & growthform == "allpft" &
                                  taxres == "stand.spec") ~ "Correction factor",
                               taxres %in% c("family","genus") ~ "Taxonomic resolution",
                               (pollination != "allpol" & growthform == "allpft") ~ "Pollination mode",
                               (growthform != "allpft" & pollination == "allpol") ~ "Growth form"),
         zone = str_remove(zone, pattern = "Scotland |Switzerland ")) %>%
  filter(country == "Scotland")

# Plot labels
labs_trait <- c("Height (cm)",
          "Leaf area (mm2)",
          "SLA (mm2/mg)")
names(labs_trait) <- c("PlantHeight", "LA", "SLA")
labs_zone <- c(zoneA = "",
               zoneB = "",
               zoneC = "")

# TRY data pollen, field data vegetation ----
windows()
cwm_range <- dfCWM %>%
  group_by(trait) %>%
  summarise(across(c("veg_mean","pollen_mean"), ~range(.))) %>%
  transmute(max_lim = max(veg_mean, pollen_mean),
         min_lim = min(veg_mean, pollen_mean)) %>%
  pivot_longer(-trait, names_to = "limits") %>%
  mutate(pollen_mean = value, veg_mean = value,
         correction = NA, taxres = NA, pollination = NA, growthform = NA)

(p_correction <-
    dfCWM %>%
    filter(country == "Scotland") %>%
    filter(treatment == "Correction factor") %>%
    ggplot(aes(x = pollen_mean, y = veg_mean, color = correction)) +
    geom_point() +
    geom_blank(data = cwm_range) +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c("darkorange", "purple"),
                       label = c("Correction", "No correction")) +
    facet_wrap(zone~trait,
               scales = "free",
               labeller = labeller(trait = labs_trait,
                                    zone = labs_zone)
               ) +
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
    ggplot(aes(x = pollen_mean, y = veg_mean, color = taxres)) +
    geom_point() +
    geom_blank(data = cwm_range) +
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
    filter(country == "Scotland") %>%
    filter(treatment == "Growth form" ) %>%
    ggplot(aes(x = pollen_mean, y = veg_mean, color = growthform)) +
    geom_point() +
    geom_blank(data = cwm_range) +
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
    ggplot(aes(x = pollen_mean, y = veg_mean, color = pollination)) +
    geom_point() +
    geom_blank(data = cwm_range) +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c("darkorange", "purple"),
                       labels = c("Not wind pollinated", "Wind pollinated")) +
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

# TRY data for vegetation and pollen ----
dfCWM_veg <- readRDS("RDS_files/04_CWM_estimates_vegetation_TRY.rds") %>%
  filter(country == "Scotland")

dfCWM <- dfCWM_veg %>%
  left_join(dfCWM_pol, by = c("country","sitename", "trait", "growthform","pollination")) %>%
  # create treatment column
  mutate(treatment = case_when((pollination == "allpol" & growthform == "allpft" &
                                  taxres == "stand.spec") ~ "Correction factor",
                               taxres %in% c("family","genus") ~ "Taxonomic resolution",
                               (pollination != "allpol" & growthform == "allpft") ~ "Pollination mode",
                               (growthform != "allpft" & pollination == "allpol") ~ "Growth form"),
         zone = str_remove(zone, pattern = "Scotland |Switzerland ")) %>%
  filter(country == "Scotland")
cwm_range <- dfCWM %>%
  group_by(trait) %>%
  summarise(across(c("veg_mean","pollen_mean"), ~range(.))) %>%
  transmute(max_lim = max(veg_mean, pollen_mean),
            min_lim = min(veg_mean, pollen_mean)) %>%
  pivot_longer(-trait, names_to = "limits") %>%
  mutate(pollen_mean = value, veg_mean = value,
         correction = NA, taxres = NA, pollination = NA, growthform = NA)

windows()
(p_correction <-
    dfCWM %>%
    filter(treatment == "Correction factor") %>%
    ggplot(aes(x = pollen_mean, y = veg_mean, color = correction)) +
    geom_point() +
    geom_blank(data = cwm_range) +
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
ggsave("Figures/Scatter_correction_TRY.png", p_correction, height = 7, width = 7)

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
ggsave("Figures/Scatter_taxonomic_resolution_TRY.png", p_taxres,
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
ggsave("Figures/Scatter_growthform_TRY.png", p_pft,
       height = 7, width = 7)

(p_pol <-
    dfCWM %>%
    filter(treatment == "Pollination mode") %>%
    ggplot(aes(x = pollen_mean, y = veg_mean, color = pollination, shape = country)) +
    geom_point() +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c("darkorange", "purple"),
                       labels = c("Not wind pollinated", "Wind pollinated")) +
    facet_wrap(zone~trait, scales = "free",
               labeller = labeller(trait = labs_trait,
                                   zone = labs_zone)) +
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()))
ggsave("Figures/Scatter_pollination_TRY.png", p_pol,
       height = 7, width = 7)
