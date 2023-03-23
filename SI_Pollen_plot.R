<<<<<<< HEAD
## ---------------------------
##
## Script name: SI_pollen_plot
##
## Purpose of script: Make pollen plot
##
## Author: Annegreet Veeken
##
## Date Created: 2022-12-19
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   
## References:
##
## ---------------------------

set.seed(1)
options(scipen = 9999)

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")

## Load data
dfPOL_Scot <- readRDS("RDS_files/01_Pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_Pollen_data_Swiss.rds") %>% 
  filter(!sitename == "unknown")
bdm_meta <- read_csv("Data/bdm_metadata.csv") %>% 
  mutate(sitename = as.character(aID_STAO))
dfPOL_Swiss <- left_join(dfPOL_Swiss, bdm_meta, by = "sitename")

# Plots
(p_scot <- dfPOL_Scot %>% 
    # tidy to allow grouped bar plot
    dplyr::select(sitename, pollentaxon, percent, adjustedpercent_mean) %>% 
    pivot_longer(cols = percent:adjustedpercent_mean, names_to = "pollen_correction", 
                 values_to = "percent") %>% 
    mutate(percent = percent*100) %>% 
    # only common pollen types
    group_by(pollentaxon) %>% 
    filter(any(percent > 5)) %>% 
    # plot
    ggplot(aes(x = sitename, y = percent, fill = pollen_correction)) +
      geom_bar(position = "dodge", stat = "identity") +
      coord_flip() +
      # scales
      scale_x_discrete("Site", labels = 16:1, limits = rev) +
      scale_y_continuous("Pollen percentage (%)", limits = c(0,100)) +
      scale_fill_manual(name = "", 
                       values = c("darkorchid", "darkorange"),
                       labels = c(`percent` = "No correction", adjustedpercent_mean = "Correction")) +
      # faceting
      facet_wrap(~pollentaxon, nrow = 1) +
      # theme
      theme_bw() +
      theme(text = element_text(size = 8),
            strip.background = element_rect(fill = "white"),
            legend.position = "top",
            legend.key.size = unit(0.1, "inch")))
ggsave("Figures/Pollen_percentage_scotland.png", p_scot, width = 7.5, height = 3)

(p_swiss <- dfPOL_Swiss %>% 
    # tidy to allow grouped bar plot
    dplyr::select(sitename, elevation, pollentaxon, percent, adjustedpercent_mean) %>% 
    pivot_longer(cols = percent:adjustedpercent_mean, names_to = "pollen_correction", 
                 values_to = "percent") %>% 
    mutate(percent = percent*100) %>% 
    # only common pollen types
    group_by(pollentaxon) %>% 
    filter(any(percent > 5)) %>% 
    # plot
    ggplot(aes(x = reorder(sitename, elevation), y = percent, fill = pollen_correction)) +
    geom_bar(position = "dodge", stat = "identity") +
    coord_flip() +
    # scales
    scale_x_discrete("Site", labels = 27:1, limits = rev) +
    scale_y_continuous("Pollen percentage (%)", limits = c(0,100)) +
    scale_fill_manual(name = "", 
                      values = c("darkorchid", "darkorange"),
                      labels = c(`percent` = "No correction", adjustedpercent_mean = "Correction")) +
    # faceting
    facet_wrap(~pollentaxon, nrow = 2) +
    # theme
    theme_bw() +
    theme(text = element_text(size = 8),
          strip.background = element_rect(fill = "white"),
          legend.position = "top",
          legend.key.size = unit(0.1, "inch")))
ggsave("Figures/Pollen_percentage_swiss.png", p_swiss, width = 7.5, height = 6)
=======
## ---------------------------
##
## Script name: SI_pollen_plot
##
## Purpose of script: Make pollen plot
##
## Author: Annegreet Veeken
##
## Date Created: 2022-12-19
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   
## References:
##
## ---------------------------

set.seed(1)
options(scipen = 9999)

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")

## Load data
dfPOL_Scot <- readRDS("RDS_files/01_Pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_Pollen_data_Swiss.rds") %>% 
  filter(!sitename == "unknown")
bdm_meta <- read_csv("Data/bdm_metadata.csv") %>% 
  mutate(sitename = as.character(aID_STAO))
dfPOL_Swiss <- left_join(dfPOL_Swiss, bdm_meta, by = "sitename")

# Plots
(p_scot <- dfPOL_Scot %>% 
    # tidy to allow grouped bar plot
    dplyr::select(sitename, pollentaxon, percent, adjustedpercent_mean) %>% 
    pivot_longer(cols = percent:adjustedpercent_mean, names_to = "pollen_correction", 
                 values_to = "percent") %>% 
    mutate(percent = percent*100) %>% 
    # only common pollen types
    group_by(pollentaxon) %>% 
    filter(any(percent > 5)) %>% 
    # plot
    ggplot(aes(x = sitename, y = percent, fill = pollen_correction)) +
      geom_bar(position = "dodge", stat = "identity") +
      coord_flip() +
      # scales
      scale_x_discrete("Site", labels = 16:1, limits = rev) +
      scale_y_continuous("Pollen percentage (%)", limits = c(0,100)) +
      scale_fill_manual(name = "", 
                       values = c("darkorchid", "darkorange"),
                       labels = c(`percent` = "No correction", adjustedpercent_mean = "Correction")) +
      # faceting
      facet_wrap(~pollentaxon, nrow = 1) +
      # theme
      theme_bw() +
      theme(text = element_text(size = 8),
            strip.background = element_rect(fill = "white"),
            legend.position = "top",
            legend.key.size = unit(0.1, "inch")))
ggsave("Figures/Pollen_percentage_scotland.png", p_scot, width = 7.5, height = 3)

(p_swiss <- dfPOL_Swiss %>% 
    # tidy to allow grouped bar plot
    dplyr::select(sitename, elevation, pollentaxon, percent, adjustedpercent_mean) %>% 
    pivot_longer(cols = percent:adjustedpercent_mean, names_to = "pollen_correction", 
                 values_to = "percent") %>% 
    mutate(percent = percent*100) %>% 
    # only common pollen types
    group_by(pollentaxon) %>% 
    filter(any(percent > 5)) %>% 
    # plot
    ggplot(aes(x = reorder(sitename, elevation), y = percent, fill = pollen_correction)) +
    geom_bar(position = "dodge", stat = "identity") +
    coord_flip() +
    # scales
    scale_x_discrete("Site", labels = 27:1, limits = rev) +
    scale_y_continuous("Pollen percentage (%)", limits = c(0,100)) +
    scale_fill_manual(name = "", 
                      values = c("darkorchid", "darkorange"),
                      labels = c(`percent` = "No correction", adjustedpercent_mean = "Correction")) +
    # faceting
    facet_wrap(~pollentaxon, nrow = 2) +
    # theme
    theme_bw() +
    theme(text = element_text(size = 8),
          strip.background = element_rect(fill = "white"),
          legend.position = "top",
          legend.key.size = unit(0.1, "inch")))
ggsave("Figures/Pollen_percentage_swiss.png", p_swiss, width = 7.5, height = 6)
>>>>>>> 25554d40a3bca064f942e712ca585f58d17d5ec9
