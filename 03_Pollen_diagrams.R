## ---------------------------
##
## Script name: Pollen diagrams
##
## Purpose of script: Plot pollen diagrams
##
## Author: Annegreet Veeken
##
## Date Created: 2022-03-24
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
if(!require(rioja)) install.packages("rioja") # plot pollen diagrams
if(!require(ggrepel)) install.packages("ggrepel") # repel point labels
if(!require(readxl)) install.packages("readxl") # repel point labels


## Prepare data ----
pollen <- readRDS("RDS_files/01_Pollen_data_Scot.rds")
polmode <- read_xlsx("Data/Pollination_modes_Reitalu_2019.xlsx")
polspec <- readRDS("RDS_files/02_PollenType_species.rds") %>% 
  select(pollentaxon, stand.spec)

# Add pollination mode and plant functional type
# Harmonize nomenclature
polmode <- polmode %>% 
  mutate(pollentaxon = recode(pollentaxon,  "Apiaceae undiff." = "Apiaceae",
                              "Caryophyllaceae undiff." = "Caryophyllaceae",
                              "Cyperaceae undiff." = "Cyperaceae",
                              "Ericacea-type" = "Ericales (tetrad)",
                              "Fraxinus excelsior" = "Fraxinus",
                              "Geranium" = "Geraniaceae",
                              "Juniperus" = "Juniperus communis",
                              "Pedicularis" = "Pedicularis palustris",
                              "Ranunculuceae undiff." = "Ranunculaceae",
                              "Rosaceae undiff." = "Rosaceae",
                              "Viola palustris-type" = "Viola")) %>% 
  add_row(pollentaxon = "Asteraceae", fam = "Asteraceae", growthform = "herb", pollination = "not wind") %>% 
  add_row(pollentaxon = "Lamiaceae", fam = "Lamiaceae", growthform = "herb", pollination = "not wind") %>% 
  add_row(pollentaxon = "Malvaceae", fam = "Malvaceae", growthform = "tree", pollination = "not wind") %>% 
  add_row(pollentaxon = "Plantago", fam = "Plantaginaceae", growthform = "herb", pollination = "wind") %>% 
  add_row(pollentaxon = "Liliaceae")  %>% 
  add_row(pollentaxon = "Pteridophyte") %>% 
  add_row(pollentaxon = "Tsuga ")

## Plot pollen data ----
pol <- pollen %>% 
  ungroup %>% 
  select(sitename, pollentaxon, percent) %>% 
  pivot_wider(names_from = pollentaxon, values_from = percent) %>% 
  mutate(across(where(is.numeric), ~.*100)) %>% 
  select(where(~ any(. > 1))) %>% # filter low values
  as.data.frame()
sitenames <- pollen %>% pull(sitename)  %>% unique
bar_colour <- polmode %>% 
  filter(pollentaxon %in% colnames(pol)) %>% 
  mutate(colour_gf = case_when(growthform == "tree" ~ "forestgreen",
                               growthform == "herb" ~ "gold2",
                               growthform == "grass" ~ "blue",
                               is.na(growthform) ~ "grey")) %>% 
  arrange(pollentaxon) %>% 
  distinct()

par(fig=c(0.01, 1, 0.07, 0.8),
    mar = c(0,0,4,0))

pol_plot <-
  strat.plot(
    pol,
    plot.bar = TRUE,
    lwd.bar = 10,
    scale.percent = TRUE,
    plot.line = FALSE,
    xSpace = 0.01,
    x.pc.lab = TRUE,
    x.pc.omit0 = TRUE,
    las = 2,
    col.bar = bar_colour$colour_gf,
    yvar = 1:16,
    y.tks = 1:16,
    y.tks.labels = sitenames
  )

## Prepare data - Switzerland ----
pollen <- readRDS("RDS_files/01_Pollen_data_Swiss.rds")
polmode <- read_xlsx("Data/Pollination_modes_Reitalu_2019.xlsx")
polspec <- readRDS("RDS_files/02_PollenType_species.rds") %>% 
  select(pollentaxon, stand.spec)

# Add pollination mode and plant functional type
# Harmonize nomenclature
polmode <- polmode %>% 
  mutate(pollentaxon = recode(pollentaxon,  "Apiaceae undiff." = "Apiaceae",
                              "Caryophyllaceae undiff." = "Caryophyllaceae",
                              "Cyperaceae undiff." = "Cyperaceae",
                              "Ericacea-type" = "Ericales (tetrad)",
                              "Fraxinus excelsior" = "Fraxinus",
                              "Geranium" = "Geraniaceae",
                              "Juniperus" = "Juniperus communis",
                              "Pedicularis" = "Pedicularis palustris",
                              "Ranunculuceae undiff." = "Ranunculaceae",
                              "Rosaceae undiff." = "Rosaceae",
                              "Viola palustris-type" = "Viola",
                              "Epilobium-type" = "Epilobium")) %>% 
  add_row(pollentaxon = "Asteraceae", fam = "Asteraceae", growthform = "herb", pollination = "not wind") %>% 
  add_row(pollentaxon = "Lamiaceae", fam = "Lamiaceae", growthform = "herb", pollination = "not wind") %>% 
  add_row(pollentaxon = "Malvaceae", fam = "Malvaceae", growthform = "tree", pollination = "not wind") %>% 
  add_row(pollentaxon = "Plantago", fam = "Plantaginaceae", growthform = "herb", pollination = "wind") %>% 
  add_row(pollentaxon = "Rumex/Oxyria", fam = "Plantaginaceae", growthform = "herb", pollination = "wind") %>% 
  add_row(pollentaxon = "Corylus", fam = NA, growthform = "shrub", pollination = "wind") %>% 
  add_row(pollentaxon = "Juglans", fam = "Juglandaceae", growthform = "tree", pollination = "wind") %>%
  add_row(pollentaxon = "Populus", fam = "Salicaceae", growthform = "tree", pollination = "wind") %>%
  add_row(pollentaxon = "Liliaceae")  %>% 
  add_row(pollentaxon = "Pteridophyte") %>% 
  add_row(pollentaxon = "Tsuga ")

## Plot pollen data ----
pol <- pollen %>% 
  ungroup %>% 
  select(sitename, pollentaxon, percent) %>% 
  pivot_wider(names_from = pollentaxon, values_from = percent) %>% 
  mutate(across(where(is.numeric), ~.*100)) %>% 
  select(where(~ any(. > 5))) %>% # filter low values
  as.data.frame()
sitenames <- pollen %>% pull(sitename) %>% unique
bar_colour <- polmode %>% 
  filter(pollentaxon %in% colnames(pol)) %>% 
  mutate(colour_gf = case_when(growthform == "tree" ~ "forestgreen",
                               growthform == "shrub" ~ "forestgreen",
                               growthform == "herb" ~ "gold2",
                               growthform == "grass" ~ "blue",
                               is.na(growthform) ~ "grey")) %>% 
  arrange(pollentaxon) %>% 
  distinct()

par(fig=c(0.01, 1, 0.07, 0.8),
    mar = c(0,0,4,0))

pol_plot_swiss <-
  strat.plot(
    pol[,-1],
    plot.bar = TRUE,
    lwd.bar = 10,
    scale.percent = TRUE,
    plot.line = FALSE,
    xSpace = 0.01,
    x.pc.lab = TRUE,
    x.pc.omit0 = TRUE,
    las = 2,
    col.bar = bar_colour$colour_gf,
    yvar = 1:nrow(pol),
    y.tks = 1:nrow(pol),
    y.tks.labels = sitenames
  )
rm(list=setdiff(ls(), c("pol_plot","pol_plot_swiss")))
