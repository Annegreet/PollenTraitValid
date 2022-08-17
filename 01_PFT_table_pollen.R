## ---------------------------
##
## Script name: Create_PFT_table
##
## Purpose of script: Create pollination mode and pft table for pollen data
##
## Author: Annegreet Veeken
##
## Date Created: 2022-07-14
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   - pollination mode comes from Reitalu (2019) J of Ecology
##
## ---------------------------


## Load packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(readxl)) install.packages("readxl")

## Load data
dfPOL_Scot <- readRDS("RDS_files/01_Pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_Pollen_data_Swiss.rds")
dfPOL <- bind_rows(dfPOL_Scot, dfPOL_Swiss)
polmode <- read_xlsx("Data/Pollination_modes_Reitalu_2019.xlsx")

pollentaxa <- dfPOL %>% 
  pull(pollentaxon) %>% 
  unique()

# Harmonize nomenclature between Reitalu study and current
polmode <- polmode %>% 
  mutate(pollentaxon_reit = pollentaxon,
         pollentaxon = dplyr::recode(pollentaxon,  
                                     "Apiaceae undiff." = "Apiaceae",
                                     "Caryophyllaceae undiff." = "Caryophyllaceae",
                                     "Cyperaceae undiff." = "Cyperaceae",
                                     "Ericacea-type" = "Ericales (tetrad)",
                                     "Fraxinus excelsior" = "Fraxinus",
                                     "Geranium" = "Geraniaceae",
                                     "Pedicularis" = "Pedicularis palustris",
                                     "Ranunculuceae undiff." = "Ranunculaceae",
                                     "Rosaceae undiff." = "Rosaceae",
                                     "Viola palustris-type" = "Viola",
                                     "Corylus/Myrica" = "Corylus",
                                     "Drosera rotundifolia-type" = "Drosera",
                                     "Epilobium-type" = "Epilobium",
                                     "Fabaceae undiff." = "Fabaceae",
                                     "Juglans regia" = "Juglans",
                                     "Populus tremula" = "Populus",
                                     "Rhinanthus-type" = "Rhinanthus",
                                     "Rumex acetosa/acetosella-type" = "Rumex/Oxyria",
                                     "Solanum nigrum-type" = "Solanum",
                                     "Tilia cordata" = "Tilia",
                                     "Plantago undiff." = "Plantago"),
         fam = dplyr::recode(fam, 
                             "Compositae" = "Asteraceae",
                             "Leguminosae" = "Fabaceae")) %>% 
  add_row(pollentaxon = "Asteraceae", fam = "Asteraceae", growthform = "herb", pollination = "not wind") %>% 
  add_row(pollentaxon = "Lamiaceae", fam = "Lamiaceae", growthform = "herb", pollination = "not wind") %>% 
  add_row(pollentaxon = "Malvaceae", fam = "Malvaceae", growthform = "tree", pollination = "not wind") %>% 
  add_row(pollentaxon = "Liliaceae", fam = "Liliaceae", growthform = "herb", pollination = "not wind")  %>% 
  add_row(pollentaxon = "Pteridophyte", fam = "Pteridophyte", growthform = "herb", pollination = "wind") %>% 
  add_row(pollentaxon = "Tsuga", fam = "Pinaceae", growthform = "tree", pollination = "wind") %>% 
  mutate(species = case_when(str_detect(pollentaxon, " ") ~ pollentaxon) %>% str_remove("-type"),
    genus = case_when(!str_detect(pollentaxon, "eae| ") ~ pollentaxon) %>% str_remove("-type"),
    genus = if_else(is.na(genus), str_remove(species, "cf.") %>% word(1), genus))

saveRDS(polmode,"RDS_files/01_Pollination_mode.rds")