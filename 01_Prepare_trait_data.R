## ---------------------------
##
## Script name: Prepare trait data
##
## Purpose of script: Cleaning and tidying, and species nomenclature harmonisation
##
## Author: Annegreet Veeken
##
## Last modified: 08-02-2022
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

# Load libraries
library(tidyverse)
library(readxl)
library(Taxonstand)

## Load data ----
plantheight_raw <- read_xlsx("Data/Plant_height.xlsx")
sla_ldmc_raw <- read_xlsx("Data/LDMC_SLA_dataset.xlsx")
harm <- read_xlsx("Data/Harmonization_table_species_names.xlsx")

## Plant height ----
# make plant height data tidy
plantheight <- plantheight_raw %>% 
  mutate_at(vars(matches("Height")), as.numeric) %>% 
  # convert to long format
  pivot_longer(cols = contains("Height"), names_to = "measurement", 
               values_to = "PlantHeight") %>% 
  # tidy up measurement column
  mutate(measurement = str_remove(measurement, "Height_")) %>% 
  # remove empty rows
  filter(!is.na(PlantHeight)) %>%
  # create flowering column based on string
  separate_rows(Flowering) %>%
  mutate(flowering = 
           ifelse(Flowering == measurement, "Flowering","Not-flowering")) %>% 
  # harmonise column names
  select(sitename = Site, species = Species, flowering,
         measurementID = measurement, PlantHeight) %>% 
  # get rid of the duplicate rows introduced by separate_rows (measurement id avoids deletion of actual duplicates)
  distinct(sitename, species, measurementID, PlantHeight, .keep_all = TRUE) %>% 
  # merge with species harmonisation list
  left_join(harm, by = "species") # check if all are species are available in list

saveRDS(plantheight, "RDS_files/01_PlantHeight.rds")

## SLA ----
sla <- sla_ldmc_raw %>% 
  # Harmonize column names
  select(sitename = Site, species = Species,
          SLA = `SLA (mm^2 mg−1)`) %>% 
  # create measurement ID
  group_by(sitename, species) %>% 
  mutate(measurementID = row_number()) %>% 
  # sort
  arrange(sitename, species, measurementID) %>% 
  # merge with species harmonisation list
  left_join(harm, by = "species") %>% 
  filter(!is.na(SLA)) # check notes from raw data file for reason why missing
saveRDS(sla, "RDS_files/01_SLA.rds")

## LDMC ----
ldmc <- sla_ldmc_raw %>% 
  select(sitename = Site, species = Species,
         LDMC = "LDMC (mg g–1)") %>% 
  distinct(sitename, species, .keep_all = TRUE) %>%  # ldmc has one measurement per species per site
  # create measurement ID
  group_by(sitename, species) %>% 
  mutate(measurementID = row_number()) %>% 
  # sort
  arrange(sitename, species, measurementID) %>% 
  # merge with species harmonisation list
  left_join(harm, by = "species") %>% 
  filter(!is.na(LDMC)) # check notes from raw data file for reason why missing
saveRDS(ldmc, "RDS_files/01_LDMC.rds")

## LA ----
la <- sla_ldmc_raw %>% 
  select(sitename = Site, species = Species,
         LA = "Leaf_area (mm^2)") %>% 
  # create measurement ID
  group_by(sitename, species) %>% 
  mutate(measurementID = row_number()) %>% 
  # sort
  arrange(sitename, species, measurementID) %>% 
  # merge with species harmonisation list
  left_join(harm, by = "species") %>% 
  filter(!is.na(LA)) # check notes from raw data file for reason why missing

saveRDS(la, "RDS_files/01_LA.rds")
    