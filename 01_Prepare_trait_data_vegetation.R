## ---------------------------
##
## Script name: Prepare trait data
##
## Purpose of script: Cleaning and tidying, and species nomenclature harmonisation
##
## Author: Annegreet Veeken
##
## Last modified: 11-03-2022
## ---------------------------
##
## Notes:
##   - make plots of traits per pollen taxa
##   - check units of the trait data and the gapfilled data
##
## ---------------------------

# Load libraries
library(tidyverse)
library(readxl)

## Load data ----
plantheight_raw <- read_xlsx("Data/Plant_height.xlsx")
sla_ldmc_raw <- read_xlsx("Data/LDMC_SLA_dataset.xlsx")
harm <- read_xlsx("Data/Harmonization_table_species_names.xlsx") %>% 
  rename(family = fam) %>% 
  mutate(stand.spec = str_remove(stand.spec, " sp."))

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
  dplyr::select(sitename = Site, species = Species, flowering,
         measurementID = measurement, PlantHeight) %>% 
  # get rid of the duplicate rows introduced by separate_rows (measurement id avoids deletion of actual duplicates)
  distinct(sitename, species, measurementID, PlantHeight, .keep_all = TRUE) %>% 
  # merge with species harmonisation list
  left_join(harm, by = "species") # check if all are species are available in list

saveRDS(plantheight, "RDS_files/01_PlantHeight.rds")

## SLA ----
sla <- sla_ldmc_raw %>% 
  # Harmonize column names
  dplyr::select(sitename = Site, species = Species,
          SLA = 'SLA (mm^2 mg−1)') %>% 
  # create measurement ID
  group_by(sitename, species) %>% 
  mutate(measurementID = row_number()) %>% 
  # sort
  arrange(sitename, species, measurementID) %>% 
  # merge with species harmonisation list
  left_join(harm, by = "species") %>% 
  filter(!is.na(SLA)) # check notes from raw data file for reason why missing
saveRDS(sla, "RDS_files/01_SLA.rds")

## LA ----
la <- sla_ldmc_raw %>% 
  dplyr::select(sitename = Site, species = Species,
         LA = "Leaf_area (mm^2)") %>% 
  # create measurement ID
  group_by(sitename, species) %>% 
  mutate(measurementID = row_number()) %>% 
  # sort
  arrange(sitename, species, measurementID) %>% 
  # merge with species harmonisation list
  left_join(harm, by = "species") %>% 
  filter(!is.na(LA)) %>%  # check notes from raw data file for reason why missing
  # convert to cm2
  mutate(LA = LA/100)
saveRDS(la, "RDS_files/01_LA.rds")
    
## Summary table per pollentaxon
bdm_traits <- readRDS("RDS_files/01_Traits_Swiss.rds")
spec <- readRDS("RDS_files/02_PollenType_species.rds")

sum_ph <- plantheight[,c("stand.spec", "PlantHeight")] %>% 
  left_join(spec, by = c("species" = "stand.spec")) %>% 
  filter(!is.na(pollentaxon)) %>% 
  group_by(pollentaxon) %>% 
  summarise(nspec = length(unique(stand.spec)),
            nobs = n(),
            LA = paste(round(mean(LA), 1), "+-", round(sd(LA), 1)))
sum_la <- la %>% 
  left_join(spec, by = c("species" = "stand.spec")) %>% 
  filter(!is.na(pollentaxon)) %>% 
  group_by(pollentaxon) %>% 
  summarise(nspec = length(unique(stand.spec)),
            nobs = n(),
            LA = paste(round(mean(LA), 1), "+-", round(sd(LA), 1)))
sum_sla <- sla %>% 
  left_join(spec, by = c("species" = "stand.spec")) %>% 
  filter(!is.na(pollentaxon)) %>% 
  group_by(pollentaxon) %>% 
  summarise(nspec = length(unique(stand.spec)),
            nobs = n(),
            SLA = paste(round(mean(SLA), 1), "+-", round(sd(SLA), 1)))

