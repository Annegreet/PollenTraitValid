## ---------------------------
##
## Script name: 01_Prepare_TRY_data
##
## Purpose of script: 
## - Prepare downloaded TRY data for further analysis
## - for trait data to reconstruct trait composition from pollen data
## - for pollination mode and plant functional type to classify vegetation data
## 
## Author: Annegreet Veeken
##
## Date Created: 2022-08-15
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
## - based on a script by Franziska Schrodt   
##
## ---------------------------

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(data.table)) install.packages("data.table")
if (!require(LCVP)) install.packages("LCVP")
if (!require(lcvplants)) install.packages("lcvplants")

# increase memory limit
memory.limit(999999)

## Traits ----
trait_raw <-
  fread(
    "Data/TRY-request/22284_16082022183406/22284.txt",
    sep = "\t",
    data.table = FALSE,
    stringsAsFactors = FALSE,
    strip.white = TRUE,
    drop = c("LastName", "FirstName", "Dataset", "Comment", 
             "V28", "Reference", "ValueKindName", "Replicates"))

#removing trait outliers based on error risk
trait <- trait_raw %>%
  dplyr::select(DatasetID,
          DataID,
          ObsDataID,
          ObservationID,
          AccSpeciesID,
          AccSpeciesName,
          TraitID,
          OriglName,
          TraitName,
          OrigObsDataID,
          OrigValueStr,
          OrigUnitStr,
          StdValue,
          UnitName,
          ErrorRisk
          ) %>%
  mutate(ErrorRisk2 = ifelse(is.na(ErrorRisk), 0, ErrorRisk)) %>%
  filter(ErrorRisk2 < 4) # error risks based on Z-scores, see TRY 2019 paper

# Get meta data
meta <- trait %>% 
  filter(is.na(TraitID))
# Filter out meta data
trait <- trait %>% 
  filter(!is.na(TraitID))

#generate list of units for ALL TRY traits
units <- trait %>%
  select(OriglName,
         OrigUnitStr, 
         TraitName, 
         UnitName) %>%
  unique()

# Subsetting out traits and naming them
# if origlname or unit is blank but there are no duplicates we are keeping it. 
trait2 <- trait %>%
  mutate(remove = case_when(TraitID == 14 & OriglName == "N amount%" ~ 1, 
                            TraitID == 14 & OriglName == "N_senesced_leaf" ~ 1,
                            TraitID == 15 & OriglName == "P_senesced_leaf" ~ 1, 
                            TraitID == 47 & OriglName == "LDMC_min" ~ 1,
                            TraitID == 47 & OriglName == "LDMC_max" ~ 1,
                            TraitID == 47 & OriglName == "WCf" ~ 1,
                            TraitID == 47 & OriglName == "Leaf dry matter concentration predicted from NIRS" ~ 1,
                            TraitID == 3107 & OriglName == "Plant_height_generative_max" ~ 1,
                            TraitID == 3107 & OriglName == "Plant_height_generative_min" ~ 1,
                            TraitID == 3110 & OriglName == "Leaf_area_min" ~ 1,
                            TraitID == 3110 & OriglName == "Leaf_area_max" ~ 1,
                            TraitID == 3116 & OriglName == "SLA_min" ~ 1,
                            TraitID == 3116 & OriglName == "SLA_max" ~ 1,
                            TraitID == 3116 & UnitName == '' ~ 1,
                            TRUE ~ 0)) %>% 
  filter(remove == 0) %>%
  select(-remove) %>%
  mutate(CleanTraitName = case_when(TraitID == 47 ~ "LDMC",
                                    TraitID == 3107 ~ "PlantHeight", 
                                    TraitID == 3116 ~ "SLA", 
                                    TraitID == 3110 ~ "LA",
                                    TraitID == 14 ~ "LeafN",
                                    TraitID == 15 ~ "LeafP")) %>%
  filter(!is.na(StdValue))


# subset out dead plants
# get list of dead plants
deadID <- meta %>%
  filter(str_detect(OrigValueStr, "Dead|dead")) %>%
  pull(ObservationID) %>% 
  unique() 

#subset plant that were not measured in natural conditions
#get list of observations that were not in natural settings
setting <- meta %>%
  mutate(drop = ifelse(OrigValueStr == "Botanical garden" |
                         OrigValueStr == "botanical garden (Bergius Botanical Garden, Stockholm, Sweden)" |
                         OrigValueStr == "Botanical gardens, greenhouses and other atypical habitats" |
                         OrigValueStr == "Chamber" |
                         OrigValueStr == "climate chamber" |
                         OrigValueStr == "Climate chamber" |
                         OrigValueStr == "Climate Chamber" |
                         OrigValueStr == "Climate chamber, non-limiting conditions, (cf. dataset reference)" |
                         OrigValueStr == "climate chambers" |
                         OrigValueStr == "Common Garden" |
                         OrigValueStr == "Controlled climate chamber" |
                         OrigValueStr == "controlled environment room" |
                         OrigValueStr == "drought treatment" |
                         OrigValueStr == "FACE" |
                         OrigValueStr == "FE" |
                         OrigValueStr == "C" |
                         OrigValueStr == "Field Experiment" |
                         OrigValueStr == "FW" |
                         OrigValueStr == "G" |
                         OrigValueStr == "GH" |
                         OrigValueStr == "Glasshouse" |
                         OrigValueStr == "Greehouse" |
                         OrigValueStr == "Green house" |
                         OrigValueStr == "greenhouse" |
                         OrigValueStr == "Greenhouse" |
                         OrigValueStr == "Greenhouse, grrowth container" |
                         OrigValueStr == "groth chamber" |
                         OrigValueStr == "growth-chamber" |
                         OrigValueStr == "growth chamber" |
                         OrigValueStr == "Growth chamber" |
                         OrigValueStr == "Growth Chamber" |
                         OrigValueStr == "growth chambers" |
                         OrigValueStr == "Growth chambers" |
                         OrigValueStr == "Growth exp" |
                         OrigValueStr == "hydroponic" |
                         OrigValueStr == "Irrigation" |
                         OrigValueStr == "Irrigation and N fertilisation (100 kg/ha)" |
                         OrigValueStr == "LAU_Ploughed/mown" |
                         OrigValueStr == "LAU_Ploughed/mown and fertilized" |
                         OrigValueStr == "mesocosm" |
                         OrigValueStr == "mini-ecosystem" |
                         OrigValueStr == "N" |
                         OrigValueStr == "natural environment, high warming +4C, preccipitation ambient" |
                         OrigValueStr == "natural environment, high warming +4C, preccipitation ambient -50%" |
                         OrigValueStr == "natural environment, high warming +4C, preccipitation ambient +50%" |
                         OrigValueStr == "natural environment, low warming +1.5C, preccipitation ambient" |
                         OrigValueStr == "natural environment, low warming +1.5C, preccipitation ambient -50%" |
                         OrigValueStr == "natural environment, low warming +1.5C, preccipitation ambient +50%" |
                         OrigValueStr == "natural environment, medium warming +2.5C, preccipitation ambient" |
                         OrigValueStr == "natural environment, medium warming +2.5C, preccipitation ambient -50%" |
                         OrigValueStr == "natural environment, medium warming +2.5C, preccipitation ambient +50%" |
                         OrigValueStr == "natural environment, no warming, preccipitation ambient -50%" |
                         OrigValueStr == "natural environment, no warming, preccipitation ambient +50%" |
                         OrigValueStr == "natural grassland, experimental nutrient NP addition" |
                         OrigValueStr == "nutrient addition experiment" |
                         OrigValueStr == "Open Top" |
                         OrigValueStr == "open-top chamber" |
                         OrigValueStr == "Open top chambers" |
                         OrigValueStr == "OTC" |
                         OrigValueStr == "plantation" |
                         OrigValueStr == "PM" |
                         OrigValueStr == "pot" |
                         OrigValueStr == "Pot-grown" |
                         OrigValueStr == "Pots outside" |
                         OrigValueStr == "pots, outside in natural environment" |
                         OrigValueStr == "shade houses" |
                         OrigValueStr == "university campus" |
                         OrigValueStr == "Uzbekistan: Irrigated desert land" |
                         OrigValueStr == "VER_permanent extensively mown meadow" |
                         OrigValueStr == "VER_permanent meadow mown and fertilized" |
                         OrigValueStr == "VER_permanent meadows mown and fertilized" |
                         OrigValueStr == "water stress experiment" |
                         OrigValueStr == "water treatment", 1, 0)) %>%
  select(ObsDataID, drop, ObservationID) %>%
  unique()

table(setting$drop)

settingID <- setting %>% 
  filter(drop == 1) %>% 
  pull(ObservationID) %>% 
  unique() 

#remove observations dead and not natural observation
trait3 <- trait2 %>% 
  filter(!ObservationID %in% c(deadID, settingID))

# dataset specific problems
#investigating problem traits
d453 <- trait3 %>% #this dataset has 3 obs per plant, but no way to link leaves so we are averaging
  filter(DatasetID == 453) %>%
  group_by(DatasetID, ObservationID, CleanTraitName) %>%
  summarise(StdValue = mean(StdValue))

# remove and replace data sets 453 and 428
trait4 <- trait3 %>%
  filter(!DatasetID == 453) %>%
  bind_rows(d453) %>% 
  dplyr::select(ObsDataID, DatasetID, OrigObsDataID, ObservationID, AccSpeciesName, 
                CleanTraitName, StdValue) %>% 
  pivot_wider(id_cols = c(ObsDataID, DatasetID, OrigObsDataID,ObservationID, 
                          AccSpeciesName),
              names_from = CleanTraitName, values_from = StdValue) %>% 
  # remove duplicated values
  filter(!(duplicated(.[4:10]) & !is.na(OrigObsDataID))) %>% 
  # convert plantheight from m to cm
  mutate(PlantHeight = PlantHeight * 100) %>% 
  # remove zeros from the measurements
  mutate(across(LeafN:LA, ~na_if(.,0)))

# add pollen to species translation table ----
spec <- readRDS("RDS_files/02_PollenType_species.rds") %>% 
  distinct()

# check for missing species in the trait data
trait_sp <- trait4$AccSpeciesName %>% unique
pol_sp <- spec$stand.spec %>% unique()
misspec <-
  pol_sp[!pol_sp %in% trait_sp] %>%
  sort %>% unique

# improve matching by standardizing species names
# lcvp <- lcvp_search(str_subset(trait_sp, pattern = " "),
#                     progress_bar = TRUE) 
# saveRDS(lcvp, "RDS_files/01_TRY_species_standardization.rds")
lcvp <- readRDS("RDS_files/01_TRY_species_standardization.rds")
lcvp_stand <- lcvp %>% 
  dplyr::select(Search, Output.Taxon, family = Family) %>% 
  mutate(stand.spec = word(Output.Taxon, 1,2),
         genus = word(Output.Taxon, 1))
trait5 <- trait4 %>% 
  left_join(lcvp_stand, by = c("AccSpeciesName" = "Search")) 

# Save clean TRY data
saveRDS(trait5, "RDS_files/01_Clean_TRY_data.rds")
# trait5 <- readRDS("RDS_files/01_Clean_TRY_data.R")

# check missing species again
trait_sp <- trait5$stand.spec %>% unique
misspec <-
  pol_sp[!pol_sp %in% trait_sp] %>%
  sort %>% unique %>% length()/length(pol_sp) * 100
misspec

# join with pollen translation table
trait6 <- trait5 %>% 
  left_join(spec, by = "stand.spec")

saveRDS(trait6, "RDS_files/01_TRY_raw_traits.rds")
# trait6 <- readRDS("RDS_files/01_TRY_raw_traits.rds")

sum_trait <- trait6  %>% 
  filter(!is.na(pollentaxon)) %>% 
  group_by(pollentaxon) %>% 
  summarise(nspec = length(unique(stand.spec)),
            nobs = across(c(SLA,LA,PlantHeight),~sum(!is.na(.)))) %>% 
  unnest(nobs) %>% 
  filter(pollentaxon %in% c("Betula", "Carpinus betulus", "Cyperaceae", 
                            "Ericales (tetrad)", "Juniperus", "Pinus", 
                            "Poaceae", "Pteridophyte", "Ranunculaceae", 
                            "Alnus", "Caryophyllaceae", "Rosaceae", 
                            "Lamiaceae", "Rubiaceae", "Malvaceae", 
                            "Plantago", "Viola", "Apiaceae", "Asteraceae", 
                            "Picea", "Salix", "Geraniaceae", "Fraxinus", "Urtica",
                            "Ulmus", "Acer", "Larix", "Abies", "Corylus", "Epilobium",
                            "Quercus", "Rumex/Oxyria", "Castanea", "Frangula alnus", 
                            "Juglans", "Tilia", "Fabaceae", "Lonicera", "Solanum",
                            "Populus", "Veronica", "Convolvulaceae" ))


saveRDS(sum_trait,"RDS_files/01_summary_table_TRY_data.rds")          

trait_val <- trait6  %>% 
  filter(!is.na(pollentaxon)) %>% 
  group_by(pollentaxon) %>% 
  summarise(across(c(SLA,LA,PlantHeight), 
                   ~paste(round(mean(., na.rm = TRUE), 1), "+-", 
                          round(sd(., na.rm = TRUE), 2)))) %>% 
  filter(pollentaxon %in% c("Betula", "Carpinus betulus", "Cyperaceae", 
                          "Ericales (tetrad)", "Juniperus", "Pinus", 
                          "Poaceae", "Pteridophyte", "Ranunculaceae", 
                          "Alnus", "Caryophyllaceae", "Rosaceae", 
                          "Lamiaceae", "Rubiaceae", "Malvaceae", 
                          "Plantago", "Viola", "Apiaceae", "Asteraceae", 
                          "Picea", "Salix", "Geraniaceae", "Fraxinus", "Urtica",
                          "Ulmus", "Acer", "Larix", "Abies", "Corylus", "Epilobium",
                          "Quercus", "Rumex/Oxyria", "Castanea", "Frangula alnus", 
                          "Juglans", "Tilia", "Fabaceae", "Lonicera", "Solanum",
                          "Populus", "Veronica", "Convolvulaceae" ))
saveRDS(trait_val,"RDS_files/01_TRY_data_values.rds")          

        