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
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(data.table)) install.packages("data.table")

## Traits ----
trait_raw <-
  fread(
    "Data/22274_15082022203634/22274.txt",
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
  filter(ErrorRisk2 < 8) 

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
  # remove problem data where there are replicates for ObservationID
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
  mutate(drop = ifelse(OrigValueStr == "Botanical garden"|
                         OrigValueStr == "botanical garden (Bergius Botanical Garden, Stockholm, Sweden)"|
                         OrigValueStr == "Botanical gardens, greenhouses and other atypical habitats"|
                         OrigValueStr == "Chamber"|
                         OrigValueStr == "climate chamber"|
                         OrigValueStr == "Climate chamber"|
                         OrigValueStr == "Climate Chamber"|
                         OrigValueStr == "Climate chamber, non-limiting conditions, (cf. dataset reference)"|
                         OrigValueStr == "climate chambers"|
                         OrigValueStr == "Common Garden"|
                         OrigValueStr == "Controlled climate chamber"|
                         OrigValueStr == "controlled environment room"|
                         OrigValueStr == "drought treatment"|
                         OrigValueStr == "FACE"|
                         OrigValueStr == "FE"|
                         OrigValueStr == "C"|
                         OrigValueStr == "Field Experiment"|
                         OrigValueStr == "FW"|
                         OrigValueStr == "G"|
                         OrigValueStr == "GH"|
                         OrigValueStr == "Glasshouse"|
                         OrigValueStr == "Greehouse"|
                         OrigValueStr == "Green house"|
                         OrigValueStr == "greenhouse"|
                         OrigValueStr == "Greenhouse"|
                         OrigValueStr == "Greenhouse, grrowth container"|
                         OrigValueStr == "groth chamber"|
                         OrigValueStr == "growth-chamber"|
                         OrigValueStr == "growth chamber"|
                         OrigValueStr == "Growth chamber"|
                         OrigValueStr == "Growth Chamber"|
                         OrigValueStr == "growth chambers"|
                         OrigValueStr == "Growth chambers"|
                         OrigValueStr == "Growth exp"|
                         OrigValueStr == "hydroponic"|
                         OrigValueStr == "Irrigation"|
                         OrigValueStr == "Irrigation and N fertilisation (100 kg/ha)"|
                         OrigValueStr == "LAU_Ploughed/mown"|
                         OrigValueStr == "LAU_Ploughed/mown and fertilized"|
                         OrigValueStr == "mesocosm"|
                         OrigValueStr == "mini-ecosystem"|
                         OrigValueStr == "N"|
                         OrigValueStr == "natural environment, high warming +4C, preccipitation ambient"|
                         OrigValueStr == "natural environment, high warming +4C, preccipitation ambient -50%"|
                         OrigValueStr == "natural environment, high warming +4C, preccipitation ambient +50%"|
                         OrigValueStr == "natural environment, low warming +1.5C, preccipitation ambient"|
                         OrigValueStr == "natural environment, low warming +1.5C, preccipitation ambient -50%"|
                         OrigValueStr == "natural environment, low warming +1.5C, preccipitation ambient +50%"|
                         OrigValueStr == "natural environment, medium warming +2.5C, preccipitation ambient"|
                         OrigValueStr == "natural environment, medium warming +2.5C, preccipitation ambient -50%"|
                         OrigValueStr == "natural environment, medium warming +2.5C, preccipitation ambient +50%"|
                         OrigValueStr == "natural environment, no warming, preccipitation ambient -50%"|
                         OrigValueStr == "natural environment, no warming, preccipitation ambient +50%"|
                         OrigValueStr == "natural grassland, experimental nutrient NP addition"|
                         OrigValueStr == "nutrient addition experiment"|
                         OrigValueStr == "Open Top"|
                         OrigValueStr == "open-top chamber"|
                         OrigValueStr == "Open top chambers"|
                         OrigValueStr == "OTC"|
                         OrigValueStr == "plantation"|
                         OrigValueStr == "PM"|
                         OrigValueStr == "pot"|
                         OrigValueStr == "Pot-grown"|
                         OrigValueStr == "Pots outside"|
                         OrigValueStr == "pots, outside in natural environment"|
                         OrigValueStr == "shade houses"|
                         OrigValueStr == "university campus"|
                         OrigValueStr == "Uzbekistan: Irrigated desert land"|
                         OrigValueStr == "VER_permanent extensively mown meadow"|
                         OrigValueStr == "VER_permanent meadow mown and fertilized"|
                         OrigValueStr == "VER_permanent meadows mown and fertilized"|
                         OrigValueStr == "water stress experiment"|
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
  dplyr::select(DatasetID, OrigObsDataID, ObservationID, AccSpeciesName, CleanTraitName, StdValue) %>% 
  pivot_wider(names_from = CleanTraitName, values_from = StdValue) %>% 
  # remove duplicated values
  filter(!(duplicated(.[4:10]) & !is.na(OrigObsDataID)))

saveRDS(trait4, "RDS_files/01_TRY_raw_traits.rds")

# Pollination and PFT ----
traits <- fread(
  "Data/21989_20072022132411/21989.txt",
  sep = "\t",
  data.table = FALSE,
  stringsAsFactors = FALSE,
  strip.white = TRUE,
  drop = c("LastName", "FirstName","Comment", 
           "V28", "Reference", "ValueKindName", "Replicates"))
spec <- readRDS("RDS_files/02_PollenType_species.rds")

pft_pol <- traits %>% 
  filter(TraitID %in% c(29, 42)) %>% 
  dplyr::select(AccSpeciesName, TraitName, OrigValueStr)

pft <- traits %>% 
  filter(Dataset %in% c("Reich-Oleksyn Global Leaf N, P Database",
                        " Plant growth form dataset for the New World",
                        "Categorical Plant Traits Database",
                        "PLANTSdata USDA")) %>%
  filter(TraitName == "Plant growth form") %>% 
  # fix duplicates and mistakes
  mutate(AccSpeciesName = recode(AccSpeciesName, 
                                 `VACCINIUM VITIS-IDAEA` = "Vaccinium vitis-idaea")) %>% 
  
  # filter(AccSpeciesName %in% c(spec$species, "VACCINIUM VITIS-IDAEA", "Senecio jacobea",
  #                              "Persicaria bistorta")) %>% 
  group_by(AccSpeciesName, Dataset) %>% 
  # group_by(AccSpeciesName) %>% 
  summarise(pft = paste(unique(OrigValueStr), collapse = ", ")) %>% 
  mutate(growthform = case_when(str_detect(pft, "tree|Tree") ~ 'tree',
                                str_detect(pft, "shrub|Shrub") ~ 'shrub',
                                str_detect(pft, "herb|Herb|Forb") ~ 'herb',
                                str_detect(pft, "grass|Grass|Graminoid|graminoid") ~ 'grass',
                                str_detect(pft, "fern|Fern") ~ 'fern')) %>% 
  filter(!is.na(growthform)) %>% 
  distinct(AccSpeciesName, growthform) %>% 
  filter(!(AccSpeciesName == "Achillea millefolium" & growthform == "shrub") ) %>% 
  filter(!(AccSpeciesName == "Agrostis capillaris" & growthform == "herb" )) %>% 
  filter(!(AccSpeciesName == "Anthoxanthum odoratum" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Arrhenatherum elatius" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Athyrium filix-femina" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Blechnum spicant" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Briza media" & growthform == "shrub")) %>% 
  filter(!(AccSpeciesName == "Carex digitata" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Carex montana" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Corylus avellana" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Crataegus laevigata" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Dactylis glomerata" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Dactylis glomerata" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Deschampsia flexuosa" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Digitalis purpurea" & growthform == "shrub")) %>% 
  filter(!(AccSpeciesName == "Dryopteris" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Empetrum nigrum" & growthform == "shrub")) %>% 
  filter(!(AccSpeciesName == "Equisetum sylvaticum" & growthform == "shrub")) %>% 
  filter(!(AccSpeciesName == "Galium rotundifolium" & growthform == "shrub")) %>% 
  filter(!(AccSpeciesName == "Equisetum sylvaticum" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Juncus bufonius" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Juniperus communis" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Lolium perenne" & growthform == "herb")) %>% 
  filter(!(str_detect(AccSpeciesName,"Luzula") & growthform == "herb")) %>% 
  filter(!(str_detect(AccSpeciesName,"Luzula") & growthform == "shrub")) %>% 
  filter(!(AccSpeciesName == "Luzula multiflora" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Milium effusum" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Nardus stricta" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Orthilia secunda" & growthform == "shrub")) %>% 
  filter(!(str_detect(AccSpeciesName, "Plantago") & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Poa nemoralis" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Pteridium aquilinum" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Pteridium aquilinum" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Ranunculus flammula" & growthform == "shrub")) %>% 
  filter(!(AccSpeciesName == "Rosa pendulina" & growthform == "shrub")) %>% 
  filter(!(AccSpeciesName == "Rubus caesius" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Selaginella selaginoides" & growthform == "herb")) %>% 
  filter(!(str_detect(AccSpeciesName, "Rumex") & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Stellaria media" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Stellaria media" & growthform == "shrub")) %>% 
  filter(!(AccSpeciesName == "Trichophorum cespitosum" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Trifolium repens" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Urtica dioica" & growthform == "shrub")) %>% 
  filter(!(AccSpeciesName == "Botrychium lunaria" & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Cardamine pratensis" & growthform == "shrub")) %>% 
  filter(!(str_detect(AccSpeciesName, "Eriophorum") & growthform == "herb")) %>% 
  filter(!(str_detect(AccSpeciesName, "Juncus") & growthform == "herb")) %>% 
  filter(!(AccSpeciesName == "Fragaria vesca" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Lonicera nigra" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Taraxacum campylodes" & growthform == "tree")) %>% 
  filter(!(AccSpeciesName == "Juncus effusus" & growthform == "tree"))

polmode <- traits %>% 
  filter(TraitName == "Pollination syndrome") %>% 
  filter(AccSpeciesName %in% c(spec$species, "VACCINIUM VITIS-IDAEA", "Senecio jacobea",
                               "Persicaria bistorta")) %>% 
  group_by(AccSpeciesName) %>% 
  # group_by(AccSpeciesName) %>% 
  summarise(polination = paste(unique(OrigValueStr), collapse = ", ")) %>% 
  mutate(polmode = case_when(str_detect(polination, "wind") ~ "wind",
                             !str_detect(polination, "wind") ~ "not wind")) %>% 
  distinct(AccSpeciesName, polmode)

df <- full_join(pft, polmode, by = "AccSpeciesName") %>% 
  full_join(spec, by = c("AccSpeciesName" = "species")) %>% 
  mutate(polmode = case_when(AccSpeciesName == "Poaceae" ~ "wind",
                             AccSpeciesName == "Asteraceae" ~ "not wind",
                             is.na(polmode) & growthform == "fern" ~ "wind",
                             is.na(polmode) & family == "Asteraceae" ~ "not wind",
                             AccSpeciesName == "Vaccinium vitis-idaea" ~ "not wind",
                             TRUE ~ polmode),
         genus = case_when(AccSpeciesName == "Poaceae" & !is.na(genus) ~ NA)) %>% 
  distinct()
saveRDS(df, "RDS_files/Polmode_pft_vegetation.rds")

