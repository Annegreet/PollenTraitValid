gap_t## ---------------------------
##
## Script name: 02_prepare_trait_data_gapfilled.R
##
## Purpose of script: Prepare gapfilled data for CWM calculation
##
## Author: Annegreet Veeken
##
## Date Created: 2022-07-19
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
if (!require(LCVP)) install.packages("LCVP")
if (!require(lcvplants)) install.packages("lcvplants")

## Add pollen type to gapdata
# Species per pollen type list
spec_pollen <- readRDS("RDS_files/02_PollenType_species.rds")
spec <- readRDS("RDS_files/01_Species_abundance_bc.rds") %>% 
  pull(stand.spec) %>% 
  c(., spec_pollen) %>% 
  unlist() %>% 
  unique

# Gap-filled trait data
gapdata <-  
  read.csv("C:/Users/lgxgv/OneDrive - The University of Nottingham/Data analysis/Data/GapfilledTraitDataFS.csv") 
gapdata <- gapdata %>% 
    # Clean species names
    mutate(Species = 
           str_replace(Spp, "\\.", " ") %>%  # replace first . occurrence by space
           str_replace("\\.", "-") %>% # replace second . occurrence by - for species names with hyphens
             str_remove(" sp")) # remover sp epiphet 
# check for missing species in the trait data
misspec <-
  spec[!spec %in% gapdata$Species] %>%
  sort %>% unique

# standardize species names gapfilled data
# gapdat_sp <- gapdata$Species %>% 
#   str_subset(pattern = "sp\\.", negate = TRUE) %>% 
#   unique() %>% 
#   lcvp_search(., progress_bar = TRUE)
# saveRDS(gapdat_sp, "RDS_files/02_LCVP_search_gapdata.rds")
gapdat_sp <- readRDS("RDS_files/02_LCVP_search_gapdata.rds")
synonyms <- gapdat_sp %>%
  mutate(
    Input.Taxon = word(Input.Taxon, 1, 2) %>% # keep binomial without authority
      str_remove("_x"),# x added for unknown reason
    Output.Taxon = if_else(!is.na(Output.Taxon), word(Output.Taxon, 1, 2) %>% str_remove("_x"),
    Input.Taxon),
    genus = word(Output.Taxon,1)
  ) %>%
  distinct() %>% 
  filter(!is.na(Output.Taxon)) %>%
  select(Species = Input.Taxon, family = Family, genus, stand.spec = Output.Taxon) %>% 
  # manually add alternatives that are missing still
  add_row(
    family = c("Asteraceae", "Asteraceae", "Poaceae", "Asteraceae", "Rosaceae",
               "Ranunculaceae", "Asteraceae", "Asteraceae", "Caryophyllaceae",
               "Asteraceae","Asteraceae", "Rosaceae", "Ranunculaceae",
               "Ericaceae","Caryophyllaceae","Caryophyllaceae","Convolvulaceae",
               "Apiaceae","Geraniaceae", "Asteraceae", "Apiaceae","Apiaceae",
               "Asteraceae","Asteraceae", "Apiaceae"),
    genus = c("Senecio", "Senecio", "Spartina", "Anthemis", "Potentilla", "Ranunculus", 
              "Chrysanthemum", "Mycelis", "Stellaria", "Leontodon", "Taraxacum", "Acaena", 
              "Hepatica", "Arctostaphylos", "Silene", "Silene", "Calystegia", "Peucedanum",
              "Geranium", "Picris", "Apium", "Caucalis", "Hieracium", "Leontodon", "Tommasini"),
    Species = c("Jacobaea vulgaris","Jacobaea aquatica","Sporobolus anglicus",
                    "Cota austriaca","Drymocallis rupestris","Ficaria verna",
                    "Glebionis segetum","Lactuca muralis","Rabelera holostea",
                    "Scorzoneroides autumnalis","Taraxacum campylodes","Acaena novae-zeelandiae",
                    "Anemone hepatica","Arctostaphylos alpinus","Atocion armeria",
                    "Atocion rupestris","Calystegia sylvatica","Cervaria rivini",
                    "Geranium macrorhizum","Helminthotheca echioides",
                    "Helosciadium nodiflorum","Orlaya platycarpos","Pilosella piloselloides",
                    "Scorzoneroides autumnalis","Tommasinia verticillaris"),
    stand.spec = c("Senecio jacobaea","Senecio aquaticus","Spartina anglica",
                     "Anthemis austriaca","Potentilla rupestris","Ranunculus ficaria",
                     "Chrysanthemum segetum","Mycelis muralis","Stellaria holostea",
                     "Leontodon autumnalis","Taraxacum officinale","Acaena novae-zelandiae",
                     "Hepatica nobilis","Arctostaphylos alpina","Silene armeria",
                     "Silene rupestris","Calystegia silvatica","Peucedanum rivini",
                     "Geranium macrorrhizum","Picris echioides","Apium nodiflorum",
                     "Caucalis platycarpos","Hieracium piloselloides","Leontodon autumnalis",
                     "Tommasini verticillaris")
  ) %>% 
  distinct() 
 
# add standardized species name
gapdata_save <- gapdata %>% 
  left_join(synonyms, by = "Species") %>% 
  dplyr::select(family, genus, stand.spec, PlantHeight, SLA, LA = Area) %>% 
  # convert to scales comparable with the collected trait data
  mutate(across(where(is.numeric), exp)) %>% 
  mutate(PlantHeight = PlantHeight * 100, # m to cm
         LA = LA/100)  # mm2 to cm2  

# check for missing species in the trait data again
misspec <- 
  spec[!spec %in% gapdata_save$stand.spec] %>% 
  sort %>% unique
# checked for synonyms, but not all found (in Missing_trait.xlsx)

saveRDS(gapdata_save, "RDS_files/02_Gapfilled_traits.rds")

## traits merged with pollen trait table
# merge trait data with spec/pol list
dfTRAIT <- gapdata %>% 
  left_join(unique(synonyms[,c("Species", "stand.spec")]), by = "Species") %>% 
  filter(stand.spec %in% spec) %>% # only the relevant species
  left_join(spec_pollen, by = c("Species" = "stand.spec")) %>%  # attach pollen taxon
  arrange(pollentaxon) %>%  # sort alphabetically
  dplyr::select(family, genus, stand.spec = Species, pollentaxon, country, PlantHeight, SLA, LA = Area) %>% 
  # convert to scales comparable with the collected trait data
  mutate(across(where(is.numeric), exp)) %>% 
  mutate(PlantHeight = PlantHeight * 100, # m to cm
         LA = LA/100) %>% # mm2 to cm2  
  filter(!is.na(pollentaxon))

# Summary table
sum_trait <- dfTRAIT %>% 
  group_by(pollentaxon) %>% 
  summarise(nspec = length(unique(stand.spec)),
            nobs = n(),
            SLA = paste(round(mean(SLA), 1), "+-", round(sd(SLA), 2)),
            LA = paste(round(mean(LA), 1), "+-", round(sd(LA), 1)),
            `Plant height` = paste(round(mean(PlantHeight), 1), "+-", round(sd(PlantHeight), 1))) %>% 
  drop_na()

saveRDS(dfTRAIT, "RDS_files/02_Gapfilled_traits_pollen.rds")
saveRDS(sum_trait, "RDS_files/02_Summary_gapdata_pollen.rds")
