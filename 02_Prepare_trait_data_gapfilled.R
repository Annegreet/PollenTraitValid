## ---------------------------
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
if(!require(tidyverse)) install.packages("tidyverse")

## Add pollen type to gapdata
# Species per pollen type list
spec <- readRDS("RDS_files/02_PollenType_species.rds")

# Gap-filled trait data
gapdata <-  
  read.csv("C:/Users/lgxgv/OneDrive - The University of Nottingham/Data analysis/PhD/Meta-analysis/Analysis/Data/GapfilledTraitDataFS.csv") 
gapdata <- gapdata %>% 
  mutate(Species = 
           str_replace(Spp, "\\.", " ") %>%  # replace first . occurrence by space
           str_replace("\\.", "-")) # replace second . occurrence by -

# check for missing species in the trait data
misspec <-
  spec$stand.spec[!spec$stand.spec %in% gapdata$Species] %>%
  sort %>% unique

# look for alternative names of the trait data
# lcvp <- lcvp_fuzzy_search(misspec, max_distance = 2, progress_bar = TRUE)
# saveRDS(lcvp, "RDS_files/03_LCVP_search.rds")
lcvp <- readRDS("RDS_files/03_LCVP_search.rds")
synonyms <- lcvp %>%
  filter(Input.Taxon != Output.Taxon) %>%
  select(Input.Taxon, Output.Taxon) %>%
  mutate(
    Input.Taxon = word(Input.Taxon, 1, 2) %>% # keep binomial without authority
      str_remove("_x"),# x added for unknown reason
    Output.Taxon = word(Output.Taxon, 1, 2) %>% str_remove("_x")
  ) %>%
  distinct() %>% 
  filter(!is.na(Output.Taxon)) %>%
  # manually add alternatives that are missing still
  add_row(
    Input.Taxon = c("Jacobaea vulgaris","Jacobaea aquatica","Sporobolus anglicus",
                    "Cota austriaca","Drymocallis rupestris","Ficaria verna",
                    "Glebionis segetum","Lactuca muralis","Rabelera holostea",
                    "Scorzoneroides autumnalis","Taraxacum campylodes","Acaena novae-zeelandiae",
                    "Anemone hepatica","Arctostaphylos alpinus","Atocion armeria",
                    "Atocion rupestris","Calystegia sylvatica","Cervaria rivini",
                    "Geranium macrorhizum","Helminthotheca echioides",
                    "Helosciadium nodiflorum","Orlaya platycarpos","Pilosella piloselloides",
                    "Scorzoneroides autumnalis","Tommasinia verticillaris"),
    Output.Taxon = c("Senecio jacobaea","Senecio aquaticus","Spartina anglica",
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

# Change to synonym
regex_pattern <-
  setNames(synonyms$Output.Taxon,
           paste0("\\b", synonyms$Input.Taxon, "\\b"))

# check for missing species in the trait data again
misspec <- 
  spec$stand.spec[!spec$stand.spec %in% gapdata$Species] %>% 
  sort %>% unique
# checked for synonyms, but not all found (in Missing_trait.xlsx)

# check for duplicates
duplicats <- spec %>% 
  group_by(country, stand.spec) %>% 
  summarise(n = n()) %>% 
  filter(n>1)

# merge trait data with spec/pol list
dfTRAIT <- gapdata %>% 
  filter(Species %in% spec$stand.spec) %>% # only the relevant species
  left_join(spec, by = c("Species" = "stand.spec")) %>%  # attach pollen taxon
  arrange(pollentaxon) %>%  # sort alphabetically
  dplyr::select(family, genus, stand.spec = Species, pollentaxon, country, PlantHeight, SLA, LA = Area,
                LDMC, 
                LeafC = CperDryMass,
                LeafN = N) %>% 
  # convert to scales comparable with the collected trait data
  mutate(across(where(is.numeric), exp)) %>% 
  mutate(PlantHeight = PlantHeight * 100) # m to cm

# Summary table
sum_trait <- dfTRAIT %>% 
  group_by(pollentaxon) %>% 
  summarise(nspec = length(unique(stand.spec)),
            nobs = n(),
            SLA = paste(round(mean(SLA), 1), "+-", round(sd(SLA), 2)),
            LA = paste(round(mean(LA), 1), "+-", round(sd(LA), 1)),
            `Plant height` = paste(round(mean(PlantHeight), 1), "+-", round(sd(PlantHeight), 1))) %>% 
  drop_na()

saveRDS(dfTRAIT, "RDS_files/02_Gapfilled_traits.rds")

rm(list=setdiff(ls(), "sum_trait"))