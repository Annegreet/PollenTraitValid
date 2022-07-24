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
spec <- spec  %>% 
  # missing species
  add_row(family = "Plantaginaceae", genus = "Digitalis", species = "Digitalis purpurea",
          count = NA, stand.spec = "Digitalis purpurea", country = "Scotland", pollentaxon = NA) %>% 
  add_row(family = "Plantaginaceae", genus = "Digitalis", species = "Digitalis purpurea",
          count = NA, stand.spec = "Digitalis purpurea", country = "Switzerland", pollentaxon = NA) %>% 
  add_row(family = "Primulaceae", genus = "Lysimachia", species = "Lysimachia nemorum",
          count = NA, stand.spec = "Lysimachia nemorum", country = "Scotland", pollentaxon = NA) %>% 
  add_row(family = "Primulaceae", genus = "Lysimachia", species = "Lysimachia nemorum",
          count = NA, stand.spec = "Lysimachia nemorum", country = "Switzerland", pollentaxon = NA) %>% 
  add_row(family = "Juncaceae", genus = "Luzula", species = "Luzula sylvatica",
          count = NA, stand.spec = "Luzula sylvatica", country = "Scotland", pollentaxon = NA) %>% 
  add_row(family = "Juncaceae", genus = "Luzula", species = "Luzula sylvatica",
          count = NA, stand.spec = "Luzula sylvatica", country = "Switzerland", pollentaxon = NA) %>% 
  add_row(family = "Brassicaceae", genus = "Cardamine", species = "Cardamine pratensis",
          count = NA, stand.spec = "Cardamine pratensis", country = "Scotland") %>% 
  add_row(family = "Brassicaceae", genus = "Cardamine", species = "Cardamine pratensis",
          count = NA, stand.spec = "Cardamine pratensis", country = "Switzerland") %>% 
  add_row(family = "Ophioglossaceae", genus = "Botrychium", species = "Botrychium lunaria",
          count = NA, stand.spec = "Botrychium lunaria", country = "Scotland") %>% 
  add_row(family = "Ophioglossaceae", genus = "Botrychium", species = "Botrychium lunaria",
          count = NA, stand.spec = "Botrychium lunaria", country = "Switzerland") %>% 
  add_row(family = "Droseraceae", genus = "Drosera", species = "Drosera rotundifolia",
          count = NA, stand.spec = "Drosera rotundifolia", country = "Scotland") %>% 
  add_row(family = "Droseraceae", genus = "Drosera", species = "Drosera rotundifolia",
          count = NA, stand.spec = "Drosera rotundifolia", country = "Switzerland") %>% 
  add_row(family = "Rosaceae", genus = "Alchemilla", species = "Alchemilla xanthochlora",
          count = NA, stand.spec = "Alchemilla xanthochlora", country = "Switzerland") %>% 
  add_row(family = "Rosaceae", genus = "Alchemilla", species = "Alchemilla xanthochlora",
          count = NA, stand.spec = "Alchemilla xanthochlora", country = "Scotland") %>% 
  add_row(family = "Cyperaceae", genus = "Carex", stand.spec= "", 
          country = "Scotland") %>% 
  add_row(family = "Juncaceae", genus = "Juncus", stand.spec= "Juncus effusus", 
          country = "Scotland") %>% 
  add_row(family = "Juncaceae", genus = "Luzula", stand.spec= "Luzula campestris", 
          country = "Scotland") %>% 
  add_row(family = "Primulaceae", genus = "Lysimachia", stand.spec= "Lysimachia europaea", 
          country = "Scotland") %>% 
  add_row(family = "Orobanchaceae", genus = "Melampyrum", stand.spec= "Melampyrum pratense", 
          country = "Scotland") %>% 
  add_row(family = "Poaceae", genus = "Molinia", stand.spec= "Poaceae", 
          country = "Scotland") %>% 
  add_row(family = "Poaceae", genus = NA, stand.spec= "Poaceae", 
          country = "Scotland") %>% 
  add_row(family = "Equisetaceae", genus = "Equisetum", stand.spec= "Equisetum sp.", 
          country = "Scotland") %>% 
  add_row(family = "Poaceae", genus = "Festuca", stand.spec= "Festuca sp.", 
          country = "Scotland") %>% 
  add_row(family = "Poaceae", genus = "Holcus", stand.spec= "Poaceae", 
          country = "Scotland") %>% 
  add_row(family = "Juncaceae", genus = "Juncus", stand.spec= "Juncus bufonius", 
          country = "Scotland") %>% 
  add_row(family = "Juncaceae", genus = "Juncus", stand.spec= "Juncus conglomeratus", 
          country = "Scotland") %>% 
  add_row(family = "Asteraceae", genus = "Scorzoneroides", stand.spec= "Scorzoneroides autumnalis", 
          country = "Scotland") %>% 
  add_row(family = "Juncaceae", genus = "Luzula", stand.spec= "Luzula multiflora", 
          country = "Scotland") %>% 
  add_row(family = "Juncaceae", genus = "Luzula", stand.spec= "Luzula pilosa", 
          country = "Scotland") %>% 
  add_row(family = "Nartheciaceae", genus = "Narthecium", stand.spec= "Narthecium ossifragum", 
          country = "Scotland") %>% 
  add_row(family = "Poaceae", genus = "Poa", stand.spec= "Poaceae", 
          country = "Scotland") %>% 
  add_row(family = "Polygalaceae", genus = "Polygala", stand.spec= "Polygala serpyllifolia", 
          country = "Scotland") %>% 
  add_row(family = "Salicaceae", genus = "Salix", stand.spec= "Salix sp.", 
          country = "Scotland") %>% 
  add_row(family = "Cyperaceae", genus = "Trichophorum", stand.spec= "Trichophorum cespitosum", 
          country = "Scotland") %>% 
  add_row(family = "Asteraceae", genus = NA, stand.spec= "Asteraceae", 
          country = "Scotland") %>% 
  add_row(family = "Orchidaceae", genus = "Dactylorhiza", stand.spec= "Dactylorhiza maculata", 
          country = "Scotland") %>% 
  add_row(family = "Poaceae", genus = "Deschampsia", stand.spec= "Poaceae", 
          country = "Scotland") %>% 
  add_row(family = "Fabaceae", genus =NA, stand.spec= "Fabaceae", 
          country = "Scotland") %>% 
  add_row(family = "Asteraceae", genus = "Hypochaeris", stand.spec= "Hypochaeris", 
          country = "Scotland") %>% 
  add_row(family = "Asteraceae", genus = "Jacobea", stand.spec= "Jacobea vulgaris", 
          country = "Scotland") %>% 
  add_row(family = "Oxalidaceae", genus = "Oxalis", stand.spec= "Oxalis acetosella", 
          country = "Scotland") %>% 
  add_row(family = "Orobanchaceae", genus = "Pedicularis", stand.spec= "Pedicularis sylvatica", 
          country = "Scotland") %>% 
  add_row(family = "Poaceae", genus = "Setaria", stand.spec= "Poaceae", 
          country = "Scotland") %>% 
  add_row(family = "Caprifoliaceae", genus = "Succisa", stand.spec= "Succisa pratensis", 
          country = "Scotland") %>% 
  add_row(family = "Poaceae", genus = "Alopecurus", stand.spec= "Alopecurus geniculatus", 
          country = "Scotland") %>% 
  add_row(family = "Boraginaceae", genus = "Myosotis", stand.spec= "Myosotis sp.", 
          country = "Scotland") %>% 
  add_row(family = "Ranunculaceae", genus = "Ranunculus", stand.spec= "Ranunculus flammula", 
          country = "Scotland") %>% 
  add_row(family = "Caryophyllaceae", genus = "Stellaria", stand.spec= "Stellaria holostea", 
          country = "Scotland") %>% 
  add_row(family = "Asteraceae", genus = "Taraxacum", stand.spec= "Taraxacum sp.", 
          country = "Scotland") %>% 
  add_row(family = "Blechnaceae", genus = "Blechnum", stand.spec= "Blechnum spicant", 
          country = "Scotland") %>% 
  add_row(family = "Asteraceae", genus = "Cirsium", stand.spec= "Cirsium sp.", 
          country = "Scotland") %>% 
  add_row(family = "Fabaceae", genus = "Lathyrus", stand.spec= "Lathyrus palustris", 
          country = "Scotland") %>% 
  add_row(family = "Orchidaceae", genus = "Neottia", stand.spec= "Neottia cordata", 
          country = "Scotland") %>% 
  add_row(family = "Poaceae", genus = "Glyceria", stand.spec= "Poaceae", 
          country = "Scotland") %>% 
  add_row(family = "Dryopteridaceae", genus = "Dryopteris", stand.spec= "Dryopteris", 
          country = "Scotland") %>% 
  add_row(family = "Primulaceae", genus = "Primula", stand.spec= "Primula vulgaris", 
          country = "Scotland") %>% 
  add_row(family = "Cyperaceae", genus = "Carex", stand.spec= "",
          country = "Switzerland") %>% 
  add_row(family = "Juncaceae", genus = "Juncus", stand.spec= "Juncus effusus",
          country = "Switzerland") %>% 
  add_row(family = "Juncaceae", genus = "Luzula", stand.spec= "Luzula campestris",
          country = "Switzerland") %>% 
  add_row(family = "Primulaceae", genus = "Lysimachia", stand.spec= "Lysimachia europaea",
          country = "Switzerland") %>% 
  add_row(family = "Orobanchaceae", genus = "Melampyrum", stand.spec= "Melampyrum pratense",
          country = "Switzerland") %>% 
  add_row(family = "Poaceae", genus = "Molinia", stand.spec= "Poaceae",
          country = "Switzerland") %>% 
  add_row(family = "Poaceae", genus = NA, stand.spec= "Poaceae",
          country = "Switzerland") %>% 
  add_row(family = "Equisetaceae", genus = "Equisetum", stand.spec= "Equisetum sp.",
          country = "Switzerland") %>% 
  add_row(family = "Poaceae", genus = "Festuca", stand.spec= "Festuca sp.",
          country = "Switzerland") %>% 
  add_row(family = "Poaceae", genus = "Holcus", stand.spec= "Poaceae",
          country = "Switzerland") %>% 
  add_row(family = "Juncaceae", genus = "Juncus", stand.spec= "Juncus bufonius",
          country = "Switzerland") %>% 
  add_row(family = "Juncaceae", genus = "Juncus", stand.spec= "Juncus conglomeratus",
          country = "Switzerland") %>% 
  add_row(family = "Asteraceae", genus = "Scorzoneroides", stand.spec= "Scorzoneroides autumnalis",
          country = "Switzerland") %>% 
  add_row(family = "Juncaceae", genus = "Luzula", stand.spec= "Luzula multiflora",
          country = "Switzerland") %>% 
  add_row(family = "Juncaceae", genus = "Luzula", stand.spec= "Luzula pilosa",
          country = "Switzerland") %>% 
  add_row(family = "Nartheciaceae", genus = "Narthecium", stand.spec= "Narthecium ossifragum",
          country = "Switzerland") %>% 
  add_row(family = "Poaceae", genus = "Poa", stand.spec= "Poaceae",
          country = "Switzerland") %>% 
  add_row(family = "Polygalaceae", genus = "Polygala", stand.spec= "Polygala serpyllifolia",
          country = "Switzerland") %>% 
  add_row(family = "Salicaceae", genus = "Salix", stand.spec= "Salix sp.",
          country = "Switzerland") %>% 
  add_row(family = "Cyperaceae", genus = "Trichophorum", stand.spec= "Trichophorum cespitosum",
          country = "Switzerland") %>% 
  add_row(family = "Asteraceae", genus = NA, stand.spec= "Asteraceae",
          country = "Switzerland") %>% 
  add_row(family = "Orchidaceae", genus = "Dactylorhiza", stand.spec= "Dactylorhiza maculata",
          country = "Switzerland") %>% 
  add_row(family = "Poaceae", genus = "Deschampsia", stand.spec= "Poaceae",
          country = "Switzerland") %>% 
  add_row(family = "Fabaceae", genus =NA, stand.spec= "Fabaceae",
          country = "Switzerland") %>% 
  add_row(family = "Asteraceae", genus = "Hypochaeris", stand.spec= "Hypochaeris",
          country = "Switzerland") %>% 
  add_row(family = "Asteraceae", genus = "Jacobea", stand.spec= "Jacobea vulgaris",
          country = "Switzerland") %>% 
  add_row(family = "Oxalidaceae", genus = "Oxalis", stand.spec= "Oxalis acetosella",
          country = "Switzerland") %>% 
  add_row(family = "Orobanchaceae", genus = "Pedicularis", stand.spec= "Pedicularis sylvatica",
          country = "Switzerland") %>% 
  add_row(family = "Poaceae", genus = "Setaria", stand.spec= "Poaceae",
          country = "Switzerland") %>% 
  add_row(family = "Caprifoliaceae", genus = "Succisa", stand.spec= "Succisa pratensis",
          country = "Switzerland") %>% 
  add_row(family = "Poaceae", genus = "Alopecurus", stand.spec= "Alopecurus geniculatus",
          country = "Switzerland") %>% 
  add_row(family = "Boraginaceae", genus = "Myosotis", stand.spec= "Myosotis sp.",
          country = "Switzerland") %>% 
  add_row(family = "Ranunculaceae", genus = "Ranunculus", stand.spec= "Ranunculus flammula",
          country = "Switzerland") %>% 
  add_row(family = "Caryophyllaceae", genus = "Stellaria", stand.spec= "Stellaria holostea",
          country = "Switzerland") %>% 
  add_row(family = "Asteraceae", genus = "Taraxacum", stand.spec= "Taraxacum sp.",
          country = "Switzerland") %>% 
  add_row(family = "Blechnaceae", genus = "Blechnum", stand.spec= "Blechnum spicant",
          country = "Switzerland") %>% 
  add_row(family = "Asteraceae", genus = "Cirsium", stand.spec= "Cirsium sp.",
          country = "Switzerland") %>% 
  add_row(family = "Fabaceae", genus = "Lathyrus", stand.spec= "Lathyrus palustris",
          country = "Switzerland") %>% 
  add_row(family = "Orchidaceae", genus = "Neottia", stand.spec= "Neottia cordata",
          country = "Switzerland") %>% 
  add_row(family = "Poaceae", genus = "Glyceria", stand.spec= "Poaceae",
          country = "Switzerland") %>% 
  add_row(family = "Dryopteridaceae", genus = "Dryopteris", stand.spec= "Dryopteris",
          country = "Switzerland") %>% 
  add_row(family = "Primulaceae", genus = "Primula", stand.spec= "Primula vulgaris",
          country = "Switzerland") 



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
            SLA = paste(round(mean(SLA), 1), "±", round(sd(SLA), 2)),
            LA = paste(round(mean(LA), 1), "±", round(sd(LA), 1)),
            `Plant height` = paste(round(mean(PlantHeight), 1), "±", round(sd(PlantHeight), 1)))

saveRDS(dfTRAIT, "RDS_files/02_Gapfilled_traits.rds")

rm(list=setdiff(ls(), "sum_trait"))