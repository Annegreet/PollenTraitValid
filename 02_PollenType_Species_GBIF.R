## ---------------------------
##
## Script name: 02_PollenType_Species_GBIF
##
## Purpose of script: This script creates a conversion table for pollen type to 
## species based on the current distribution of terrestrial species belonging to 
## the taxon in the UK/Switzerland
##
##
## Author: Annegreet Veeken
##
## Date Created: 2022-06-30
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
## 
## 
##   
## ---------------------------


## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rgbif)) install.packages("rgbif")
if (!require(LCVP)) devtools::install_github("idiv-biodiversity/LCVP")
if (!require(lcvplants)) devtools::install_github("idiv-biodiversity/lcvplants")

# Load in pollen data
dfPOL_scot <- readRDS("RDS_files/01_Pollen_data_Scot.rds")
dfPOL_swiss <- readRDS("RDS_files/01_Pollen_data_Swiss.rds")
dfPOL <- bind_rows(dfPOL_scot, dfPOL_swiss)

# extract taxon column for GBIF search
taxa <- dfPOL %>% 
  pull(pollentaxon) %>% 
  unique() 
taxa <- taxa[!taxa == "unindentified"]
taxa <- c(taxa, "Oxyria", "Rumex")

# divide in genera, families and subfamily, to facilitate GBIF search
sp <- taxa %>% str_subset("[:space:]") 
sp <- sp[!sp == "Ericales (tetrad)"] # replaced by relevant tetrad ericales Ericaceae and empetraceae
gn <- taxa %>% str_subset("^(?!.*eae)") %>% 
  str_subset("[:space:]", negate = TRUE)
gn <- gn[!gn == "Pteridophyte"] # replaced by polypodiopsida
fam <- taxa %>% str_subset("ceae") %>% 
  c(., "Ericaceae", "Empetraceae" )
class <- "Polypodiopsida"

# find gbif id's
# for species
splist <-
  purrr::map(sp, ~ name_suggest(., rank = 'SPECIES', limit = 100))
spkeys <- splist %>% 
  unlist(recursive = FALSE) %>% 
  keep(names(.) == "data") %>% 
  compact() %>% 
  purrr::map(., ~filter(., !is.na(canonicalName))) %>% 
  bind_rows() %>% 
  # only include keys for exact taxon matches, not partial matches (possibly evaluate limit setting in name_suggest)
  filter(canonicalName %in% sp) 

# for genera
gnlist <-
  purrr::map(gn, ~ name_suggest(., rank = 'GENUS', limit = 100))
# check this list for missing returns

gnkeys <- gnlist %>% 
  unlist(recursive = FALSE) %>% 
  keep(names(.) == "data") %>% 
  compact() %>% 
  purrr::map(., ~filter(., !is.na(canonicalName))) %>% 
  bind_rows() %>% 
  # only include keys for exact taxon matches, not partial matches (possibly evaluate limit setting in name_suggest)
  filter(canonicalName %in% gn) 

# for family
famlist <-
  purrr::map(fam, ~ name_suggest(., rank = 'FAMILY', limit = 100))
famkeys <- famlist %>% 
  unlist(recursive = FALSE) %>% 
  keep(names(.) == "data") %>% 
  compact() %>% 
  purrr::map(., ~filter(., !is.na(canonicalName))) %>% 
  bind_rows() %>% 
  # only include keys for exact taxon matches, not partial matches (possibly evaluate limit setting in name_suggest)
  filter(canonicalName %in% fam) 

classlist <-
    purrr::map(class, ~ name_suggest(.,limit = 100))
classkeys <- classlist %>% 
  unlist(recursive = FALSE) %>% 
  keep(names(.) == "data") %>% 
  compact() %>% 
  purrr::map(., ~filter(., !is.na(canonicalName))) %>% 
  bind_rows() %>% 
  # only include keys for exact taxon matches, not partial matches (possibly evaluate limit setting in name_suggest)
  filter(canonicalName %in% class) 

# concatenate all taxon keys
taxonkeys <-  bind_rows(spkeys, gnkeys, famkeys, classkeys) %>% 
  pull(key)
  
# get country keys
countries <- c("United Kingdom", "Switzerland")
countrykeys <- isocodes %>% filter(name %in% countries) %>% pull(code)

countrykey <- countrykeys[1]
taxonkey <- taxonkeys[1:2]

## create function  that downloads, cleans, saves and creates conversion table for pollen type to species
if (1) {
gbifsaveclean <- function(taxonkey, countrykey){
  # perform gbif download request
  occsearch <- occ_download(
    pred_in("basisOfRecord", "HUMAN_OBSERVATION"),
    pred_in("taxonKey", taxonkey),
    pred_in("country", countrykey),
    pred_gte("year", 2012), # observation from 2012 onwards
    pred("hasCoordinate", TRUE),
    format = "SPECIES_LIST") # only download species list, not occurrence points
  
  occdata <- occ_download_wait(occsearch)
  
  occdata <- occ_download_get(occsearch, path = "Data/GBIF-requests/",
                              overwrite = TRUE) %>% 
    occ_download_import()
  # for already downloaded data
  occdata <- occ_download_get(str_remove(list.files(path = "Data/GBIF-requests/")[1], ".zip"),
                              overwrite = TRUE) %>% 
    occ_download_import()
  
  # create species per pollen taxon list
  spec <-
    occdata %>% 
    filter(taxonRank == "SPECIES") %>% 
    dplyr::select_if(colnames(.) %in% c("family", "genus", "species")) %>% 
    drop_na() %>% 
    arrange(species) 
  
  spec <- spec %>% distinct()
  
  # standardize species names according to the leipzig plant catalogue
  spec.stand <-
    lcvp_search(spec$species, progress_bar = TRUE) 
  
  spec <-
    spec %>% mutate(stand.spec = word(spec.stand$Output.Taxon, 1,2) %>% # only keep binomial without authority 
                      str_remove("_x"),# x added for unknown reason,
                    genus = word(spec.stand$Output.Taxon, 1),
                    family = spec.stand$Family,
                    country = countrykey)
 spec
 }

spec_gb <- gbifsaveclean(taxonkey = taxonkeys, countrykey = "GB")

spec_ch <- gbifsaveclean(taxonkey = taxonkeys, countrykey = "CH")

spec <- bind_rows("Scotland" = spec_gb, "Switzerland" = spec_ch, .id = "country")
saveRDS(spec, "RDS_files/02_gbif_raw.rds")
}

# combine to make table, add pollen taxon
specpol <- spec %>% 
  mutate(pollentaxon = case_when(species %in% sp ~ species,
                                 genus %in% c("Rumex", "Oxyria") ~ "Rumex/Oxyria",
                                 genus %in% gn ~ genus,
                                 family == "Ericaceae" ~ "Ericales (tetrad)",
                                 family %in% c("Pteridaceae", "Aspleniaceae",
                                               "Athyriaceae", "Salviniaceae",
                                               "Dryopteridaceae", "Hymenophyllaceae",
                                               "Polypodiaceae", "Dryopteridaceae",
                                               "Dennstaedtiaceae", "Blechnaceae",
                                               "Thelypteridaceae", "Equisetaceae",
                                               "Osmundaceae","Cystopteridaceae") ~ "Pteridophyte",
                                 family %in% fam ~ family
                                 )
         )  %>% 
  distinct(stand.spec, pollentaxon, country, .keep_all = TRUE)

# check for duplicate values (species assigned to two or more pollen type) 
duplicats <- specpol %>%
  filter(duplicated(stand.spec),
         duplicated(stand.spec, fromLast = TRUE))

specpol <- specpol %>% 
  # missing alternative name
  mutate(stand.spec = replace(.$stand.spec,
                              .$species == "Rabelera holostea",
                              values = "Stellaria holostea")) %>% 
  # lcvp didn't retrieve synonym
  filter(!species == "Aria edulis") %>% 
  add_row(family = "Rosaceae", genus = "Sorbus", species = "Sorbus aria", 
          stand.spec = "Sorbus aria", country = "Switzerland", pollentaxon = "Rosaceae") %>% 
  # remove na's (after checking for synonyms manually, all cultivars)
  filter(!is.na(stand.spec)) %>% 
  as_tibble()

# some species recorded in the field not retrieved by gbif 
# (mostly because pollen not represented in pollen record)
dfABUN_a <- readRDS("RDS_files/01_Species_abundance_a.rds")
dfABUN_bc <- readRDS("RDS_files/01_Species_abundance_bc.rds") %>% 
  ungroup
bdm_abun <- readRDS("RDS_files/01_Species_abundance_Swiss.rds") 

# Scotland
sp_extra <- c(dfABUN_a$stand.spec, dfABUN_bc$stand.spec) %>% 
  unique %>% 
  .[!. %in% unique(specpol$stand.spec)] # check for missing species
# find synonyms
lcvp_out <- lcvp_search(str_subset(sp_extra, " ")) %>% as_tibble()
# create table with pollentaxon
lcvp_scot <- lcvp_out %>% 
  mutate(genus = word(Output.Taxon, 1),
         stand.spec = word(Output.Taxon, 1,2)) %>% 
  dplyr::select(family = Family, genus , species = Search, stand.spec) %>% 
  drop_na() %>% 
  mutate(pollentaxon = case_when(species %in% sp ~ species,
                                 genus %in% c("Rumex", "Oxyria") ~ "Rumex/Oxyria",
                                 genus %in% gn ~ genus,
                                 family == "Ericaceae" ~ "Ericales (tetrad)",
                                 family %in% c("Pteridaceae", "Aspleniaceae",
                                               "Athyriaceae", "Salviniaceae",
                                               "Dryopteridaceae", "Hymenophyllaceae",
                                               "Polypodiaceae", "Dryopteridaceae",
                                               "Dennstaedtiaceae", "Blechnaceae",
                                               "Thelypteridaceae", "Equisetaceae",
                                               "Osmundaceae","Cystopteridaceae") ~ "Pteridophyte",
                                 family %in% fam ~ family)
  ) 
# Switzerland
sp_extra <- bdm_abun$stand.spec %>% 
  unique %>% 
  .[!. %in% unique(specpol$stand.spec)] # check for missing species
# find synonyms
lcvp_out <- lcvp_search(str_subset(sp_extra, " ")) %>% as_tibble()
# create table with pollentaxon
lcvp_swiss <- lcvp_out %>% 
  mutate(genus = word(Output.Taxon, 1),
         stand.spec = word(Output.Taxon, 1,2)) %>% 
  dplyr::select(family = Family, genus , species = Search, stand.spec) %>% 
  drop_na() %>% 
  mutate(pollentaxon = case_when(species %in% sp ~ species,
                                 genus %in% c("Rumex", "Oxyria") ~ "Rumex/Oxyria",
                                 genus %in% gn ~ genus,
                                 family == "Ericaceae" ~ "Ericales (tetrad)",
                                 family %in% c("Pteridaceae", "Aspleniaceae",
                                               "Athyriaceae", "Salviniaceae",
                                               "Dryopteridaceae", "Hymenophyllaceae",
                                               "Polypodiaceae", "Dryopteridaceae",
                                               "Dennstaedtiaceae", "Blechnaceae",
                                               "Thelypteridaceae", "Equisetaceae",
                                               "Osmundaceae","Cystopteridaceae") ~ "Pteridophyte",
                                 family %in% fam ~ family)
  ) 

# add to pollen species table
specpol <- specpol  %>% 
  add_row(lcvp_scot, country = "Scotland") %>%
  add_row(lcvp_swiss, country = "Switzerland") %>% 
  # missing species
  add_row(family = "Cyperaceae", genus = "Carex", stand.spec = "Carex", 
             country = "Switzerland", pollentaxon = "Cyperaceae") %>% 
  add_row(family = "Juncaceae", genus = "Juncus", stand.spec = "Juncus", 
          country = "Switzerland", pollentaxon = "Cyperaceae") %>%
  add_row(family = "Caryophyllaceae", genus = "Cerastium", stand.spec = "Cerastium", 
          country = "Switzerland", pollentaxon = "Caryophyllaceae") %>% 
  add_row(family = "Ranunculaceae", genus = "Ranunculus", stand.spec = "Ranunculus", 
          country = "Switzerland", pollentaxon = "Ranunculaceae") %>% 
  add_row(family = "Poaceae", genus = "Festuca", stand.spec = "Festuca", 
          country = "Switzerland", pollentaxon = "Poaceae") %>% 
  add_row(family = "Poaceae", genus = NA, stand.spec = "Poaceae",
          country = "Switzerland", pollentaxon = "Poaceae") %>% 
  add_row(family = "Equisetaceae", genus = "Equisetum", stand.spec = "Equisetum",
          country = "Switzerland", pollentaxon = "Pteridophyte") %>% 
  add_row(family = "Salicaceae", genus = "Salix", stand.spec = "Salix",
          country = "Switzerland", pollentaxon = "Salix") %>% 
  add_row(family = "Asteraceae", genus = NA, stand.spec = "Asteraceae",
          country = "Switzerland", pollentaxon = "Asteraceae") %>% 
  add_row(family = "Fabaceae", genus = NA, stand.spec = "Fabaceae",
          country = "Switzerland", pollentaxon = "Fabaceae") %>% 
  add_row(family = "Asteraceae", genus = "Hypochaeris", stand.spec = "Hypochaeris",
          country = "Switzerland", pollentaxon = "Asteraceae") %>% 
  add_row(family = "Boraginaceae", genus = "Myosotis", stand.spec = "Myosotis",
          country = "Switzerland", pollentaxon = NA) %>% 
  add_row(family = "Asteraceae", genus = "Taraxacum", stand.spec = "Taraxacum",
          country = "Switzerland", pollentaxon = "Asteraceae") %>% 
  add_row(family = "Asteraceae", genus = "Cirsium", stand.spec = "Cirsium",
          country = "Switzerland", pollentaxon = "Asteraceae") %>% 
  add_row(family = "Poaceae", genus = "Glyceria", stand.spec = "Poaceae",
          country = "Switzerland", pollentaxon = "Poaceae") %>% 
  add_row(family = "Dryopteridaceae", genus = "Dryopteris", stand.spec = "Dryopteris",
          country = "Switzerland", pollentaxon = "Pteridophyte") %>% 
  add_row(family = "Blechnaceae", genus = "Blechnum", stand.spec = "Blechnum spicant",
          country = "Switzerland", pollentaxon = "Pteridophyte") %>% 
  add_row(family = "Asteraceae", genus = "Leontondon", stand.spec = "Leontondon",
          country = "Switzerland", pollentaxon = "Asteraceae") %>% 
  add_row(family = "Betulaceae", genus = "Betula", stand.spec = "Betula",
          country = "Switzerland", pollentaxon = "Betula") %>% 
  add_row(family = "Cupressaceae", genus = "Cupressus", stand.spec = "Cupressus",
          country = "Switzerland", pollentaxon = NA) %>% 
  add_row(family = "Nothofagaceae", genus = "Nothofagus", stand.spec = "Nothofagus",
          country = "Switzerland", pollentaxon = NA) %>% 
  add_row(family = "Pinaceae", genus = "Pinus", stand.spec = "Pinus",
          country = "Switzerland", pollentaxon = "Pinus") %>% 
  add_row(family = "Salicaceae", genus = "Populus", stand.spec = "Populus",
          country = "Switzerland", pollentaxon = "Populus") %>% 
  add_row(family = "Fagaceae", genus = "Quercus", stand.spec = "Quercus",
          country = "Switzerland", pollentaxon = "Quercus") %>% 
  add_row(family = "Malvaceae", genus = "Tilia", stand.spec = "Tilia",
          country = "Switzerland", pollentaxon = "Tilia") %>% 
  ##
  add_row(family = "Cyperaceae", genus = "Carex", stand.spec = "Carex", 
        country = "Scotland", pollentaxon = "Cyperaceae") %>% 
  add_row(family = "Juncaceae", genus = "Juncus", stand.spec = "Juncus", 
          country = "Scotland", pollentaxon = "Cyperaceae") %>%
  add_row(family = "Caryophyllaceae", genus = "Cerastium", stand.spec = "Cerastium", 
          country = "Scotland", pollentaxon = "Caryophyllaceae") %>% 
  add_row(family = "Ranunculaceae", genus = "Ranunculus", stand.spec = "Ranunculus", 
          country = "Scotland", pollentaxon = "Ranunculaceae") %>% 
  add_row(family = "Poaceae", genus = "Festuca", stand.spec = "Festuca", 
          country = "Scotland", pollentaxon = "Poaceae") %>% 
  add_row(family = "Poaceae", genus = NA, stand.spec = "Poaceae",
          country = "Scotland", pollentaxon = "Poaceae") %>% 
  add_row(family = "Equisetaceae", genus = "Equisetum", stand.spec = "Equisetum",
          country = "Scotland", pollentaxon = "Pteridophyte") %>% 
  add_row(family = "Salicaceae", genus = "Salix", stand.spec = "Salix",
          country = "Scotland", pollentaxon = "Salix") %>% 
  add_row(family = "Asteraceae", genus = NA, stand.spec = "Asteraceae",
          country = "Scotland", pollentaxon = "Asteraceae") %>% 
  add_row(family = "Fabaceae", genus = NA, stand.spec = "Fabaceae",
          country = "Scotland", pollentaxon = "Fabaceae") %>% 
  add_row(family = "Asteraceae", genus = "Hypochaeris", stand.spec = "Hypochaeris",
          country = "Scotland", pollentaxon = "Asteraceae") %>% 
  add_row(family = "Boraginaceae", genus = "Myosotis", stand.spec = "Myosotis",
          country = "Scotland", pollentaxon = NA) %>% 
  add_row(family = "Asteraceae", genus = "Taraxacum", stand.spec = "Taraxacum",
          country = "Scotland", pollentaxon = "Asteraceae") %>% 
  add_row(family = "Asteraceae", genus = "Cirsium", stand.spec = "Cirsium",
          country = "Scotland", pollentaxon = "Asteraceae") %>% 
  add_row(family = "Poaceae", genus = "Glyceria", stand.spec = "Poaceae",
          country = "Scotland", pollentaxon = "Poaceae") %>% 
  add_row(family = "Dryopteridaceae", genus = "Dryopteris", stand.spec = "Dryopteris",
          country = "Scotland", pollentaxon = "Pteridophyte") %>% 
  add_row(family = "Blechnaceae", genus = "Blechnum", stand.spec = "Blechnum spicant",
          country = "Scotland", pollentaxon = "Pteridophyte") %>% 
  add_row(family = "Asteraceae", genus = "Leontondon", stand.spec = "Leontondon",
          country = "Scotland", pollentaxon = "Asteraceae") %>%
  add_row(family = "Betulaceae", genus = "Betula", stand.spec = "Betula",
          country = "Scotland", pollentaxon = "Betula") %>% 
  add_row(family = "Cupressaceae", genus = "Cupressus", stand.spec = "Cupressus",
          country = "Scotland", pollentaxon = NA) %>% 
  add_row(family = "Nothofagaceae", genus = "Nothofagus", stand.spec = "Nothofagus",
          country = "Scotland", pollentaxon = NA) %>% 
  add_row(family = "Pinaceae", genus = "Pinus", stand.spec = "Pinus",
          country = "Scotland", pollentaxon = "Pinus") %>% 
  add_row(family = "Salicaceae", genus = "Populus", stand.spec = "Populus",
          country = "Scotland", pollentaxon = "Populus") %>% 
  add_row(family = "Fagaceae", genus = "Quercus", stand.spec = "Quercus",
          country = "Scotland", pollentaxon = "Quercus") %>% 
  add_row(family = "Malvaceae", genus = "Tilia", stand.spec = "Tilia",
          country = "Scotland", pollentaxon = "Tilia") %>% 
  distinct()


saveRDS(specpol, "RDS_files/02_PollenType_species.rds")
