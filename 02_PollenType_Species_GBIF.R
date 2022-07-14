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
## - area of species selection
## - evaluate limit setting occ_search
##   
## ---------------------------


## Load packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(rgbif)) install.packages("rgbif")
if(!require(LCVP)) devtools::install_github("idiv-biodiversity/LCVP")
if(!require(lcvplants)) devtools::install_github("idiv-biodiversity/lcvplants")
if(!require(CoordinateCleaner)) install.packages("CoordinateCleaner")


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

## create function  that downloads, cleans, saves and creates conversion table for pollen type to species

gbifsaveclean <- function(taxonkey, countrykey){
  # perform gbif search
  occsearch <- occ_search(taxonKey = taxonkey,
                        country = countrykey,
                        limit = 1000,  # evaluate limit setting
                        year = "1990, 2020",
                        hasCoordinate = T)     
  
  # proceed to next steps of the function (cleaning) if the search retrieved data  
  if(!is.null(occsearch$data)) {
  # save raw data
  #saveRDS(occsearch, paste("RDS_files/RawGBIF/GBIFoccurances_", countrykey, taxonkey, ".rds", sep = ""))
  
  # only select data element, discard meta
  occdat <- occsearch$data 
  
  # clear empty elements (key didn't retrieve data)
  occdat <- occdat %>% compact()
  
  # clean occurrences
  occdat <- occdat %>% #purrr::map(. 
    data.frame() %>%
      cc_val(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records with empty coordinates
      cc_equ(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records with equal latitude and longitude coordinates
      cc_cap(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records within a certain radius around country capitals
      cc_cen(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records within a radius around the geographic centroids of political countries and provinces.
      cc_inst(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records assigned to the location of zoos, botanical gardens, herbaria, universities and museums
      cc_zero(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes records with zero  long or lat
      cc_outl(lon = "decimalLongitude", lat = "decimalLatitude") %>% # removes geographical outliers
      cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude") # removes duplicates
     
  # save cleaned data frame
  #saveRDS(occdat, paste("RDS_files/CleanGBIF/CleanedGBIFOccurances_", countrykey, taxonkey, ".rds", sep = ""))

  # create species per pollen taxon list
  spec <-
    occdat %>% 
    dplyr::select_if(colnames(.) %in% c("family", "genus", "species")) %>% 
    drop_na() %>% 
    arrange(species) 
  
  spec_num <-
    spec %>% 
    group_by(species) %>% 
    summarise(count = n()) %>% 
    drop_na() %>% 
    arrange(species)
  spec <- spec %>% distinct()
  spec <- cbind(spec, spec_num)
  spec <- spec[,-4] # remove duplicate species column
  
  # standardize species names according to the leipzig plant catalogue
  spec.stand <-
    lcvp_search(spec$species,progress_bar = TRUE) 
  spec <-
    spec %>% mutate(stand.spec = word(spec.stand$Output.Taxon, 1,2) %>% # only keep binomial without authority 
                      str_remove("_x"),# x added for unknown reason,
                    genus = word(spec.stand$Output.Taxon, 1),
                    family = spec.stand$Family,
                    country = countrykey)
 spec
  }
}

spec_gb <- 
  purrr::map(taxonkeys,
           ~gbifsaveclean(taxonkey = .x, countrykey = "GB")) %>% 
  # remove empty 
  keep(., ~is.data.frame(.)) %>% 
  bind_rows() 
spec_ch <- purrr::map(taxonkeys,
                   ~gbifsaveclean(taxonkey = .x, countrykey = "CH")) %>% 
  # remove empty 
  keep(., ~is.data.frame(.)) %>% 
  bind_rows() 

spec <- bind_rows("Scotland" = spec_gb, "Switzerland" = spec_ch, .id = "country")
saveRDS(spec, "RDS_files/02_gbif_raw.rds")

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
                              .$spec == "Rabelera holostea",
                              values = "Stellaria holostea")) %>% 
  # remove na's (after checking for synonyms manually, all cultivars)
  filter(!is.na(stand.spec))  
  
saveRDS(specpol, "RDS_files/02_PollenType_species.rds")
