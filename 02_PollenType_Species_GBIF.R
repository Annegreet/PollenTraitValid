# This script creates a conversion table for pollen type to species based on the current 
# distribution of terrestrial species belonging to the taxon in the UK
# Annegreet Veeken

# to do
# -replace if else by case_when

# libraries
library(rgbif)
library(tidyverse)
library(Taxonstand)
library(CoordinateCleaner)
library(readxl)

# Load in pollen data
dfPOL <- readRDS("RDS_files/01_Pollen_data_Scot.rds")

# extract taxon column for GBIF search
taxa <- dfPOL %>% 
  pull(pollentaxon) %>% 
  unique() 
taxa <- taxa[!taxa == "unindentified"]

# divide in genera, families and subfamily, to facilitate GBIF search
sp <- taxa %>% str_subset("[:space:]") 
sp <- sp[!sp == "Ericales (tetrad)"]
gn <- taxa %>% str_subset("^(?!.*eae)") %>% 
  str_subset("[:space:]", negate = TRUE)
gn <- gn[!gn == "Pteridophyte"]
fam <- taxa %>% str_subset("ceae")
ord <-  "Ericales"
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

ordlist <-
  purrr::map(ord, ~ name_suggest(., rank = 'ORDER', limit = 100))
ordkeys <- ordlist %>% 
  unlist(recursive = FALSE) %>% 
  keep(names(.) == "data") %>% 
  compact() %>% 
  purrr::map(., ~filter(., !is.na(canonicalName))) %>% 
  bind_rows() %>% 
  # only include keys for exact taxon matches, not partial matches (possibly evaluate limit setting in name_suggest)
  filter(canonicalName %in% ord) 

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
taxonkeys <-  bind_rows(spkeys, gnkeys, famkeys,ordkeys, classkeys) %>% 
  pull(key)
  
# get country keys
countries <- c("United Kingdom")
countrykeys <- isocodes %>% filter(name %in% countries) %>% pull(code)

## create function  that downloads, cleans, saves and creates conversion table for pollen type to species

gbifsaveclean <- function(taxonkey, countrykey){
  ## perform gbif search
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
  
  ## clean occurrences
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
    #purrr::map_dfr(. %>% 
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
  
  # standardise species names according to the Plant List
  spec.stand <-
    TPL(spec$species) 
  spec <-
    spec %>% mutate(stand.spec = str_c(spec.stand$New.Genus, 
                                       spec.stand$New.Species, sep = " "))
 spec
  }
}

spec <- purrr::map(taxonkeys,
           ~gbifsaveclean(taxonkey = .x, countrykey = "GB"))

# combine to make table, add pollen taxon

specpol <- spec %>% 
  # remove
  keep(., ~is.data.frame(.)) %>% 
  bind_rows() %>% 
  mutate(pollentaxon = ifelse(species %in% sp, species,
                              ifelse(genus %in% gn, genus, 
                                     ifelse(family %in% fam, family, 
                                            ifelse(family == "Ericaceae", "Ericales (tetrad)",
                                                   ifelse(family %in% c("Pteridaceae", "Aspleniaceae",
                                                                     "Athyriaceae", "Salviniaceae",
                                                                     "Dryopteridaceae", "Hymenophyllaceae",
                                                                     "Polypodiaceae", "Dryopteridaceae",
                                                                     "Dennstaedtiaceae", "Blechnaceae",
                                                                     "Thelypteridaceae"), "Pteridophyte", NA)))))) %>% 

  filter(!is.na(pollentaxon))

saveRDS(specpol, "RDS_files/02_PollenType_species.rds")
