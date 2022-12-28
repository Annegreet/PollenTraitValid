## ---------------------------
##
## Script name: Swiss data compilation
##
## Purpose of script: Compile swiss trait and species information
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

# Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(LCVP)) install.packages("LCVP")
if (!require(lcvplants)) install.packages("lcvplants")

# Load data
bdm_meta <- read_csv("Data/bdm_metadata.csv") # contains plot codes, X/Y coordinates (on the Swiss grid), elevation, some basic climate data (MAT, TAP) and a habitat classification
bdm_sp <- read_csv("Data/bdm_spp_comp.csv") # contains the relative abundances of all species (columns) in each plot/subplot (rows)
bdm_sp_codes <- read_csv("Data/bdm_spp_codes.csv") # contains the species codes for species IDs in the traits and species composition datasets
bdm_traits <- read_csv("Data/bdm_traits.csv") #  contains the dry mass, leaf area, SLA and height of plants in each plot/subplot 
bdm_pollen <- readRDS("RDS_files/01_Pollen_data_Swiss.rds")

# extract plot ids 
plot_id <- bdm_pollen %>% 
  pull(sitename) %>% 
  str_remove(pattern = "X") %>% 
  unique()

# check sites 
bdm_meta <- bdm_meta %>% 
  filter(aID_STAO %in% plot_id) 

# saveRDS(bdm_meta, "RDS_files/01_Meta_data_BDM.rds")
# bdm_meta <- bdm_meta %>%
#   filter(aID_STAO %in% bdm_sp$plot)
# write.csv(bdm_meta, file = "Data/Swiss_LC/all_BDM_plots_coord.csv")

# Prepare sp code for merge with trait and abundance data
bdm_sp_codes <- bdm_sp_codes %>% 
  # create 1 species column
  mutate(binom = case_when(!is.na(Species) ~ paste(Genus,Species),
                           !is.na(Genus) ~ paste(Genus, "sp."),
                           !is.na(Family) ~ paste(Family))
         ) %>% 
  # correct spelling mistakes
  mutate(binom = recode(binom, 
                        `Trifolum reppens` = "Trifolium repens",
                        `Cragaetus laevigata` = "Crataegus laevigata",
                        `Helix hedera` = "Hedera helix",
                        `Leontondon helveticus` = "Leontodon helveticus",
                        `Leotondon Autumnalis` = "Leontodon autumnalis",
                        `Cirisium spinosissimum` = "Cirsium spinosissimum",
                        `vaccinium myrtillus` = "Vaccinium myrtillus",
                        `vaccinium vitis idaea` = "Vaccinium vitis-idaea",
                        `Anthoxanthum odoratum aggr.` = "Anthoxanthum odoratum"),
         Species_code = recode(Species_code, 
                               soldal = "solalp", # Soldanella alpinia idem
                               soraur = "sorauc", # Sorbus aucuparia
                               lonneg = "lonnig", # Lonicera nigra
                               giumon = "geumon" # Geum montanum
                               ) 
         ) %>% 
  # missing codes
  add_row(Species_code = c("Rubfru", "Alcxan"),
          Species = c("Rubus fruticosus", "Alchemilla xanthochlora"))

## Abundance ----
# Join species codes with abundance data
bdm_sp <- bdm_sp %>% 
  # convert to long
  pivot_longer(!plot:subplot, names_to = "species_code", values_to = "abun") %>% 
  # remove 0 occurrences
  filter(!abun == 0) %>% 
  # join species codes
  left_join(bdm_sp_codes, by = c("species_code" = "Species_code")) %>% 
  dplyr::select(plot, subplot, abun, binom, Species, species_code)

# standardize species names
spec <- bdm_sp %>% 
  filter(!is.na(Species)) %>% 
  pull(binom) %>% 
  unique()

lcvp <- lcvp_search(spec) # find synonyms

synonyms <- lcvp %>% 
  dplyr::select(binom = Search, stand.spec = Output.Taxon, family = Family) %>% 
  mutate(stand.spec = word(stand.spec, 1, 2),# drop authorship
         genus = word(stand.spec,1)) 
regex_pattern <-
  setNames(synonyms$stand.spec,
           paste0("\\b", synonyms$binom, "\\b"))

bdm_sp <- bdm_sp %>% 
  # replace by standardized name
  mutate(stand.spec = str_replace_all(binom, regex_pattern)) %>% 
  # remove sp. epiphet
  mutate(stand.spec = str_remove(stand.spec, pattern = " sp\\.")) %>% 
  # code unidentified and no plant categories
  mutate(stand.spec = case_when(!is.na(stand.spec) ~ stand.spec,
                           species_code %in% c("FO","GR", "CR") ~ "Unindentified",
                           species_code %in% c("litter", "moss", "deadwood",
                                                   "dead wood", "dead trunk",
                                                   "stump", "fungi", "rock",
                                                   "root", "other", "bareground") ~ "No plants",
                           (str_detect(species_code, pattern = "[:digit:]") & !is.na(binom)) ~ binom,
                           (str_detect(species_code, pattern = "[:digit:]") & is.na(binom)) ~ "Unindentified"
                                )
         ) %>% 
  # select relevant columns, rename
  dplyr::select(sitename = plot, subplot,stand.spec, abun)

# calculate plant abundance
bdm_sp <- bdm_sp %>% 
  # filter stuff that is not a plant
  filter(!stand.spec == "No plants") %>%
  # calculate species abundance on the plot level
  group_by(sitename,stand.spec) %>% 
  summarise(abun = sum(abun)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  # add family and genus info again
  left_join(synonyms, by = c("stand.spec")) %>% 
  select(-binom) %>% 
  mutate(genus = case_when(is.na(genus) & !stand.spec == "Unindentified" & !str_ends(stand.spec, pattern = "ceae") &
                           str_detect(stand.spec, pattern = " ", negate = TRUE) ~ stand.spec,
                           TRUE ~ genus),
         family = case_when(is.na(family) & str_ends(stand.spec, pattern = "ceae") ~ stand.spec,
                            is.na(family) & genus == "Carex" ~ "Cyperaceae",
                            is.na(family) & genus == "Juncus" ~ "Juncaceae",
                            is.na(family) & genus == "Leontondon" ~ "Asteraceae",
                            is.na(family) & genus == "Dryopteris" ~ "Dryopteridaceae",
                            is.na(family) & genus == "Cerastium" ~ "Caryophyllaceae",
                            is.na(family) & genus == "Ranunculus" ~ "Ranunculaceae",
                            TRUE ~ family))

saveRDS(bdm_sp, "RDS_files/01_Species_abundance_Swiss.rds")

## Traits ----
bdm_traits <- bdm_traits %>% 
  left_join(bdm_sp_codes, by = c("species_code" = "Species_code"))
bdm_traits <- bdm_traits %>% 
  # replace by standardized name
  mutate(stand.spec = str_replace_all(binom, regex_pattern)) %>% 
  # code unidentified and no plant categories
  mutate(stand.spec = case_when(!is.na(stand.spec) ~ stand.spec,
                                species_code %in% c("FO","GR", "CR") ~ "Unindentified",
                                species_code %in% c("litter", "moss", "deadwood",
                                                    "dead wood", "dead trunk",
                                                    "stump", "fungi", "rock",
                                                    "root", "other", "bareground") ~ "No plants",
                                str_detect(species_code, pattern = "[:digit:]") ~ "Unindentified") # missing species code
  ) %>% 
  # select relevant columns, rename to common name 
  dplyr::select(species_code, family = Family, genus = Genus, 
                stand.spec, PFT = Functional_group,  
                LA = area_mg, SLA = SLA_mm2_mg, PlantHeight = height_cm) %>% 
  mutate(LA = LA * 100 #cm2 to mm2
         )
saveRDS(bdm_traits, "RDS_files/01_Traits_Swiss.rds")


  