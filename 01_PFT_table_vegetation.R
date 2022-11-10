## ---------------------------
##
## Script name: PFT_table_vegetation
##
## Purpose of script: Create pollination mode and pft table for pollen data
##
## Author: Annegreet Veeken
##
## Date Created: 2022-07-14
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   -
##
## ---------------------------

# Libraries
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(data.table)) install.packages("data.table")

# Pollination and PFT for vegetation----
traits <- fread(
  "Data/TRY-request/22386_31082022154852/22386.txt",
  sep = "\t",
  data.table = FALSE,
  stringsAsFactors = FALSE,
  strip.white = TRUE,
  drop = c("LastName", "FirstName","Comment", 
           "V28", "Reference", "ValueKindName", "Replicates"))

# get species list of vegetation
dfABUN_a <- readRDS("RDS_files/01_Species_abundance_a.rds")
dfABUN_bc <- readRDS("RDS_files/01_Species_abundance_bc.rds") %>% 
  ungroup
bdm_abun <- readRDS("RDS_files/01_Species_abundance_Swiss.rds") 

spec <- c(dfABUN_a$stand.spec, dfABUN_bc$stand.spec, bdm_abun$stand.spec) %>% 
  unique 

# create table with plant functional type data
pft <- traits %>% 
  filter(Dataset %in% c("Reich-Oleksyn Global Leaf N, P Database",
                        " Plant growth form dataset for the New World",
                        "Categorical Plant Traits Database",
                        "PLANTSdata USDA")) %>%
  filter(TraitName == "Plant growth form") %>% 
  # fix duplicates and mistakes
  mutate(AccSpeciesName = recode(AccSpeciesName, 
                                 `VACCINIUM VITIS-IDAEA` = "Vaccinium vitis-idaea")) %>% 
  group_by(AccSpeciesName) %>% 
  # generalize pft values
  summarise(pft = paste(unique(OrigValueStr), collapse = ", ")) %>% 
  mutate(growthform = case_when(str_detect(pft, "tree|Tree") ~ 'tree',
                                str_detect(pft, "herb|Herb|Forb") ~ 'herb',
                                str_detect(pft, "shrub|Shrub") ~ 'shrub',
                                str_detect(pft, "grass|Grass|Graminoid|graminoid") ~ 'grass',
                                str_detect(pft, "fern|Fern") ~ 'fern')) %>% 
  filter(!is.na(growthform)) %>% 
  distinct(AccSpeciesName, growthform) %>% 
  # fix mistakes
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
  filter(!(AccSpeciesName == "Juncus effusus" & growthform == "tree")) %>% 
  # add missing
  add_row(AccSpeciesName = c("Holcus lanatus","Molinia caerulea",
                             "Carex","Festuca",
                             "Jacobaea vulgaris","Trichophorum caespitosum",
                             "Trientalis europaea","Salix",
                             "Taraxacum","Primula acaulis",
                             "Equisetum","Hypochaeris",
                             "Myosotis","Cirsium",
                             "Equisetum sylvaticum","Glyceria fluitans",
                             "Lonicera nigra","Unindentified",
                             "Hedera helix","Pilosella cymosa",
                             "Anemone hepatica","Vicia sepium",
                             "Ranunculus lanuginosus","Alchemilla vulgaris",
                             "Leontondon","Mutellina purpurea",
                             "Cerastium","Juncus",
                             "Ranunculus tuberosus","Soldanella alpina",
                             "Clematis vitalba","Bistorta officinalis",
                             "Vaccinium nitens","Euphrasia minima",
                             "Ranunculus","Polygaloides chamaebuxus",
                             "Agrostis capillaris", "Anthoxanthum odoratum",
                             "Arrhenatherum elatius", "Athyrium filix-femina",
                             "Corylus avellana","Crataegus laevigata", 
                             "Dactylis glomerata", "Deschampsia flexuosa", 
                             "Juniperus communis", "Lolium perenne","Luzula sylvatica",
                             "Orthilia secunda","Pteridium aquilinum", "Rubus caesius",
                             "Juncus bufonius", "Larix decidua", "Lathyrus palustris",
                             "Galium odoratum","Galium verum", "Geranium robertianum",
                             "Clematis vitalba", "Helianthemum salicifolium",
                             "Rhinanthus glacialis", "Potentilla brevifolia",
                             "Blechnum spicant","Sambucus nigra" ,"Abies grandis",  "Abies procera", 
                             "Acer campestre"  , "Aesculus hippocastanum", "Alnus glutinosa","Alnus incana",
                             "Betula pubescens", "Betula", "Botrychium lunaria",  "Carpinus betulus",
                             "Castanea sativa", "Chamaecyparis lawsoniana", "Crataegus monogyna", "Cryptomeria japonica",
                             "Cupressus", "Ilex aquifolium", "Larix kaempferi","Malus sylvestris",
                             "Nothofagus", "Picea omorika", "Picea sitchensis", "Pinus contorta",
                             "Pinus nigra", "Pinus", "Populus", "Populus tremula",
                             "Prunus avium"   , "Prunus padus"  , "Pseudotsuga menziesii",
                             "Quercus robur","Quercus rubra", "Quercus", "Salix caprea",
                             "Salix cinerea","Sequoia sempervirens", "Sorbus aria", "Sorbus torminalis",
                             "Taxus baccata","Thuja plicata", "Tilia", "Tsuga heterophylla",
                             "Ulmus glabra","Ulmus procera"), 
          growthform = c("grass","grass", "grass", "grass", "herb",
                         "grass","herb", "shrub", "herb", "herb",
                         "fern", "herb", "herb", "herb", "fern",
                         "grass", "shrub", NA, "shrub", "herb",
                         "herb","herb","herb","herb","herb",
                         "herb","herb","grass","herb","herb",
                         "herb","herb","shrub","herb","herb",
                         "herb", "grass", "grass", "grass",
                         "fern", "shrub", "shrub", "grass",
                         "grass", "shrub","grass","grass", "herb",
                         "fern", "shrub", "grass", "tree", "herb",
                         "herb","herb","herb","herb", "herb", "herb",
                         "herb",
                         "fern", "shrub", "tree", "tree", 
                         "tree", "tree", "fern", "tree",
                         "tree", "tree", "shrub", "tree",
                         "tree","tree","tree","tree",
                         "tree","tree","tree","tree",
                         "tree","tree","tree","tree",
                         "tree","tree","tree",
                         "tree","tree","tree","tree",
                         "tree","tree","tree","tree",
                         "tree","tree","tree","tree",
                         "tree","tree","tree","tree",
                         "tree","tree"
                         )
  ) %>% 
  # Ericaceae specified as herb in pollen data -> do the same with the vegetation
  mutate(growthform = case_when(AccSpeciesName %in% c("Calluna vulgaris", 
                                                      "Erica cinerea",
                                                      "Erica tetralix",
                                                      "Vaccinium myrtillus",
                                                      "Empetrum nigrum",
                                                      "Vaccinium vitis-idaea",
                                                      "Vaccinium nitens") ~ "herb",
                                TRUE ~ growthform))

polmode <- traits %>% 
  filter(TraitName == "Pollination syndrome") %>% 
  mutate(AccSpeciesName = recode(AccSpeciesName, 
                                 `VACCINIUM VITIS-IDAEA` = "Vaccinium vitis-idaea")) %>% 
  group_by(AccSpeciesName) %>% 
  summarise(polination = paste(unique(OrigValueStr), collapse = ", ")) %>% 
  mutate(polmode = case_when(str_detect(polination, "wind") ~ "wind",
                             !str_detect(polination, "wind") ~ "not wind")) %>% 
  distinct(AccSpeciesName, polmode) %>% 
  # add missing
  add_row(AccSpeciesName = c("Holcus lanatus","Molinia caerulea",
                             "Carex","Festuca",
                             "Jacobaea vulgaris","Trichophorum caespitosum",
                             "Trientalis europaea","Salix",
                             "Taraxacum","Primula acaulis",
                             "Equisetum","Hypochaeris",
                             "Myosotis","Cirsium",
                             "Equisetum sylvaticum","Glyceria fluitans",
                             "Lonicera nigra","Unindentified",
                             "Hedera helix","Pilosella cymosa",
                             "Anemone hepatica","Vicia sepium",
                             "Ranunculus lanuginosus","Alchemilla vulgaris",
                             "Leontondon","Mutellina purpurea",
                             "Cerastium","Juncus",
                             "Ranunculus tuberosus","Soldanella alpina",
                             "Clematis vitalba","Bistorta officinalis",
                             "Vaccinium nitens","Euphrasia minima",
                             "Ranunculus","Polygaloides chamaebuxus",
                             "Asteraceae", "Dryopteris","Poaceae",
                             "Selaginella selaginoides",
                             "Juncus bufonius", "Larix decidua", "Lathyrus palustris",
                             "Potentilla brevifolia","Carex curvula", "Lycopodium annotinum",
                             "Trifolium alpinum",
                             "Blechnum spicant","Sambucus nigra" ,"Abies grandis",  "Abies procera", 
                             "Acer campestre"  , "Aesculus hippocastanum", "Alnus glutinosa","Alnus incana",
                             "Betula pubescens", "Betula", "Botrychium lunaria",  "Carpinus betulus",
                             "Castanea sativa", "Chamaecyparis lawsoniana", "Crataegus monogyna", "Cryptomeria japonica",
                             "Cupressus", "Ilex aquifolium", "Larix kaempferi","Malus sylvestris",
                             "Nothofagus", "Picea omorika", "Picea sitchensis", "Pinus contorta",
                             "Pinus nigra", "Pinus", "Populus", "Populus tremula",
                             "Prunus avium"   , "Prunus padus"  , "Pseudotsuga menziesii",
                             "Quercus robur","Quercus rubra", "Quercus", "Salix caprea",
                             "Salix cinerea","Sequoia sempervirens", "Sorbus aria", "Sorbus torminalis",
                             "Taxus baccata","Thuja plicata", "Tilia", "Tsuga heterophylla",
                             "Ulmus glabra","Ulmus procera"),
          polmode = c("wind","wind", "wind", "wind", "not wind",
                      "wind","not wind", "wind", "not wind", "not wind",
                      "wind", "not wind", "not wind", "not wind", "wind",
                      "wind", "not wind", NA, "not wind", "not wind",
                      "not wind","not wind","not wind","not wind","not wind",
                      "not wind","not wind","wind","not wind","not wind",
                      "not wind","not wind","not wind","not wind","not wind",
                      "not wind", "not wind","wind","wind","wind", "wind",
                      "wind", "not wind", "not wind", "wind", "wind","not wind",
                      "wind","wind","not wind","not wind","wind","wind","wind",
                      "wind","wind","wind","wind","wind","wind","not wind","wind",
                      "wind","not wind","wind","not wind","wind","wind","wind",
                      "wind","wind","wind","wind","wind","not wind","not wind",
                      "wind","wind","wind","wind","not wind","not wind","not wind",
                      "wind","not wind","not wind","wind","wind","not wind","wind",
                      "not wind","not wind"
          )) %>% 
  # Ericaceae specified as non wind pollinated in pollen data
  mutate(polmode = case_when(AccSpeciesName %in% c("Calluna vulgaris") ~ "not wind",
                                TRUE ~ polmode))

dfPFT <- full_join(pft, polmode, by = "AccSpeciesName") %>% 
  distinct()
saveRDS(dfPFT, "RDS_files/Polmode_pft_vegetation.rds")

# Make summary table per site -----
dfPFT <- readRDS("RDS_files/Polmode_pft_vegetation.rds")
zoneA <- dfABUN_a %>%
  rename(family = fam) %>% 
  ungroup %>% 
  mutate(zone = "zoneA",
         country = "Scotland")
zoneB <- dfABUN_bc %>% 
  ungroup() %>% 
  dplyr::select(sitename, stand.spec, genus, family = fam,
                abun = spec_abun_b) %>% 
  mutate(zone = "zoneB",
         country = "Scotland")
zoneC <- dfABUN_bc %>% 
  ungroup %>% 
  dplyr::select(sitename, stand.spec, genus, family = fam, 
                abun = spec_abun_c) %>% 
  mutate(zone = "zoneC",
         country = "Scotland")
bdm_zoneA <- bdm_abun %>%
  dplyr::select(sitename, stand.spec, genus, family, abun) %>%
  mutate(sitename = as.character(sitename)) %>%
  mutate(zone = "zoneA",
         country = "Switzerland")

dfABUN  <- bind_rows(zoneA, zoneB, zoneC,
                     bdm_zoneA)

sum_polmode <- dfABUN %>% 
  left_join(dfPFT, by = c("stand.spec" = "AccSpeciesName")) %>% 
  group_by(zone, sitename, polmode) %>% 
  # percentage of pollination type
  summarise(abun = sum(abun, na.rm = T)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  pivot_wider(names_from = polmode, values_from = abun)
  
sum <- dfABUN %>% 
  left_join(dfPFT, by = c("stand.spec" = "AccSpeciesName")) %>% 
  mutate(growthform = case_when(growthform %in% c("tree", "shrub") ~ "Woody",
                                growthform %in% c("herb", "grass", "fern") ~ "Non-Woody")
  ) %>% 
  group_by(zone, sitename, growthform) %>% 
  # percentage of growthform type
  summarise(abun = sum(abun)) %>% 
  mutate(abun = abun/sum(abun, na.rm = T)) %>% 
  pivot_wider(names_from = growthform, values_from = abun) %>% 
  # bind with pollination mode summary
  left_join(sum_polmode, by = c("zone", "sitename"), suffix = c("_pft", "_polmode")) %>% 
  mutate(across(where(is.numeric), ~round(.,2)*100))

saveRDS(sum, "RDS_files/01_Percentage_cover_pft.rds")

