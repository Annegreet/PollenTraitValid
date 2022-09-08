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
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(data.table)) install.packages("data.table")

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
                             "Rhinanthus glacialis", "Potentilla brevifolia"),
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
                         "herb")
  ) %>% 
  # Ericaceae specified as herb in pollen data
  mutate(growthform = case_when(AccSpeciesName %in% c("Calluna vulgaris", 
                                                      "Erica cinerea",
                                                      "Erica tetralix") ~ "herb",
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
                             "Trifolium alpinum"),
          polmode = c("wind","wind", "wind", "wind", "not wind",
                      "wind","not wind", "wind", "not wind", "not wind",
                      "wind", "not wind", "not wind", "not wind", "wind",
                      "wind", "not wind", NA, "not wind", "not wind",
                      "not wind","not wind","not wind","not wind","not wind",
                      "not wind","not wind","wind","not wind","not wind",
                      "not wind","not wind","not wind","not wind","not wind",
                      "not wind", "not wind","wind","wind","wind", "wind",
                      "wind", "not wind", "not wind", "wind", "wind","not wind")) %>% 
  # Ericaceae specified as non wind pollinated in pollen data
  mutate(polmode = case_when(AccSpeciesName %in% c("Calluna vulgaris")~ "not wind",
                                TRUE ~ polmode))


dfPFT <- full_join(pft, polmode, by = "AccSpeciesName") %>% 
  distinct()
saveRDS(dfPFT, "RDS_files/Polmode_pft_vegetation.rds")

# Make summary table per site -----
zoneA <- dfABUN_a %>%
  arrange(sitename, stand.spec) %>% 
  # select plots that are comparable with Switzerland
  filter(distance %in% c("0 meter", "1.5-3 meter")) %>%
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  rename(family = fam) %>% 
  ungroup %>% 
  mutate(zone = "zoneA",
         country = "Scotland")
zoneB <- dfABUN_bc %>% 
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  ungroup() %>% 
  dplyr::select(sitename, stand.spec, genus, family = fam,
                abun = spec_abun_b, growthform, polmode)%>% 
  mutate(zone = "zoneB",
         country = "Scotland")
zoneC <- dfABUN_bc %>% 
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  ungroup %>% 
  dplyr::select(sitename, stand.spec, genus, family = fam, 
                abun = spec_abun_c, growthform, polmode) %>% 
  mutate(zone = "zoneC",
         country = "Scotland")
bdm_zoneA <- bdm_abun %>%
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  dplyr::select(sitename, stand.spec, genus, family, abun,
                growthform, polmode) %>%
  mutate(sitename = as.character(sitename)) %>%
  mutate(zone = "zoneA",
         country = "Switzerland")
bdm_abunbc <- readRDS("RDS_files/01_Species_abundance_zoneBC_swiss.rds")
bdm_zoneB <- bdm_abunbc %>%
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  ungroup() %>% 
  dplyr::select(sitename, stand.spec, genus, family, abun = spec_abun_b,
                growthform, polmode) %>%
  mutate(sitename = as.character(sitename)) %>%
  mutate(zone = "zoneB",
         country = "Switzerland")
bdm_zoneC <- bdm_abunbc %>%
  # join with pft and polmode data
  left_join(dfPFT, 
            by = c("stand.spec" = "AccSpeciesName")) %>% 
  ungroup() %>% 
  dplyr::select(sitename, stand.spec, genus, family, abun = spec_abun_c,
                growthform, polmode) %>%
  mutate(sitename = as.character(sitename)) %>%
  mutate(zone = "zoneC",
         country = "Switzerland")

dfABUN  <- bind_rows(zoneA, zoneB, zoneC,
                     bdm_zoneA, bdm_zoneB,bdm_zoneC)

sum_polmode <- dfABUN %>% 
  group_by(zone, sitename, stand.spec) %>% 
  summarise(abun = sum(abun)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  left_join(dfPFT, c("stand.spec" = "AccSpeciesName")) %>% 
  dplyr::select(zone, sitename, polmode,abun) %>% 
  drop_na() %>% 
  group_by(zone, sitename, polmode) %>% 
  summarise(percent = sum(abun, na.rm = TRUE)) %>% 
  pivot_wider(names_from = polmode, values_from = percent)
  
sum <- dfABUN %>% 
  group_by(zone,sitename, stand.spec) %>% 
  summarise(abun = sum(abun)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  left_join(dfPFT, c("stand.spec" = "AccSpeciesName")) %>% 
  drop_na %>% 
  mutate(growthform = case_when(growthform %in% c("tree", "shrub") ~ "Woody",
                                growthform %in% c("herb", "grass", "fern") ~ "Non-Woody")
         ) %>% 
  group_by(zone, sitename, growthform) %>% 
  summarise(percent = sum(abun, na.rm =TRUE)) %>% 
  pivot_wider(names_from = growthform, values_from = percent) %>%  
  left_join(sum_polmode, by = c("zone", "sitename")) %>% 
  replace_na(list(Woody = 0)) %>% 
  mutate(across(where(is.numeric), ~round(.,2)*100))
saveRDS(sum, "RDS_files/01_Percentage_cover_pft.rds")

