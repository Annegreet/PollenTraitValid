## ---------------------------
##
## Script name: Vegetation_pollen_comparison
##
## Purpose of script: Compare representation of species in the vegetation in the pollen
##
## Author: Annegreet Veeken
##
## Date Created: 2022-03-21
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   - 
##     
## ---------------------------

## Load packages
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(rioja)) install.packages("rioja") # plot pollen diagrams
if(!require(ggrepel)) install.packages("ggrepel") # repel point labels
if(!require(readxl)) install.packages("readxl") # repel point labels

## Prepare data ----
dfPOL_Scot <- readRDS("RDS_files/01_Pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_Pollen_data_Swiss.rds")
dfPOL <- bind_rows(dfPOL_Scot, dfPOL_Swiss)
polmode <- readRDS("RDS_files/01_Pollination_mode.rds")
polspec <- readRDS("RDS_files/02_PollenType_species.rds") %>% 
  dplyr::select(pollentaxon, stand.spec) %>% distinct() %>% as_tibble()
dfABUN_a <- readRDS("RDS_files/01_Species_abundance_a.rds")
dfABUN_bc <- readRDS("RDS_files/01_Species_abundance_bc.rds") %>% ungroup
bdm_abun <- readRDS("RDS_files/01_Species_abundance_Swiss.rds") 

# Summarize vegetation data by pollen type
zoneA <- dfABUN_a %>%
  # join with pollen to species table
  left_join(polspec, by = "stand.spec") %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  # select plots that are comparable with Switzerland
  filter(distance %in% c("0 meter", "1.5-3 meter")) %>% 
  # calculate species abundance on the plot level
  group_by(sitename, pollentaxon) %>% 
  summarise(abun = sum(abun)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  filter(!is.nan(abun)) %>%  
  mutate(zone = "zoneA",
         country = "Scotland")
zoneA_swiss <- bdm_abun %>% 
  # join with pollen to species table
  left_join(polspec, by = c("stand.spec" )) %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  # calculate species abundance on the plot level
  group_by(sitename, pollentaxon) %>% 
  summarise(abun = sum(abun)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  filter(!is.nan(abun)) %>% 
  mutate(sitename = as.character(sitename),
         zone = "zoneA",
         country = "Switzerland")
zoneB <- dfABUN_bc %>% 
  # join with pollen to species table
  left_join(polspec, by = c("stand.spec" )) %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  # calculate species abundance on the plot level
  group_by(sitename, pollentaxon) %>% 
  summarise(abun = sum(spec_abun_b)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  filter(!is.nan(abun)) %>% 
  mutate(zone = "zoneB",
         country = "Scotland")
zoneC <- dfABUN_bc %>% 
  # join with pollen to species table
  left_join(polspec, by = c("stand.spec" )) %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  # calculate species abundance on the plot level
  group_by(sitename, pollentaxon) %>% 
  summarise(abun = sum(spec_abun_c)) %>% 
  mutate(abun = abun/sum(abun)) %>% 
  filter(!is.nan(abun)) %>% 
  mutate(zone = "zoneC",
         country = "Scotland")

veg <- bind_rows(zoneA,zoneA_swiss, zoneB, zoneC) %>% 
  left_join(dfPOL_Scot, by = c("sitename","pollentaxon")) 

# Add pollination mode and plant functional type
# Harmonize nomenclature
polmode <- polmode %>% 
  mutate(pollentaxon = dplyr::recode(pollentaxon,  
                              "Apiaceae undiff." = "Apiaceae",
                              "Caryophyllaceae undiff." = "Caryophyllaceae",
                              "Cyperaceae undiff." = "Cyperaceae",
                              "Ericacea-type" = "Ericales (tetrad)",
                              "Fraxinus excelsior" = "Fraxinus",
                              "Geranium" = "Geraniaceae",
                              "Juniperus" = "Juniperus communis",
                              "Pedicularis" = "Pedicularis palustris",
                              "Ranunculuceae undiff." = "Ranunculaceae",
                              "Rosaceae undiff." = "Rosaceae",
                              "Viola palustris-type" = "Viola")) %>% 
  add_row(pollentaxon = "Asteraceae", fam = "Asteraceae", growthform = "herb", pollination = "not wind") %>% 
  add_row(pollentaxon = "Lamiaceae", fam = "Lamiaceae", growthform = "herb", pollination = "not wind") %>% 
  add_row(pollentaxon = "Malvaceae", fam = "Malvaceae", growthform = "tree", pollination = "not wind") %>% 
  add_row(pollentaxon = "Plantago", fam = "Plantaginaceae", growthform = "herb", pollination = "wind") %>% 
  add_row(pollentaxon = "Liliaceae")  %>% 
  add_row(pollentaxon = "Pteridophyte") %>% 
  add_row(pollentaxon = "Tsuga ")


## Calculate pollen representation value (Julier 2018, https://doi.org/10.1080/01916122.2017.1356392)
## Plot pollen representation per zone ----
zone_lab <- c("Zone A (0 - 3.4 m)", "Zone B (10 - 100 m)", "Zone C (100 - 1000 m)")
names(zone_lab) <- c("A", "B", "C")
(pol_veg_plot <- veg %>% 
  group_by(zone, pollentaxon) %>% 
  summarise(veg = mean(abun,na.rm = T), pol = median(percent,na.rm = T)) %>% 
  left_join(polmode, by = "pollentaxon") %>% 
ggplot(aes(x = veg, y = pol)) +
  geom_point(aes(col = pollination)) +
  # annotations
  geom_abline(intercept = 0, slope = 1) + 
  geom_text_repel(aes(label = pollentaxon), max.overlaps = 50, size = 3, segment.color = '#999999') +
  ggtitle("Pollen percentages") +
  scale_x_continuous("Vegetation basal area (%)") +
  scale_y_continuous("Pollen abundance (%)") + 
  scale_color_manual(name = "Pollination mode", 
                     values = c("darkorchid", "darkorange", "cyan4")) +
  # Faceting
  facet_wrap(~zone) +
  # Theme
  theme_bw(base_size = 14))

(pol_veg_plot_adj <- veg %>% 
  group_by(zone, pollentaxon) %>% 
  summarise(veg = mean(abun,na.rm = T), pol = median(adjustedpercent_mean,na.rm = T)) %>% 
  left_join(polmode, by = "pollentaxon")%>% 
ggplot(aes(x = veg, y = pol)) +
  geom_point(aes(col = pollination)) +
  # annotations
  geom_abline(intercept = 0, slope = 1) + 
  geom_text_repel(aes(label = pollentaxon), max.overlaps = 50, size = 3, segment.color = '#999999') +
  ggtitle("Adjusted pollen percentages") +
  scale_x_continuous("Vegetation basal area (%)") +
  scale_y_continuous("Pollen abundance (%)") + 
  scale_color_manual(name = "Pollination mode", 
                     values = c("darkorchid", "darkorange", "cyan4")) +
  # Faceting
  facet_wrap(~zone) +
  # Theme
  theme_bw(base_size = 14)
)
rm(list=setdiff(ls(), c("pol_veg_plot", "pol_veg_plot_adj")))