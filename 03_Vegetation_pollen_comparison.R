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
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(rioja)) install.packages("rioja") # plot pollen diagrams
if (!require(ggrepel)) install.packages("ggrepel") # repel point labels
if (!require(readxl)) install.packages("readxl") # repel point labels

## Prepare data ----
dfPOL_Scot <- readRDS("RDS_files/01_Pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_Pollen_data_Swiss.rds")
dfPOL <- bind_rows("Scotland" = dfPOL_Scot, "Switzerland" = dfPOL_Swiss, .id = "country") %>% 
  dplyr::select(country,sitename, pollentaxon, percent, adjustedpercent_mean) %>% 
  mutate(sitename = str_remove(sitename, pattern = "X"))
polmode <- readRDS("RDS_files/01_Pollination_mode.rds")
polspec <- readRDS("RDS_files/02_PollenType_species.rds") %>% 
  dplyr::select(pollentaxon, stand.spec) %>% distinct() %>% as_tibble()
dfABUN_a <- readRDS("RDS_files/01_Species_abundance_a.rds")
dfABUN_bc <- readRDS("RDS_files/01_Species_abundance_bc.rds") %>% ungroup
bdm_abun <- readRDS("RDS_files/01_Species_abundance_Swiss.rds") 
bdm_abun_bc <- readRDS("RDS_files/01_Species_abundance_zoneBC_swiss.rds")

# Summarize vegetation data by pollen type
zoneA <- dfABUN_a %>%
  # join with pollen to species table
  left_join(polspec, by = "stand.spec") %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  # select plots that are comparable with Switzerland
  filter(distance %in% c("0 meter", "1.5-3 meter")) %>% 
  # calculate species abundance on the plot level
  group_by(sitename, pollentaxon) %>% 
  summarise(abun = sum(abun, na.rm = TRUE)) %>% 
  mutate(abun = abun/sum(abun, na.rm = TRUE)) %>% 
  mutate(zone = "zoneA",
         country = "Scotland")
zoneA_swiss <- bdm_abun %>% 
  # join with pollen to species table
  left_join(polspec, by = c("stand.spec" )) %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  # calculate species abundance on the plot level
  group_by(sitename, pollentaxon) %>% 
  summarise(abun = sum(abun, na.rm = TRUE)) %>% 
  mutate(abun = abun/sum(abun, na.rm = TRUE)) %>% 
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
  summarise(abun = sum(spec_abun_b, na.rm = TRUE)) %>% 
  mutate(abun = abun/sum(abun, na.rm = TRUE)) %>% 
  filter(!is.nan(abun)) %>% 
  mutate(zone = "zoneB",
         country = "Scotland")
zoneB_swiss <- bdm_abun_bc %>% 
  # join with pollen to species table
  left_join(polspec, by = c("stand.spec" )) %>% 
  ungroup() %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  dplyr::select(sitename, pollentaxon, abun = spec_abun_b) %>% 
  mutate(sitename = as.character(sitename),
         zone = "zoneB",
         country = "Switzerland")

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
zoneC_swiss <- bdm_abun_bc %>% 
  # join with pollen to species table
  left_join(polspec, by = c("stand.spec" )) %>% 
  ungroup() %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  dplyr::select(sitename, pollentaxon, abun = spec_abun_c) %>% 
  mutate(sitename = as.character(sitename),
         zone = "zoneC",
         country = "Switzerland")

veg <- bind_rows(zoneA, zoneA_swiss, zoneB, zoneB_swiss, zoneC, zoneC_swiss) %>% 
  left_join(dfPOL, by = c("country","sitename","pollentaxon")) 

# Add pollination mode and plant functional type
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
zone_lab <- c("Inner ring (3.4 m)", "Middle ring (100 m)", "Outer ring (1km)")
names(zone_lab) <- c("zoneA", "zoneB", "zoneC")
windows()
nopol <- veg %>% 
  filter(pollentaxon == "Not in pollen data") %>% 
  group_by(country, zone) %>% 
  summarise(nopol_mean = round(mean(abun, na.rm = T),2)*100,
            nopol_sd = round(sd(abun, na.rm = T),2)*100) %>% 
  ungroup() %>% 
  transmute(veg = 0.75, pol = rep(c(1,0.9), 3), zone = factor(zone),
            country = country,
            lab = paste0(country,": ", nopol_mean,"%"," ± ",nopol_sd)) 
pol_veg <- veg %>% 
  group_by(country, zone, pollentaxon) %>% 
  summarise(veg = mean(abun,na.rm = T)*100, 
            pol = mean(percent,na.rm = T)*100) %>% 
  left_join(polmode, by = "pollentaxon") %>% 
  filter(!is.na(pollination)) %>% 
  group_by(pollentaxon) %>% 
  filter(any(pol > 10)) 
p <- ggplot(pol_veg, aes(x = veg, y = pol, shape = country)) +
  geom_point(aes(col = pollination)) +
  # annotations
  geom_abline(intercept = 0, slope = 1) + 
  geom_text_repel(aes(label = pollentaxon),
                  max.overlaps = 50, size = 3, segment.color = '#999999') +
  ggtitle("Pollen percentages") +
  scale_x_continuous("Vegetation basal area (%)", limits = c(0,100)) +
  scale_y_continuous("Pollen abundance (%)", limits = c(0,100)) + 
  scale_color_manual(name = "Pollination mode", 
                     values = c("darkorchid", "darkorange"),
                     labels = c("Not wind pollinated", "Wind pollinated")) +
  # Faceting
  facet_wrap(~zone, labeller = labeller(zone = zone_lab)) +
  # Theme
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) 
ggsave("Figures/Pollen_rep_percent.png", p, height = 5, width = 7)

pol_veg_adj <- veg %>% 
  group_by(country, zone, pollentaxon) %>% 
  summarise(veg = mean(abun,na.rm = T)*100, 
            pol = mean(adjustedpercent_mean,na.rm = T)*100) %>% 
  left_join(polmode, by = "pollentaxon") %>% 
  filter(!is.na(pollination)) %>% 
  group_by(pollentaxon) %>% 
  filter(any(pol > 10)) 
p2 <- ggplot(pol_veg_adj, aes(x = veg, y = pol, shape = country)) +
  geom_point(aes(col = pollination)) +
  # annotations
  geom_abline(intercept = 0, slope = 1) + 
  geom_text_repel(aes(label = pollentaxon),
                  max.overlaps = 50, size = 3, segment.color = '#999999') +
  ggtitle("Adjusted pollen percentages") +
  scale_x_continuous("Vegetation basal area (%)", limits = c(0,100)) +
  scale_y_continuous("Pollen abundance (%)", limits = c(0,100)) + 
  scale_color_manual(name = "Pollination mode", 
                     values = c("darkorchid", "darkorange"),
                     labels = c("Not wind pollinated", "Wind pollinated")) +
  # Faceting
  facet_wrap(~zone, labeller = labeller(zone = zone_lab)) +
  # Theme
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        legend.title = element_blank())  
p2
ggsave("Figures/Pollen_rep_adjustedpercent.png", p2, height = 5, width = 7)
