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
## References: 
## - pollen representation plot (Julier 2018, https://doi.org/10.1080/01916122.2017.1356392)
## Notes:
##  
##     
## ---------------------------

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggrepel)) install.packages("ggrepel") # repel point labels
if (!require(readxl)) install.packages("readxl") 

## Prepare data ----
dfPOL_Scot <- readRDS("RDS_files/01_Pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_Pollen_data_Swiss.rds")
dfPOL <- bind_rows("Scotland" = dfPOL_Scot, "Switzerland" = dfPOL_Swiss, .id = "country") %>% 
  mutate(sitename = str_remove(sitename, pattern = "X"))
polmode <- readRDS("RDS_files/01_Pollination_mode.rds")
polspec <- readRDS("RDS_files/02_PollenType_species.rds") %>% 
  dplyr::select(pollentaxon, stand.spec) %>% distinct() %>% as_tibble()
dfABUN_a <- readRDS("RDS_files/01_Species_abundance_a.rds")
dfABUN_bc <- readRDS("RDS_files/01_Species_abundance_bc.rds") %>% ungroup
bdm_abun <- readRDS("RDS_files/01_Species_abundance_Swiss.rds") 

# Summarize vegetation data by pollen type
zoneA <- dfABUN_a %>%
  # join with pollen to species table
  full_join(polspec, by = "stand.spec") %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  # calculate species abundance on the plot level
  group_by(sitename, pollentaxon) %>% 
  summarise(abun = sum(abun, na.rm = TRUE)) %>% 
  mutate(abun = abun/sum(abun, na.rm = TRUE)) %>% 
  mutate(zone = "zoneA",
         country = "Scotland")
zoneA_swiss <- bdm_abun %>% 
  # join with pollen to species table
  full_join(polspec, by = c("stand.spec" )) %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  # calculate species abundance on the plot level
  group_by(sitename, pollentaxon) %>% 
  summarise(abun = sum(abun, na.rm = TRUE)) %>% 
  mutate(abun = abun/sum(abun, na.rm = TRUE)) %>% 
  filter(!is.nan(abun)) %>% 
  mutate(sitename = as.character(sitename),
         zone = "Switzerland",
         country = "Switzerland")
zoneB <- dfABUN_bc %>% 
  # join with pollen to species table
  full_join(polspec, by = c("stand.spec" )) %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  # calculate species abundance on the plot level
  group_by(sitename, pollentaxon) %>% 
  summarise(abun = sum(spec_abun_b, na.rm = TRUE)) %>% 
  mutate(abun = abun/sum(abun, na.rm = TRUE)) %>% 
  filter(!is.nan(abun)) %>% 
  mutate(zone = "zoneB",
         country = "Scotland")
zoneC <- dfABUN_bc %>% 
  # join with pollen to species table
  full_join(polspec, by = c("stand.spec" )) %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  # calculate species abundance on the plot level
  group_by(sitename, pollentaxon) %>% 
  summarise(abun = sum(spec_abun_c, na.rm = T)) %>% 
  mutate(abun = abun/sum(abun, na.rm = T)) %>% 
  filter(!is.nan(abun)) %>% 
  mutate(zone = "zoneC",
         country = "Scotland")

veg <- bind_rows(zoneA, zoneA_swiss, zoneB, zoneC) %>% 
  left_join(dfPOL, by = c("country","sitename","pollentaxon")) 

## Plot pollen representation per zone ----
zone_lab <- c("Switzerland (1.8 m)", "Inner ring (10 m)", "Middle ring (100 m)", "Outer ring (1km)")
names(zone_lab) <- c("Switzerland","zoneA", "zoneB", "zoneC")

nopol <- veg %>% 
  filter(pollentaxon == "Not in pollen data") %>% 
  group_by(country, zone) %>% 
  summarise(nopol_mean = round(mean(abun, na.rm = T), 2)*100,
            nopol_sd = round(sd(abun, na.rm = T), 2)*100) %>% 
  ungroup() 
pol_veg <- veg %>% 
  group_by(country, zone, pollentaxon) %>% 
  summarise(veg = mean(abun, na.rm = T)*100, 
            pol = mean(percent, na.rm = T)*100) %>% 
  filter(!(is.nan(veg) & is.nan(pol))) %>% 
  mutate(r_rel = pol/veg,
         r_rel_abs =  case_when(r_rel == Inf ~ "not in vegetation data",
                            is.nan(r_rel) ~ "not in pollen data")
            ) %>% 
  left_join(polmode, by = "pollentaxon") %>% 
  filter(!is.na(pollination)) 

   
(p <- pol_veg %>% 
  group_by(pollentaxon) %>% 
  filter(any(pol > 5)) %>% 
  ggplot(aes(x = veg, y = pol)) +
  geom_point(aes(col = pollination)) +
  # annotations
  geom_abline(intercept = 0, slope = 1) + 
  geom_text_repel(aes(label = pollentaxon),
                  max.overlaps = 50, size = 3, segment.color = '#999999') +
  scale_x_continuous("Mean vegetation basal area (%)", limits = c(0,100)) +
  scale_y_continuous("Mean pollen abundance (%)", limits = c(0,100)) + 
  scale_color_manual(name = "Pollination mode", 
                     values = c("darkorchid", "darkorange"),
                     labels = c(`not wind` = "Not wind pollinated", wind = "Wind pollinated")) +
  # Faceting
  facet_wrap(~zone, labeller = labeller(zone = zone_lab)) +
  # Theme
  theme_bw() +
  theme(text = element_text(size = 10),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom"
  ))
ggsave("Figures/Pollen_rep_percent.png", p, height = 7, width = 7)


pol_veg_adj <- veg %>% 
  group_by(country, zone, pollentaxon) %>% 
  summarise(veg = mean(abun,na.rm = T)*100, 
            pol = mean(adjustedpercent_mean, na.rm = T)*100) %>% 
  left_join(polmode, by = "pollentaxon") %>% 
  filter(!is.na(pollination)) %>% 
  group_by(pollentaxon) %>% 
  filter(any(pol > 5)) 

pol_veg_draw <- veg %>% 
  group_by(country, zone, pollentaxon) %>% 
  summarise(veg = mean(abun,na.rm = T)*100, 
            across(contains("adjustedpercent_draw"), ~mean(., na.rm = T)*100,  .names = "{.col}")) %>% 
  left_join(polmode, by = "pollentaxon") %>% 
  filter(!is.na(pollination)) %>% 
  ungroup() %>% 
  dplyr::select(zone, pollentaxon, veg, contains("adjusted"), pollination) %>%
  pivot_longer(cols = contains("adjusted"), names_to = "draw", values_to = "pol") %>% 
  filter(pollentaxon %in% unique(pol_veg_adj$pollentaxon))

(p2 <- ggplot() +
  # draws
  geom_point(data = pol_veg_draw, aes(x = veg, y = pol), alpha = 0.2) +
  # mean adjusted values
  geom_point(data = pol_veg_adj, aes(x = veg, y = pol, color = pollination), size = 1.75) +
  # annotations
  geom_abline(intercept = 0, slope = 1) + 
  geom_text_repel(data = pol_veg_adj, aes(x = veg, y = pol, label = pollentaxon),
                  max.overlaps = 50, size = 3, segment.color = '#999999') +
  scale_x_continuous("Mean vegetation basal area (%)", limits = c(0,100)) +
  scale_y_continuous("Mean pollen abundance (%)", limits = c(0,100)) + 
  scale_color_manual(name = "Pollination mode", 
                     values = c("darkorchid", "darkorange"),
                     labels = c(`not wind` = "Not wind pollinated", wind = "Wind pollinated")) +
  # Faceting
  facet_wrap(~zone, labeller = labeller(zone = zone_lab)) +
  # Theme
  theme_bw() +
  theme(text = element_text(size = 10),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom"
    ))

ggsave("Figures/Pollen_rep_adjustedpercent.png", p2, height = 7, width = 7)

## Trait values and representation of taxa  ----
dfTRAIT <- readRDS("RDS_files/02_Gapfilled_traits.rds")
la <- readRDS("RDS_files/03_taxon_trait_estimates_Scotland_percent_LA.rds") 

# calculate r-rel values
r_rel <- zoneB %>%
  full_join(dfPOL_Scot, by = c("sitename", "pollentaxon")) %>% 
  filter(!abun == 0) %>% 
  filter(pollentaxon == "Not in pollen data") %>% 
  # calculate average representation 
  group_by(pollentaxon) %>% 
  summarise(veg = mean(abun, na.rm = T), 
            pol = mean(percent, na.rm = T)) %>% 
  mutate(r_rel = pol/veg) %>% 
  mutate(repres = case_when(is.nan(r_rel)  ~ "Not in pollen data",
                            r_rel > 1 ~ "Overrepresented in pollen",
                            r_rel < 1 ~"Underrepresented in pollen"))
## add taxa that are not in pollen data
notinpollen <- dfABUN_bc %>% 
  # join with pollen to species table
  full_join(polspec, by = c("stand.spec" )) %>% 
  mutate(pollentaxon = if_else(is.na(pollentaxon), "Not in pollen data", pollentaxon)) %>% 
  filter(pollentaxon == "Not in pollen data") %>% 
  filter(!spec_abun_b == 0)

ggplot(r_rel, aes(x = repres, y = Mean)) +
  geom_boxplot() +
  facet_grid(~zone)
r_rel %>% 
  mutate(r_rel = if_else(is.na(r_rel), 0, r_rel),
         r_rel = if_else(r_rel == Inf, 1000, r_rel)) %>% 
ggplot( aes(x = r_rel, y = Mean)) +
  geom_point() +
  facet_grid(~zone)
