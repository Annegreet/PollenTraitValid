## ---------------------------
##
## Script name: 05_Scatterplots
##
## Purpose of script: Create scatterplots to go in to the appendix
##
## Author: Annegreet Veeken
##
## Date Created: 2022-09-08
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
if (!require(tidyverse)) install.packages("tidyverse")

## Load data
# load PFT summary table
# to filter sites with to little vegetation cover of particular type
sum_pft <- readRDS("RDS_files/01_Percentage_cover_pft.rds") %>% 
  replace_na(list(`Non-Woody` = 0,
                  Woody = 0,
                  `not wind` = 0,
                  wind = 0)) 

# collate data
traits <- c("LA","SLA","PlantHeight")

pollen_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "03_CWM_estimates_") %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(pattern = "pollen") %>% 
  str_subset(pattern = "notinpollen", negate = TRUE) %>% 
  str_subset(pattern = "draw", negate = TRUE)

pollen_lab <- pollen_files %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_remove(., "03_CWM_estimates_") %>% 
  str_remove(., ".rds")

veg_files <- list.files("RDS_files/") %>% 
  str_subset(pattern = "03_CWM_estimates_") %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_subset(pattern = "pollen", negate = TRUE) 

veg_lab <- veg_files %>% 
  str_subset(paste(traits, collapse = "|")) %>% 
  str_remove(., "03_CWM_estimates_") %>% 
  str_remove(., ".rds")

folderpath.fun <- function(x)
{paste("RDS_files/", x, sep = "/")}

dfCWM_pol <- pollen_files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(pollen_lab, ~readRDS(.) %>% 
                mutate(label = .y)) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(correction = case_when(str_detect(label, pattern = "adjustedpercent_draw") ~ str_extract(label, pattern = "draw[0-9]*$"),
                                str_detect(label, pattern = "adjustedpercent_mean") ~ "correction",
                                str_detect(label, pattern = "percent") ~ "no correction",
                                str_detect(label, pattern = "adjusted_helinger") ~ "helinger correction",
                                TRUE ~ "no correction"),
         pollination = replace_na(pollination, "allpol"),
         growthform = replace_na(growthform, "allpft"),
         sitename = str_remove(sitename, "X")
  ) %>%
  # don't include cwm values with removed species (supplementary analysis)
  # filter(is.na(removed_spec)) %>% 
  dplyr::select(country, sitename, pollen_mean = Mean, 
                pollen_sd = SD, trait, correction, pollination, growthform
  ) 


dfCWM_veg <- veg_files %>% 
  folderpath.fun(.) %>% 
  purrr::map2(veg_lab, ~readRDS(.) %>% 
                as.data.frame() %>% 
                mutate(label = .y) %>% 
                mutate(across(one_of("sitename"),  as.character))) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(country = str_extract(label, pattern = "Scotland|Switzerland"),
         pollination = dplyr::recode(pollination, `no wind` = "not wind"),
         notinpollen = if_else(str_detect(label, "notinpollen"), "notinpollen", "inpollen")) %>%
  # only Scottish data
  filter(country == "Scotland") %>% 
  # remove sites with to little cover of pft (5% or less) 
  left_join(sum_pft, by = c("sitename", "zone")) %>% 
  filter((growthform == "allpft" & pollination == "allpol") |
           (growthform == "trsh" & Woody >= 5) |
           (growthform == "herb" & `Non-Woody` >= 5) |
           (pollination == "wind" & wind >= 5) |
           (pollination == "not wind" & `not wind` >= 5)) %>% 
  dplyr::select(country, sitename, trait, veg_mean = Mean, veg_sd = SD,
                growthform, pollination, taxres, zone) 

saveRDS(dfCWM_pol,"RDS_files/04_CWM_estimates_pollen_Scotland.rds")
saveRDS(dfCWM_veg,"RDS_files/04_CWM_estimates_vegetation_Scotland.rds")

# join data
dfCWM <- dfCWM_veg %>%
  left_join(dfCWM_pol, by = c("country","sitename", "trait", "growthform","pollination")) %>%
  # create treatment column
  mutate(treatment = case_when((pollination == "allpol" & growthform == "allpft" &
                                  taxres == "stand.spec") ~ "Correction factor",
                               (taxres %in% c("family","genus") & correction == "no correction") ~ "Taxonomic resolution",
                               (pollination != "allpol" & growthform == "allpft" & correction == "no correction") ~ "Pollination mode",
                               (growthform != "allpft" & pollination == "allpol" & correction == "no correction") ~ "Growth form"),
         zone = str_remove(zone, pattern = "Scotland |Switzerland ")) %>%
  filter(country == "Scotland") 

# Plot labels
labs_trait <- as_labeller(c(PlantHeight = "Height~(cm)(log)",
                            LA = "Leaf~area~(cm^{2})(log)",
                            SLA ="SLA~(mm^{2}/mg)(log)"),
                          default = label_parsed)
labs_zone <- c(zoneA = "",
               zoneB = "",
               zoneC = "")

# gapfilled data pollen, field data + gapfilled vegetation ----
cwm_range <- dfCWM %>%
  group_by(trait) %>%
  summarise(across(c("veg_mean","pollen_mean"), ~range(., na.rm = TRUE))) %>%
  transmute(max_lim = max(veg_mean, pollen_mean, na.rm = TRUE),
            min_lim = min(veg_mean, pollen_mean, na.rm = TRUE)) %>%
  pivot_longer(-trait, names_to = "limits") %>%
  mutate(pollen_mean = value, veg_mean = value,
         correction = NA, taxres = NA, pollination = NA, growthform = NA)

(p_correction <-
    dfCWM %>%
    filter(country == "Scotland") %>%
    filter(treatment == "Correction factor") %>%
  ggplot(aes(x = pollen_mean, y = veg_mean, color = correction)) +
    geom_point(size = 1) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.1) +
    geom_blank(data = cwm_range) +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c("no correction" = "darkorange","correction" = "purple", 
                                  "helinger correction" = "cyan4"),
                       label = c("no correction" = "No correction","correction" = "Pollen productivity estimates", 
                                 "helinger correction" = "Helinger transformation")) +
    facet_wrap(zone~trait,
               scales = "free",
               labeller = labeller(trait = labs_trait,
                                   zone = labs_zone)
    ) +
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()))

ggsave("Figures/Scatter_correction.png", p_correction, height = 7, width = 7)


(p_pft <-
    dfCWM %>%
    filter(country == "Scotland") %>%
    filter(treatment == "Growth form" ) %>%
    ggplot(aes(x = pollen_mean, y = veg_mean, color = growthform)) +
    geom_point(size = 1) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.1) +
    geom_blank(data = cwm_range) +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c(herb = "darkorange", trsh = "purple"),
                       labels = c(herb = "Non-woody", trsh ="Woody")) +
    facet_wrap(zone~trait, scales = "free",
               labeller = labeller(trait = labs_trait,
                                   zone = labs_zone)) +
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()))
ggsave("Figures/Scatter_growthform.png", p_pft,
       height = 7, width = 7)

(p_pol <-
    dfCWM %>%
    filter(treatment == "Pollination mode") %>%
    ggplot(aes(x = pollen_mean, y = veg_mean, color = pollination)) +
    geom_point(size = 1) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.1) +
    geom_blank(data = cwm_range) +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c("not wind" = "darkorange", wind = "purple"),
                       labels = c("not wind" = "Not wind pollinated", wind = "Wind pollinated")) +
    facet_wrap(zone~trait, scales = "free",
               labeller = labeller(trait = labs_trait,
                                   zone = labs_zone)) +
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()))
ggsave("Figures/Scatter_pollination.png", p_pol,
       height = 7, width = 7)

(p_taxres <-
    dfCWM %>%
    filter(treatment == "Taxonomic resolution" |
             (treatment == "Correction factor" & correction == "no correction")) %>%
    ggplot(aes(x = pollen_mean, y = veg_mean, color = taxres)) +
    geom_point(size = 1) +
    geom_abline(intercept = 0, slope = 1, linewidth = 0.1) +
    geom_blank(data = cwm_range) +
    scale_x_continuous("Pollen") +
    scale_y_continuous("Vegetation") +
    scale_color_manual(values = c(family = "darkorange", genus = "purple", stand.spec = "cyan4"),
                       labels = c(family = "Family", genus = "Genus", stand.spec = "Species")) +
    facet_wrap(zone~trait, scales = "free",
               labeller = labeller(trait = labs_trait,
                                   zone = labs_zone)) +
    theme_bw() +
    theme(legend.position =  "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.title = element_blank()))
ggsave("Figures/Scatter_taxonomic_resolution.png", p_taxres,
       height = 7, width = 7)
