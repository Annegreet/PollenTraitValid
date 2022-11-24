## ---------------------------
##
## Script name: Evaluation of trait data
##
## Purpose of script: Comparing trait data collected in the field with TRY data
##
## Author: Annegreet Veeken
##
## Date Created: 2022-09-06
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

# Trait data from TRY for pollen based reconstruction
TRY <- readRDS("RDS_files/01_TRY_raw_traits.rds") %>% 
  dplyr::select(stand.spec, pollentaxon, PlantHeight, SLA, LA) %>% 
  filter(!is.na(pollentaxon)) %>% 
  distinct() %>% 
  mutate(source = "TRY") 
gapfilled <- readRDS("RDS_files/02_Gapfilled_traits_pollen.rds") %>% 
  dplyr::select(stand.spec, pollentaxon, PlantHeight, SLA, LA) %>% 
  filter(!is.na(pollentaxon)) %>% 
  distinct() %>% 
  mutate(source = "Gapfilled") 
spec <- readRDS("RDS_files/02_PollenType_species.rds") %>% 
  dplyr::select(-country) %>% 
  distinct()

# Field collected traits for vegetation reconstruction
plantheight <- readRDS("RDS_files/01_PlantHeight.rds") %>% 
  ungroup() %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec,pollentaxon, PlantHeight) %>% 
  mutate(source = "Field")
sla <- readRDS("RDS_files/01_SLA.rds") %>% 
  ungroup() %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec,pollentaxon, SLA) %>% 
  mutate(source = "Field")
la <- readRDS("RDS_files/01_LA.rds") %>% 
  ungroup() %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec,pollentaxon, LA) %>% 
  mutate(source = "Field")

lTrait <- list(PlantHeight = plantheight,
               LDMC = ldmc,
               SLA = sla,
               LA = la,
               LeafC = c,
               LeafN = n)
# Load swiss data
bdm_ph <- readRDS("RDS_files/01_Traits_Swiss.rds") %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec, pollentaxon, PlantHeight) %>% 
  mutate(source = "Field swiss")
bdm_la <- readRDS("RDS_files/01_Traits_Swiss.rds") %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec, pollentaxon, LA) %>% 
  mutate(source = "Field swiss")
bdm_sla <- readRDS("RDS_files/01_Traits_Swiss.rds") %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec, pollentaxon, SLA) %>% 
  mutate(source = "Field swiss")


pollentaxa <- c("Betula", "Carpinus betulus", "Cyperaceae", 
                                          "Ericales (tetrad)", "Juniperus", "Pinus", 
                                          "Poaceae", "Pteridophyte", "Ranunculaceae", 
                                          "Alnus", "Caryophyllaceae", "Rosaceae", 
                                          "Lamiaceae", "Rubiaceae", "Malvaceae", 
                                          "Plantago", "Viola", "Apiaceae", "Asteraceae", 
                                          "Picea", "Salix", "Geraniaceae", "Fraxinus", "Urtica",
                                          "Ulmus", "Acer", "Larix", "Abies", "Corylus", "Epilobium",
                                          "Quercus", "Rumex/Oxyria", "Castanea", "Frangula alnus", 
                                          "Juglans", "Tilia", "Fabaceae", "Lonicera", "Solanum",
                                          "Populus", "Veronica", "Convolvulaceae" )
# plot
(ph_plot <- TRY %>% 
  bind_rows(gapfilled, plantheight, bdm_ph) %>% 
  dplyr::select(stand.spec, pollentaxon, PlantHeight, source) %>% 
  drop_na() %>% 
  filter(pollentaxon %in% pollentaxa) %>% 
  ggplot(aes(x = pollentaxon, y = log(PlantHeight), fill = source)) + 
  # geoms
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = c("darkorange", "purple", "cyan4", "white")) +
  labs(fill = "Trait source") +
  coord_flip() +
  # # Faceting
  # facet_wrap(~trait, labeller = labeller(trait = trait_lab)) +
  # Theme
  theme_bw())
ggsave("Figures/PlantHeight_trait_values.png", ph_plot, height = 9, width = 5)

(la_plot <-  TRY %>% 
    bind_rows(gapfilled, la, bdm_la) %>% 
    dplyr::select(stand.spec, pollentaxon, LA, source) %>% 
    drop_na() %>% 
    filter(pollentaxon %in% pollentaxa) %>% 
  ggplot(aes(x = pollentaxon, y = log(LA), fill = source)) + 
  # geoms
  geom_boxplot(outlier.size =  0.1) +
  scale_fill_manual(values = c("darkorange", "purple", "cyan4", "white")) +
  labs(fill = "Trait source") +
  coord_flip() +
  # # Faceting
  # facet_wrap(~trait, labeller = labeller(trait = trait_lab)) +
  # Theme
  theme_bw())
ggsave("Figures/LA_trait_values.png", la_plot, height = 9,width = 5)
(sla_plot <-  TRY %>% 
    bind_rows(gapfilled, sla, bdm_sla) %>% 
    dplyr::select(stand.spec, pollentaxon, SLA, source) %>% 
    drop_na() %>% 
    filter(pollentaxon %in% pollentaxa) %>% 
  ggplot(aes(x = pollentaxon, y = log(SLA), fill = source)) + 
  # geoms
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = c("darkorange", "purple", "cyan4", "white")) +
  labs(fill = "Trait source") +
  coord_flip() +
  # # Faceting
  # facet_wrap(~trait, labeller = labeller(trait = trait_lab)) +
  # Theme
  theme_bw())
ggsave("Figures/SLA_trait_values.png", sla_plot, height = 9,width = 5)

