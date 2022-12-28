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

gapfilled <- readRDS("RDS_files/02_Gapfilled_traits_pollen.rds") %>% 
  filter(country == "Switzerland") %>% 
  dplyr::select(stand.spec, pollentaxon, PlantHeight, SLA, LA) %>% 
  filter(!is.na(pollentaxon)) %>% 
  distinct() %>% 
  mutate(source = "Gapfilled TRY data") 
spec <- readRDS("RDS_files/02_PollenType_species.rds") %>% 
  dplyr::select(-country) %>% 
  distinct()

# Field collected traits for vegetation reconstruction
plantheight <- readRDS("RDS_files/01_PlantHeight.rds") %>% 
  ungroup() %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec,pollentaxon, PlantHeight) %>% 
  mutate(source = "Field Scotland")
sla <- readRDS("RDS_files/01_SLA.rds") %>% 
  ungroup() %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec,pollentaxon, SLA) %>% 
  mutate(source = "Field Scotland")
la <- readRDS("RDS_files/01_LA.rds") %>% 
  ungroup() %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec,pollentaxon, LA) %>% 
  mutate(source = "Field Scotland")

lTrait <- list(PlantHeight = plantheight,
               SLA = sla,
               LA = la)
# Load swiss data
bdm_ph <- readRDS("RDS_files/01_Traits_Swiss.rds") %>% 
  # filter trees, only measured up to 200 cm
  filter(!PFT == "TR") %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec, pollentaxon, PlantHeight) %>% 
  mutate(source = "Field Switzerland")
bdm_la <- readRDS("RDS_files/01_Traits_Swiss.rds") %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec, pollentaxon, LA) %>% 
  mutate(source = "Field Switzerland")
bdm_sla <- readRDS("RDS_files/01_Traits_Swiss.rds") %>% 
  left_join(spec, by = "stand.spec") %>% 
  select(stand.spec, pollentaxon, SLA) %>% 
  mutate(source = "Field Switzerland")


pollentaxa <- c("Betula","Ericales (tetrad)","Pinus","Poaceae","Pteridophyte","Ranunculaceae",
  "Rumex/Oxyria","Abies","Corylus","Fraxinus","Picea", "Asteraceae", "Apiaceae",
  "Cyperaceae","Rosaceae")

# plot
(ph_plot <- 
  # prepare data
    gapfilled %>% 
  bind_rows(plantheight, bdm_ph) %>% 
  dplyr::select(stand.spec, pollentaxon, PlantHeight, source) %>% 
  drop_na() %>% 
  filter(pollentaxon %in% pollentaxa) %>% 
    # Plot
  ggplot(aes(x = pollentaxon, y = log(PlantHeight), fill = source)) + 
  # geoms
  geom_boxplot(outlier.size = 0.1) +
  # scales and labelling
  scale_x_discrete("") +
  scale_y_continuous("Plant height (cm)(log)") +
  scale_fill_manual(values = c("darkorange", "purple", "cyan4")) +
  labs(fill = "Trait source") +
  coord_flip() +
  # Theme
  theme_bw())
ggsave("Figures/PlantHeight_trait_values.png", ph_plot, height = 5, width = 5)

(la_plot <- gapfilled %>% 
    bind_rows(la, bdm_la) %>% 
    dplyr::select(stand.spec, pollentaxon, LA, source) %>% 
    drop_na() %>% 
    filter(pollentaxon %in% pollentaxa) %>% 
  ggplot(aes(x = pollentaxon, y = log(LA), fill = source)) + 
  # geoms
  geom_boxplot(outlier.size =  0.1) +
  # scales and labelling
  scale_x_discrete("") +
  scale_y_continuous(expression(paste("Leaf area (",mm^2, ")(log)"))) + 
  scale_fill_manual(values = c("darkorange", "purple", "cyan4")) +
  labs(fill = "Trait source") +
  coord_flip() +
  # Theme
  theme_bw())
ggsave("Figures/LA_trait_values.png", la_plot, height = 5, width = 5)

(sla_plot <- gapfilled %>% 
    bind_rows(sla, bdm_sla) %>% 
    dplyr::select(stand.spec, pollentaxon, SLA, source) %>% 
    drop_na() %>% 
    filter(pollentaxon %in% pollentaxa) %>% 
  ggplot(aes(x = pollentaxon, y = log(SLA), fill = source)) + 
  # geoms
  geom_boxplot(outlier.size = 0.1) +
  # scales and labelling
  scale_x_discrete("") +
  scale_y_continuous(expression(paste("Specific leaf area (",mm^2, "/mg)(log)"))) + 
  scale_fill_manual(values = c("darkorange", "purple", "cyan4")) +
  labs(fill = "Trait source") +
  coord_flip() +
  # Theme
  theme_bw())
ggsave("Figures/SLA_trait_values.png", sla_plot, height = 5, width = 5)

