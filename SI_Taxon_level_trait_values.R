## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Annegreet Veeken
##
## Date Created: 2022-11-03
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   
## References:
##
## ---------------------------

set.seed(1)
options(scipen = 9999)

## Load packages
if (!require(tidyverse)) install.packages("tidyverse")

## Load data
traits <- c("PlantHeight", "SLA", "LA")
files <- list.files("RDS_files/", full.names = TRUE) %>% 
  str_subset(pattern = "03_Species|03_taxon") 

label <- files %>% 
  str_remove("RDS_files/03_")
  
dfESTIM <- files %>% 
  purrr::map2(label, ~readRDS(.) %>% 
                mutate(label = .y)) %>% 
  bind_rows() %>% 
  as_tibble() %>% 
  mutate(taxon = ifelse(str_detect(label, "zoneB"), "vegetation", "pollen"),
         trait = str_extract(label, pattern = paste(traits, collapse = "|")),
         country = str_extract(label, pattern = "Scotland|Switzerland"),
         reduced_list = ifelse(str_detect(label, pattern = "field"), "field", "GBIF"))

dfPOL <- readRDS("RDS_files/01_Pollen_data_Scot.rds")
dfTRAIT <- readRDS("RDS_files/02_Gapfilled_traits_pollen.rds") %>% 
  filter(country == "Scotland")

pollentax <- dfPOL %>% 
  group_by(pollentaxon) %>% 
  filter(any(percent > 0.08)) %>% 
  pull(pollentaxon) %>% 
  unique()

pollen_lab <- as_labeller(c("Betula" = "italic(Betula)",
                            "Cyperaceae" = "Cyperaceae",
                            "Ericales (tetrad)" = "Ericales~(tetrad)",
                            "Pinus" = "italic(Pinus)",
                            "Poaceae" = "Poaceae",
                            "Pteridophyte" = "Pteridophyte"), default = label_parsed)
# Pollen ----
## LA
estimated <- dfESTIM %>% 
  filter(trait == "LA") %>% 
  filter(taxon == "pollen") %>% 
  filter(country == "Scotland") %>% 
  group_split(tax, reduced_list)
names(estimated) <- dfESTIM %>% filter(taxon == "pollen") %>% filter(country == "Scotland") %>% 
  mutate(label = paste(.$tax, .$reduced_list, sep = "_")) %>% pull(label) %>% unique

df <- estimated %>% 
  purrr::map_dfr(~rnorm(20000, mean = .$mean.tax, sd = .$sd.tax)) %>% 
  pivot_longer(everything(), names_to = "pollentaxon", values_to = "LA") %>% 
  separate(pollentaxon, into = c("pollentaxon", "sp_list"), sep = "_")

(p <- ggplot(data = df[df$pollentaxon %in% pollentax,], aes(x = LA, fill = sp_list, color = sp_list)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual("", values = c("darkorange", "purple"), labels = c("field" = "Field species list",
                                                                       "GBIF" = "GBIF species list")) +
  scale_color_manual("", values = c("darkorange", "purple"), labels = c("field" = "Field species list",
                                                                         "GBIF" = "GBIF species list")) +
  xlab(parse(text = "Leaf~area~(cm^2)(log)")) +
  facet_wrap(~pollentaxon, scales = "free_y", labeller = pollen_lab) + 
  theme_bw() +
  theme(text = element_text(size = 10),
          strip.background = element_rect(fill = "white"),
        legend.position = "bottom")
)
ggsave("Figures/Trait_values_pollen_LA.png", p, height = 4, width = 7)

## SLA
estimated <- dfESTIM %>% 
  filter(trait == "SLA") %>% 
  filter(taxon == "pollen") %>% 
  filter(country == "Scotland") %>% 
  group_split(tax, reduced_list)
names(estimated) <- dfESTIM %>% filter(taxon == "pollen") %>% filter(country == "Scotland") %>% 
  mutate(label = paste(.$tax, .$reduced_list, sep = "_")) %>% pull(label) %>% unique

df <- estimated %>% 
  purrr::map_dfr(~rnorm(20000, mean = .$mean.tax, sd = .$sd.tax)) %>% 
  pivot_longer(everything(), names_to = "pollentaxon", values_to = "SLA") %>% 
  separate(pollentaxon, into = c("pollentaxon", "sp_list"), sep = "_")

(p <- ggplot(data = df[df$pollentaxon %in% pollentax,], aes(x = SLA, fill = sp_list, color = sp_list)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual("", values = c("darkorange", "purple"), labels = c("field" = "Field species list",
                                                                         "GBIF" = "GBIF species list")) +
    scale_color_manual("", values = c("darkorange", "purple"), labels = c("field" = "Field species list",
                                                                          "GBIF" = "GBIF species list")) +
    xlab(parse(text = "SLA~(mm^2/mg)(log)")) +
    facet_wrap(~pollentaxon, scales = "free_y", labeller = pollen_lab) + 
    theme_bw() +
    theme(text = element_text(size = 10),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom")
)

ggsave("Figures/Trait_values_pollen_SLA.png", p, height = 4, width = 7)

## Plant Height
estimated <- dfESTIM %>% 
  filter(trait == "PlantHeight") %>% 
  filter(taxon == "pollen") %>% 
  filter(country == "Scotland") %>% 
  group_split(tax, reduced_list)
names(estimated) <- dfESTIM %>% filter(taxon == "pollen") %>% filter(country == "Scotland") %>% 
  mutate(label = paste(.$tax, .$reduced_list, sep = "_")) %>% pull(label) %>% unique

df <- estimated %>% 
  purrr::map_dfr(~rnorm(20000, mean = .$mean.tax, sd = .$sd.tax)) %>% 
  pivot_longer(everything(), names_to = "pollentaxon", values_to = "PlantHeight") %>% 
  separate(pollentaxon, into = c("pollentaxon", "sp_list"), sep = "_")

(p <- ggplot(data = df[df$pollentaxon %in% pollentax,], aes(x = PlantHeight, fill = sp_list, color = sp_list)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual("", values = c("darkorange", "purple"), labels = c("field" = "Field species list",
                                                                         "GBIF" = "GBIF species list")) +
    scale_color_manual("", values = c("darkorange", "purple"), labels = c("field" = "Field species list",
                                                                          "GBIF" = "GBIF species list")) +
    xlab(parse(text = "Heigh (cm)(log)")) +
    facet_wrap(~pollentaxon, scales = "free_y", labeller = pollen_lab) + 
    theme_bw() +
    theme(text = element_text(size = 10),
          strip.background = element_rect(fill = "white"),
          legend.position = "bottom")
)

ggsave("Figures/Trait_values_pollen_PH.png", p, height = 4, width = 7)


# Vegetation ----
dfABUN_bc <- readRDS("RDS_files/01_Species_abundance_bc.rds")
veg_sp <- dfABUN_bc %>% 
  group_by(stand.spec) %>% 
  filter(any(spec_abun_b > 10)) %>% 
  pull(stand.spec) %>% 
  unique()

plantheight <- readRDS("RDS_files/01_PlantHeight.rds") %>% 
  ungroup()  %>% 
  filter(stand.spec %in% veg_sp)
sla <- readRDS("RDS_files/01_SLA.rds") %>% 
  ungroup() %>% 
  filter(stand.spec %in% veg_sp)
la <- readRDS("RDS_files/01_LA.rds") %>% 
  ungroup() %>% 
  filter(stand.spec %in% veg_sp)

lTrait <- list(PlantHeight = plantheight,
               SLA = sla,
               LA = la)

dfTRY_raw <- readRDS("RDS_files/02_Gapfilled_traits.rds")
lTrait <-
  lTrait %>% 
  purrr::map(., ~mutate(., stand.spec = 
                          case_when(stand.spec %in% veg_sp ~ stand.spec,
                                    genus %in% veg_sp ~ genus,
                                    family %in% veg_sp ~ family)
  )
  ) 

dfTRY <- 
  dfTRY_raw %>% 
  mutate(., stand.spec = case_when(stand.spec %in% veg_sp ~ stand.spec,
                                   genus %in% veg_sp ~ genus,
                                   family %in% veg_sp ~ family)
  )

## LA
estimated <- dfESTIM %>% 
  filter(trait == "LA") %>% 
  filter(taxon == "vegetation") %>% 
  filter(country == "Scotland") %>% 
  group_split(tax)
names(estimated) <- dfESTIM %>% filter(taxon == "vegetation") %>% filter(country == "Scotland") %>% pull(tax) %>% unique

df <- estimated %>% 
  purrr::map_dfr(~rnorm(20000, mean = .$mean.tax, sd = .$sd.tax)) %>% 
  pivot_longer(everything(), names_to = "stand.spec", values_to = "LA")

# Trait data
traittax <- lTrait$LA %>% 
  pull(stand.spec) %>% 
  unique()
missingtax <- veg_sp[!veg_sp %in% traittax] 

# filter species with too little observations (10 or less)
nobs <- lTrait$LA %>%
  group_by(stand.spec) %>%
  summarise(n = n()) %>%
  filter(n <= 10) %>%
  pull(stand.spec)

# find taxa with not or too little data in TRY gapfilled data
dfTRY_sel <- dfTRY %>% 
  filter(stand.spec %in% c(missingtax, nobs)) %>% 
  dplyr::select(stand.spec, LA)

dfTRAIT <- lTrait$LA %>% 
  dplyr::select(stand.spec, LA) %>% 
  bind_rows(dfTRY_sel)

(p <- ggplot(data = df[df$stand.spec %in% veg_sp,], aes(x = LA)) +
    geom_density() +
    geom_histogram(data = dfTRAIT[dfTRAIT$stand.spec %in% veg_sp,],
                   aes(x = log(LA), y = ..density..)) +
    xlab("LA (cm2)(log)") +
    facet_wrap(~stand.spec, scales = "free", nrow = 7) +
    theme_bw() +
    theme(strip.text = element_text(face = "italic"))
)
ggsave("Figures/Trait_values_vegetation_LA.png", p, height = 8, width = 6)

## SLA
estimated <- dfESTIM %>% 
  filter(trait == "SLA") %>% 
  filter(taxon == "vegetation") %>% 
  filter(country == "Scotland") %>% 
  group_split(tax)
names(estimated) <- dfESTIM %>% filter(taxon == "vegetation") %>% filter(country == "Scotland") %>% pull(tax) %>% unique

df <- estimated %>% 
  purrr::map_dfr(~rnorm(20000, mean = .$mean.tax, sd = .$sd.tax)) %>% 
  pivot_longer(everything(), names_to = "stand.spec", values_to = "SLA")

# Trait data
traittax <- lTrait$SLA %>% 
  pull(stand.spec) %>% 
  unique()
missingtax <- veg_sp[!veg_sp %in% traittax] 

# filter species with too little observations (10 or less)
nobs <- lTrait$SLA %>%
  group_by(stand.spec) %>%
  summarise(n = n()) %>%
  filter(n <= 10) %>%
  pull(stand.spec)

# find taxa with not or too little data in TRY gapfilled data
dfTRY_sel <- dfTRY %>% 
  filter(stand.spec %in% c(missingtax, nobs)) %>% 
  dplyr::select(stand.spec, SLA)

dfTRAIT <- lTrait$SLA %>% 
  dplyr::select(stand.spec, SLA) %>% 
  bind_rows(dfTRY_sel)

(p <- ggplot(data = df[df$stand.spec %in% veg_sp,], aes(x = SLA)) +
    geom_density() +
    geom_histogram(data = dfTRAIT[dfTRAIT$stand.spec %in% veg_sp,],
                   aes(x = log(SLA), y = ..density..)) +
    xlab("SLA (mm2/mg)(log)") +
    facet_wrap(~stand.spec, scales = "free", nrow = 7) +
    theme_bw() +
    theme(strip.text = element_text(face = "italic")))
ggsave("Figures/Trait_values_vegetation_SLA.png", p, height = 8, width = 6)

## Plant Height
estimated <- dfESTIM %>% 
  filter(trait == "PlantHeight") %>% 
  filter(taxon == "vegetation") %>% 
  filter(country == "Scotland") %>% 
  group_split(tax)
names(estimated) <- dfESTIM %>% filter(taxon == "vegetation") %>% filter(country == "Scotland") %>% pull(tax) %>% unique

df <- estimated %>% 
  purrr::map_dfr(~rnorm(20000, mean = .$mean.tax, sd = .$sd.tax)) %>% 
  pivot_longer(everything(), names_to = "stand.spec", values_to = "PlantHeight")

# Trait data
traittax <- lTrait$PlantHeight %>% 
  pull(stand.spec) %>% 
  unique()
missingtax <- veg_sp[!veg_sp %in% traittax] 

# filter species with too little observations (10 or less)
nobs <- lTrait$PlantHeight %>%
  group_by(stand.spec) %>%
  summarise(n = n()) %>%
  filter(n <= 10) %>%
  pull(stand.spec)

# find taxa with not or too little data in TRY gapfilled data
dfTRY_sel <- dfTRY %>% 
  filter(stand.spec %in% c(missingtax, nobs)) %>% 
  dplyr::select(stand.spec, PlantHeight)

dfTRAIT <- lTrait$PlantHeight %>% 
  dplyr::select(stand.spec, PlantHeight) %>% 
  bind_rows(dfTRY_sel)

(p <- ggplot(data = df[df$stand.spec %in% veg_sp,], aes(x = PlantHeight)) +
    geom_density() +
    geom_histogram(data = dfTRAIT[dfTRAIT$stand.spec %in% veg_sp,],
                   aes(x = log(PlantHeight), y = ..density..)) +
    xlab("Height (cm)(log)") +
    facet_wrap(~stand.spec, scales = "free", nrow = 7) +
    theme_bw() +
    theme(strip.text = element_text(face = "italic"))
)
ggsave("Figures/Trait_values_vegetation_PH.png", p, height = 8, width = 5)
