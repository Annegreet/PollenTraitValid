## ---------------------------
##
## Script name: SI_Remove_taxa
##
## Purpose of script: Testing the effect of a single species or taxon on the CMW
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
if (!require(rjags))  install.packages("rjags")
if (!require(runjags)) install.packages("runjags")
if (!require(mcmc)) install.packages("mcmc")
if (!require(coda)) install.packages("coda")
if (!require(furrr))  install.packages("furrr")

source("Cwm_functions.R")

## Load and prepare data -----
# Trait data ----
dfTRAIT <- readRDS("RDS_files/02_Gapfilled_traits_pollen.rds")

# Pollen data ----
dfPOL_Scot <- readRDS("RDS_files/01_pollen_data_Scot.rds")
dfPOL_Swiss <- readRDS("RDS_files/01_pollen_data_Swiss.rds")
dfPOL <- bind_rows("Scotland" = dfPOL_Scot, "Switzerland" =
                     dfPOL_Swiss, .id = "country") %>% 
  ungroup()


selectedcountry <- "Scotland"
selectedtrait <- "LA"
selectedabun <- str_subset(colnames(dfPOL), "percent")[1]
removedspec <- unique(dfPOL$pollentaxon)[3]
normaltraits <- c("PlantHeight","SLA","LA")

pollentax <- dfPOL_Scot %>% 
  group_by(pollentaxon) %>% 
  filter(any(percent > 0.08)) %>% 
  pull(pollentaxon) %>% 
  unique()

# Pollen  ----
labs_trait <- as_labeller(c(PlantHeight = "Height~(cm)(log)",
                            LA = "Leaf~area~(cm^{2})(log)",
                            SLA ="SLA~(mm^{2}/mg)(log)"),
                          default = label_parsed)

# Run model
plan(multisession(workers = 2))
# furrr::future_map(pollentax,
#                     ~cwm_pol_remove(selectedtrait = "LA",
#                                 selectedcountry = "Scotland",
#                                 selectedabun = "percent",
#                                 removedspec = .x),
#                   .options = furrr_options(seed = TRUE))
furrr::future_map(pollentax,
                  ~cwm_pol_remove(selectedtrait = "SLA",
                              selectedcountry = "Scotland",
                              selectedabun = "percent",
                              removedspec = .x),
                  .options = furrr_options(seed = TRUE))
furrr::future_map(pollentax,
                  ~cwm_pol_remove(selectedtrait = "PlantHeight",
                              selectedcountry = "Scotland",
                              selectedabun = "percent",
                              removedspec = .x),
                  .options = furrr_options(seed = TRUE))  


if(1){
  # Compile estimates
pollen_files <- 
  list.files("RDS_files/",full.names = TRUE) %>% 
  str_subset(pattern = "pollen") %>% 
  str_subset(pattern = paste(paste(unique(dfPOL_Scot$pollentaxon), collapse = "|"), "Ericales", sep = "|")) 

pollen_files_all <-  list.files("RDS_files/",full.names = TRUE) %>% 
  str_subset(pattern = "03_CWM_estimates_pollen_Scotland_LA_percent.rds|03_CWM_estimates_pollen_Scotland_SLA_percent.rds|03_CWM_estimates_pollen_Scotland_PlantHeight_percent.rds")

dfCWM_pol <- c(pollen_files, pollen_files_all) %>% 
  purrr::map2(., str_remove(c(pollen_files, pollen_files_all), "RDS_files/"), ~readRDS(.x) %>% 
                as.data.frame())  %>% 
  bind_rows()  %>% 
  mutate(removed_spec = ifelse(is.na(removed_spec), "All taxa", removed_spec))

labs_trait <- as_labeller(c(PlantHeight = "Plant~height~(cm)(log)",
                            LA = "Leaf~area~(cm^{2})(log)",
                            SLA ="SLA~(mm^{2}/mg)(log)"),
                          default = label_parsed)

pollen_lab <- as_labeller(c("Betula" = "italic(Betula)",
                            "Cyperaceae" = "Cyperaceae",
                            "Ericales (tetrad)" = "Ericales~(tetrad)",
                            "Pinus" = "italic(Pinus)",
                            "Poaceae" = "Poaceae",
                            "Pteridophyte" = "Pteridophyte"), default = label_parsed)

(p <- dfCWM_pol %>% 
  filter(removed_spec %in% c(pollentax,
                             "All taxa")) %>%
  ggplot(aes(x = Mean, y = removed_spec)) +
  geom_boxplot(fill = "grey") +
  scale_y_discrete("", labels = expression("Betula" = italic(Betula),
                                  "Cyperaceae" = "Cyperaceae",
                                  "Ericales (tetrad)" = "Ericales (tetrad)",
                                  "Pinus" = italic(Pinus),
                                  "Poaceae" = "Poaceae",
                                  "Pteridophyte" = "Pteridophyte",
                                  "All taxa" = "All taxa")) +
  xlab("") +
  facet_grid(~trait, labeller = labs_trait, scales = "free") +
  theme_bw() +
  theme(text = element_text(size = 10),
        strip.background = element_rect(fill = "white"),
        legend.position = "none")
)
ggsave("Figures/CWM_pollen_removed_taxa.png", height = 4, width = 7)
}
