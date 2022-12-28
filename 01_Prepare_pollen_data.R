## ---------------------------
##
## Script name: 
##
## Purpose of script:
##
## Author: Annegreet Veeken
##
## Date Created: 2022-01
##
## Email: veeken.g.a@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

options(scipen = 999)
set.seed(1)
## Load packages
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(readxl)) install.packages("readxl")


# Scotland ----
pollen_raw <- read.csv("Data/Pollen_count_Scotland-checked.csv", skip = 1, sep = ";")
ppe <- read_xlsx("Data/Githumbi_2022_RPP.xlsx")
polmode <- readRDS("RDS_files/01_Pollination_mode.rds")

# make tidy
dfPOL <- pollen_raw %>% 
  # filter empty rows
  filter(!Key.index == "---") %>% 
  filter(!str_detect(Sample.name,"---")) %>% 
  # dplyr relevant columns
  dplyr::select(pollentaxon = Sample.name, contains("C0")) %>% 
  # convert to long format
  pivot_longer(cols = contains("C0"), names_to = "sitename",
               values_to = "count") %>% 
  # correct miss named site name in polycounter
  filter(!sitename == "C014") %>% 
  mutate(sitename = recode(sitename, C014.1 = "C014"),
         # convert count to numeric       
         count = as.numeric(count)) %>% 
  # remove zeros
  filter(!count == 0) %>% 
  mutate(pollentaxon = recode(pollentaxon, Tsuga = "unindentified", 
                              Liliaceae = "unindentified", Drosera = "unindentified")) %>%  # probably misidentified
  #correct pollen taxon names
  mutate(pollentaxon = 
           # remove "-type" to facilitate GBIF search
           str_remove(pollentaxon, "-type") %>% 
           # Typos
           recode(Caryophyllaccea = "Caryophyllaceae",
                  "Juniperus communis" = "Juniperus")
         ) %>% 
  # filter out unidentified
  filter(!pollentaxon %in% c("unindentified")) %>% 
  filter(!is.na(count))

# percentage data with missing ppe's 
missing_ppe_scot <- dfPOL %>%
  # join with ppe data 
  left_join(ppe, by = c("pollentaxon" = "joining_taxon")) %>% 
  group_by(sitename) %>% 
  mutate(percent = count/sum(count, na.rm = TRUE)) %>% 
  filter(is.na(Taxon)) %>% 
  summarise(percent_missing = sum(percent))

# calculate pollen percentage - not corrected
dfPOL_nc <- dfPOL %>%
  # calculate percentages
  group_by(sitename) %>% 
  mutate(percent = count/sum(count, na.rm = TRUE)) %>%
  group_by(sitename, pollentaxon) %>%
  summarise(percent = sum(percent)) 

# make function to continue draw from rnorm till RPP is positive (otherwise negative percentages)
norm_positive <- function(RPP, RPP_sd){
  success <- FALSE
  while (!success) {
    x <- rnorm(1, RPP, RPP_sd)
    success <- x > 0
    }
  return(x)
}

n <- 100 # number of draws from RPP distribution
dfPOL_cor <- dfPOL %>%
  # join with ppe data 
  left_join(ppe, by = c("pollentaxon" = "joining_taxon")) %>% 
  # remove taxa without correction factors
  filter(!is.na(Taxon)) %>% 
  # correction with mean RPP
  mutate(mean = count/RPP) %>%  
  # correction with random draw from RPP distribution
  rowwise() %>%
  mutate(col_rnorm = list(setNames(count/replicate(n,norm_positive(RPP, RPP_sd)), paste0("draw",1:n)))) %>%
  unnest_wider(col_rnorm) %>% 
  # calculate percentages
  group_by(sitename) %>% 
  mutate(across(contains(c("mean", "draw")), ~./sum(.))) %>%
  group_by(sitename, pollentaxon) %>%
  summarise(across(contains(c("mean", "draw")), sum,
            .names = "adjustedpercent_{col}")) 

dfPOL <- dfPOL_nc %>% 
  left_join(dfPOL_cor, by = c("sitename", "pollentaxon")) %>% 
  # add PFT and pollination mode
  left_join(polmode, by = "pollentaxon") %>% 
  ungroup() 

saveRDS(dfPOL, "RDS_files/01_Pollen_data_Scot.rds")

# Switzerland ----
# read in pollen data
dfPOL <- read.csv("Data/Pollen_count_Switserland-checked.csv", skip = 1, sep = ";") %>% 
  filter(!Key.index == "---") %>% # filter empty rows
  filter(!Sample.name == "---") %>% 
  pivot_longer(cols = c(starts_with("X"),"unknown"), # convert to long format
               names_to = "sitename",
               values_to = "count") %>%
  mutate(count = as.numeric(count)) %>% # change variable types
  dplyr::select(sitename, pollentaxon = Sample.name, count) %>%  # select relevant rows
  mutate(pollentaxon = recode(pollentaxon, Tsuga = "unindentified", 
         Liliaceae = "unindentified", Drosera = "unindentified")) %>%  # probably misidentified
  # remove zeros
  filter(!count == 0) %>% 
  # filter out unidentified
  filter(!pollentaxon %in% c("unindentified")) %>% 
  #correct pollentaxon names
  mutate(pollentaxon = 
           # remove "-type" to facilitate GBIF search
           str_remove(pollentaxon, "-type") %>% 
           # Typos
           recode(Caryophyllaccea = "Caryophyllaceae")
  ) %>% 
  mutate(sitename = recode(sitename, "X563186" = "X563166" )) # miss named site

# percentage data with missing ppe's 
missing_ppe_swiss <- dfPOL %>%
  # join with ppe data 
  left_join(ppe, by = c("pollentaxon" = "joining_taxon")) %>% 
  group_by(sitename) %>% 
  mutate(percent = count/sum(count, na.rm = TRUE)) %>% 
  filter(is.na(Taxon)) %>% 
  summarise(percent_missing = sum(percent))

# calculate pollen percentage - not corrected
dfPOL_nc <- dfPOL %>%
  # calculate percentages
  group_by(sitename) %>% 
  mutate(percent = count/sum(count, na.rm = TRUE)) %>%
  group_by(sitename, pollentaxon) %>%
  summarise(percent = sum(percent)) 

dfPOL_cor <- dfPOL %>%
  # join with ppe data 
  left_join(ppe, by = c("pollentaxon" = "joining_taxon")) %>% 
  # remove taxa without correction factors
  filter(!is.na(Taxon)) %>% 
  # correction with mean RPP
  mutate(mean = count/RPP) %>%  
  # correction with random draw from RPP distribution
  rowwise() %>%
  mutate(col_rnorm = list(setNames(count/replicate(n,norm_positive(RPP, RPP_sd)), paste0("draw",1:n)))) %>%
  unnest_wider(col_rnorm) %>% 
  # calculate percentages
  group_by(sitename) %>% 
  mutate(across(contains(c("mean", "draw")), ~./sum(.))) %>%
  group_by(sitename, pollentaxon) %>%
  summarise(across(contains(c("mean", "draw")), sum,
                   .names = "adjustedpercent_{col}")) 

dfPOL <- dfPOL_nc %>% 
  left_join(dfPOL_cor, by = c("sitename", "pollentaxon")) %>% 
  # add PFT and pollination mode
  left_join(polmode, by = "pollentaxon") %>% 
  ungroup()

saveRDS(dfPOL, "RDS_files/01_Pollen_data_Swiss.rds")
