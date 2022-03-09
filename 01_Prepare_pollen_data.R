library(tidyverse)
library(readxl)

# Scotland ----
pollen_raw <- read.csv("Data/Pollen_count.csv", skip = 1)
ppe <- read_xlsx("Data/RPP_Dataset_v1_Table_3_4.xlsx")
ppe <- tibble(ppe.taxon = ppe$`Target taxon`,
              PPE = ppe$`RPP v1 (Northern Hemisphere)`)


# make tidy
dfPOL <- pollen_raw %>% 
  # filter empty rows
  filter(!Key.index == "---") %>% 
  filter(!str_detect(Sample.name,"---")) %>% 
  # select relevant columns
  select(pollentaxon = Sample.name, contains("C0")) %>% 
  # convert to long format
  pivot_longer(cols = contains("C0"), names_to = "sitename",
               values_to = "count") %>% 
  # correct miss named site name in polycounter
  filter(!sitename == "C014") %>% 
  mutate(sitename = recode(sitename, C014.1 = "C014"),
         # convert count to numeric       
         count = as.numeric(count)) %>% 
  #correct pollentaxon names
  mutate(pollentaxon = 
           # remove "-type" to facilitate GBIF search
           str_remove(pollentaxon, "-type") %>% 
           # Typos
           recode(Caryophyllaccea = "Caryophyllaceae",
                  Lilaceae = "Liliaceae")
         )
# join with ppe data
dfPOL <- dfPOL %>% 
  mutate(ppe.taxon = case_when(pollentaxon %in% ppe$ppe.taxon ~ pollentaxon,
                               pollentaxon == "Juniperus communis" ~ "Juniperus",
                               pollentaxon == "Ericales (tetrad)" ~ "Ericales",
                               pollentaxon == "Carpinus betulus" ~ "Carpinus",
                               pollentaxon == "Plantago" ~ "Plantaginaceae",
                               pollentaxon %in% c("Pteridophyte", "Geraniaceae", "Viola","Pedicularis palustris",
                                                  "Malvaceae","Tsuga" ) ~ "Not applicable"))
# calculate pollen percentage
dfPOL <- dfPOL %>%
  filter(!pollentaxon == "unindentified") %>% 
  left_join(ppe, by = "ppe.taxon") %>% 
  replace_na(list(PPE = 1)) %>% # if no ppe data is available apply no correction?
  group_by(sitename) %>% 
  mutate(adjustedcount = count*PPE) %>%
  # calculate percentages
  mutate(adjustedpercent = adjustedcount/sum(adjustedcount),
         percent = count/sum(count, na.rm = TRUE)) %>%
  group_by(sitename, pollentaxon) %>%
  summarise(adjustedpercent = sum(adjustedpercent, na.rm = TRUE),
            percent = sum(percent)) 
saveRDS(dfPOL, "RDS_files/01_Pollen_data_Scot.rds")

# Switzerland ----
# read in pollen data
pollen <- read.csv("Pollen_count_Switserland.csv", skip = 1) %>% 
  filter(!Key.index == "---") %>% # filter empty rows
  filter(!Sample.name == "---") %>% 
  pivot_longer(cols = c(starts_with("X"),"unknown"), # convert to long format
               names_to = "site.name",
               values_to = "count") %>%
  mutate(count = as.numeric(count)) %>% # change variable types
  select(site.name, taxa = Sample.name, count)  # select relevant rows

# pollen sum
pollen %>%   
  group_by(site.name) %>% 
  summarise(polsum = sum(count)) %>%
  filter(polsum < 600) %>% 
  arrange(polsum) %>% 
  print(n = nrow(.))

# pollen percentage
pollen %>% 
  group_by(site.name) %>% 
  mutate(polsum = sum(count)) %>% 
  group_by(site.name, taxa) %>% 
  mutate(prop = count / polsum * 100) %>% 
  View()

# How many pollen left to count
pollen %>%   
  group_by(site.name) %>% 
  summarise(polsum = sum(count),
            polleft = 600 - polsum) %>%
  ungroup %>% 
  summarise(sum(polleft))




