# install.packages(c("tidyverse", "xlsx"))
library(tidyverse)
library(xlsx)

# Accepted pollen types
p_vars <- read.csv("Data/EPD data tables/P_VARS.csv", sep = ";")
# Ecological group names
p_group <- read.csv("Data/EPD data tables/P_GROUP.csv", sep = ";")
# Group codes
group_code <- read.csv("Data/EPD data tables/GROUP_CODES.csv", sep = ";")

# Filter for accepted pollen types names
acc_varname <- p_vars %>% 
  select(AccVar., VarName) %>% 
  distinct(AccVar., .keep_all = TRUE) %>%
  rename(AccVarName = VarName) %>% 
# add originale pollen types names
  left_join(p_vars, by = "AccVar.") %>% 
  select(AccVar. , AccVarName, Var., VarName) 

# Add ecological groups
ecol_group <- p_group %>% 
  left_join(acc_varname, by = "Var.") %>% 
  left_join(group_code, by = "GroupId") %>% 
  # select relevant columns and rename
  select(VarID = Var., GroupID = GroupId, AccVarID = AccVar.,
        AccVarName, VarName, GroupName)

# write to excel
write.xlsx(ecol_group, file = "Data/EPD data tables/EPD_ecological_groups.xlsx",
           row.names = FALSE)
