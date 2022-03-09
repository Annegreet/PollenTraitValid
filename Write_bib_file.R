library(bibtex)
library(tidyverse)

# Compile .bib file multiple references
ref_files <- list.files("References/")

ref_list <- list()

for(i in 1:length(ref_files)){
  ref <- read.bib(paste0("References/", ref_files[i]))
  ref_list <- append(ref_list, ref)
}

write.bib(ref_list, file ="Validation_study_references")
