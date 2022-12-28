# This script writes a reference list from the bib files in the reference
# folder

library(tidyverse)
library(bibtex)
library(rbibutils)

# Compile .bib file multiple references
ref_files <- list.files("References/", full.names = TRUE)

ref_list <- list()

for (i in 1:length(ref_files)) {
  ref <- readBib(ref_files[i], direct = TRUE, texChars = "export")
  ref_list <- append(ref_list, ref)
}

write.bib(ref_list, file = "Validation_study_references.bib")
