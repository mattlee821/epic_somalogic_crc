rm(list = ls())

# enviroment ====
library(dplyr)
library(data.table)

# complete ====
data <- fread("analysis/002_metaboprep/complete/metaboprep_release_2023_09_19/filtered_data/epic_proteomics_somalogic_complete_2023_09_19_Filtered_metabolite_data.txt")
colnames(data)[1] <- "idepic"
id_data <- data$idepic
phenofile <- fread("analysis/001_phenofile/complete/phenofile.txt")
phenofile <- phenofile[phenofile$idepic %in% id_data, ]
phenofile <- phenofile %>%
  select(-contains("seq"))
data <- left_join(phenofile, data, by = "idepic")
##
write.table(data, "analysis/003_data-for-analysis/complete/data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# cancer ====
data <- fread("analysis/002_metaboprep/cancer/metaboprep_release_2023_09_21/filtered_data/epic_proteomics_somalogic_complete_2023_09_21_Filtered_metabolite_data.txt")
colnames(data)[1] <- "idepic"
id_data <- data$idepic
phenofile <- fread("analysis/001_phenofile/cancer/phenofile.txt")
phenofile <- phenofile[phenofile$idepic %in% id_data, ]
phenofile <- phenofile %>%
  select(-contains("seq"))
data <- left_join(phenofile, data, by = "idepic")
##
write.table(data, "analysis/003_data-for-analysis/cancer/data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# non-cancer ====
data <- fread("analysis/002_metaboprep/non-cancer/metaboprep_release_2023_09_21/filtered_data/epic_proteomics_somalogic_complete_2023_09_21_Filtered_metabolite_data.txt")
colnames(data)[1] <- "idepic"
id_data <- data$idepic
phenofile <- fread("analysis/001_phenofile/non-cancer/phenofile.txt")
phenofile <- phenofile[phenofile$idepic %in% id_data, ]
phenofile <- phenofile %>%
  select(-contains("seq"))
data <- left_join(phenofile, data, by = "idepic")
##
write.table(data, "analysis/003_data-for-analysis/non-cancer/data.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
