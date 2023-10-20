rm(list=ls())
set.seed(821)

# enviroment ====
library(dplyr)
library(data.table)
library(SomaDataIO)

# data ====
info_analyte <- read_adat("/data/Epic/subprojects/Depot_Somalogic/sources/SS-2230719.SS-2215915.SS-2230718.SS-2230720_v4.1_CitratePlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
info_analyte <- getAnalyteInfo(info_analyte)
info_analyte <- select(info_analyte, AptName, SeqId, SeqIdVersion, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, Organism, Units, Type, Dilution)
file_list <- list.files("analysis/003_analysis/complete/", pattern = "results", full.names = TRUE, recursive = TRUE)
table <- data.frame()

for (files in file_list) {
  cat("Processing file:", files, "\n")
  # make labels ====  
  path_in <- gsub("results.txt", "", files)
  path_label <- unlist(strsplit(path_in, "/"))
  sex <- path_label[length(path_label) - 1]
  cancer <- path_label[length(path_label)]

  # read data ====
  data <- fread(files, header = TRUE)
  data$sex <- sex
  data$cancer <- cancer
  data <- left_join(data, info_analyte, by = c("exposure" = "AptName"))
  table <- bind_rows(table, data)
}

# save ====
write.table(table, "analysis/003_analysis/complete/results_formatted.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

