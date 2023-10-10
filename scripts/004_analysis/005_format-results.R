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

data <- fread("analysis/004_analysis/results.txt")
data <- left_join(data, info_analyte, by = c("exposure" = "AptName"))

