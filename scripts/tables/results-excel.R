rm(list=ls())
set.seed(821)

# enviroment ====
library(dplyr)
library(data.table)
library(tidyr)
library(openxlsx)

# data ====
data <- fread("analysis/003_analysis/complete/results_formatted.txt")
## rearrange 
data <- select(data,
               "Target", "SeqId", "cancer",
               "sex", "analysis", "followup", "model", "covariates",
               "n", "nevent",
               "coef_exp", "se", "ci_lower_exp", "ci_upper_exp", "pval", "FDR_BH", "FDR_bonferroni")

# calculate change in HR for models ====
data_change <- data %>%
  mutate(ID = paste(Target, SeqId, cancer, sex, analysis, followup, sep = ";")) %>%
  select(ID, model, coef_exp) %>%
  pivot_wider(names_from = model, values_from = coef_exp)

data_change <- data_change %>%
  mutate(
    percent_change_1_2 = ((`1` - `2`) / `1`) * 100,
    percent_change_1_3 = ((`1` - `3`) / `1`) * 100,
    percent_change_1_4 = ((`1` - `4`) / `1`) * 100,
    percent_change_1_5 = ((`1` - `5`) / `1`) * 100
  )

# make dataframes ====
## Get unique combinations of levels in the specified columns
combinations <- unique(data %>%
                         select(cancer, sex, analysis, followup))
for (col in names(combinations)) {
  combinations[[col]] <- as.factor(combinations[[col]])
}

combinations$cancer <- factor(combinations$cancer, levels = c("overall", "colon", "rectum", "early_onset"))
combinations$sex <- factor(combinations$sex, levels = c("combined", "female", "male"))
combinations$analysis <- factor(combinations$analysis, levels = c("complete", "excl2year"))
combinations$followup <- factor(combinations$followup, levels = c("complete", "below", "above"))
setorder(combinations, cancer, sex, analysis, followup)
combinations <- paste0(combinations$cancer, ".", combinations$sex, ".", combinations$analysis, ".", combinations$followup)

data_frames <- split(data, list(data$cancer, data$sex, data$analysis, data$followup))
data_frames <- data_frames[match(combinations, names(data_frames))]

# rename ====
data_list <- c(list(data_change), data_frames)
names(data_list)[1] <- "percent_change" # rename first dataframe
## rename dataframes for string limit in excel
for (name in names(data_list)) { 
  new_name <- gsub("overall", "all", name)
  names(data_list)[names(data_list) == name] <- new_name
}
for (name in names(data_list)) { 
  new_name <- gsub("early_onset", "eo", name)
  names(data_list)[names(data_list) == name] <- new_name
}
for (name in names(data_list)) { 
  new_name <- gsub("combined", "sc", name)
  names(data_list)[names(data_list) == name] <- new_name
}
for (name in names(data_list)) { 
  new_name <- gsub("excl2year", "2year", name)
  names(data_list)[names(data_list) == name] <- new_name
}

# save ====
## all ====
wb <- createWorkbook()
for (i in seq_along(data_list)) {
  addWorksheet(wb, sheetName = names(data_list)[i])
  writeData(wb, sheet = i, x = data_list[[i]])
}
saveWorkbook(wb, file = "analysis/tables/results.xlsx")

## model 1 ====
data_list1 <- lapply(seq_along(data_list)[-1], function(i) {
  subset(data_list[[i]], model == 1)
})
names(data_list1) <- names(data_list)[-1] # Assign names back to the modified dataframes
wb <- createWorkbook()
for (i in seq_along(data_list1)) {
  addWorksheet(wb, sheetName = names(data_list1)[i])
  writeData(wb, sheet = i, x = data_list1[[i]])
}
saveWorkbook(wb, file = "analysis/tables/results-model1.xlsx")

## model 2 ====
data_list1 <- lapply(seq_along(data_list)[-1], function(i) {
  subset(data_list[[i]], model == 2)
})
names(data_list1) <- names(data_list)[-1] # Assign names back to the modified dataframes
wb <- createWorkbook()
for (i in seq_along(data_list1)) {
  addWorksheet(wb, sheetName = names(data_list1)[i])
  writeData(wb, sheet = i, x = data_list1[[i]])
}
saveWorkbook(wb, file = "analysis/tables/results-model2.xlsx")

## model 3 ====
data_list1 <- lapply(seq_along(data_list)[-1], function(i) {
  subset(data_list[[i]], model == 3)
})
names(data_list1) <- names(data_list)[-1] # Assign names back to the modified dataframes
wb <- createWorkbook()
for (i in seq_along(data_list1)) {
  addWorksheet(wb, sheetName = names(data_list1)[i])
  writeData(wb, sheet = i, x = data_list1[[i]])
}
saveWorkbook(wb, file = "analysis/tables/results-model3.xlsx")

## model 4 ====
data_list1 <- lapply(seq_along(data_list)[-1], function(i) {
  subset(data_list[[i]], model == 4)
})
names(data_list1) <- names(data_list)[-1] # Assign names back to the modified dataframes
wb <- createWorkbook()
for (i in seq_along(data_list1)) {
  addWorksheet(wb, sheetName = names(data_list1)[i])
  writeData(wb, sheet = i, x = data_list1[[i]])
}
saveWorkbook(wb, file = "analysis/tables/results-model4.xlsx")

## model 5 ====
data_list1 <- lapply(seq_along(data_list)[-1], function(i) {
  subset(data_list[[i]], model == 5)
})
names(data_list1) <- names(data_list)[-1] # Assign names back to the modified dataframes
wb <- createWorkbook()
for (i in seq_along(data_list1)) {
  addWorksheet(wb, sheetName = names(data_list1)[i])
  writeData(wb, sheet = i, x = data_list1[[i]])
}
saveWorkbook(wb, file = "analysis/tables/results-model5.xlsx")

## specific table ====
data1 <- subset(data, model == 5)
data1 <- subset(data1, analysis == "complete")
data1 <- subset(data1, followup == "complete")
data_frames <- split(data1, list(data1$cancer, data1$sex))
wb <- createWorkbook()
for (i in seq_along(data_frames)) {
  addWorksheet(wb, sheetName = names(data_frames)[i])
  writeData(wb, sheet = i, x = data_frames[[i]])
}
saveWorkbook(wb, file = "analysis/tables/results-specific.xlsx")


