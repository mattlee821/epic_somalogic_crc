rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(haven)

# data ====
data <- readRDS("/data/Epic/subprojects/Somalogic/work/Vivian/MultiEndpoints/Data/ForAnalyses/cancer_NormalizedSoma_InvRank_noSMP.rds")
data <- data %>%
  mutate(across(where(~ is.labelled(.)), as_factor))

## exclusions ====
nrow(data) - nrow(subset(data, rowcheck == "PASS"))
data <- subset(data, rowcheck == "PASS") # remove samples that are FLAGGED = 0

nrow(data) - sum(is.na(data$samplenotes))
data <- data[is.na(data$samplenotes), ] # remove samples with sample notes = 61

nrow(data) - sum(is.na(data$assaynotes))
data <- data[is.na(data$assaynotes), ] # remove samples with assay notes = 119

nrow(data) - nrow(subset(data, why_del == ""))
data <- subset(data, why_del == "") # remove samples with a reason for exclusion = 192

# select columns ====
columns_covariates <- c(
  "idepic", "idepic_bio", "idepic_aliq",
  "country", "center",
  "sex",
  "age", # age at recruitment
  "age_rcr", # age at recruitment
  "age_blood", # age at blood collection
  "fasting_c", # fasting status 
  "menop_bld", # menopausal status at blood collection 
  "height_adj", "weight_adj", "hip_adj", "waist_adj", "whr_adj", "bmi_adj",
  "l_school", # highest education
  "smoke_stat", # never, former, smoker, unknown
  "a_1_per_aggr", # age at first period
  "hyster", # hysterectomy
  "ovariect", # ovariectomy
  "menopause", "a_menpoause", # menopause and age at menopause
  "ever_pill", "ever_horm", # ever used pill or hormone
  "use_horm", "use_pill", "use_phrt", # current use of pill or hormone
  "pa_mets", "pa_index", # recreational and household METS and cambridge activity index
  "alc_re", # alcohol at recruitment g/d
  "qe_us_energy", # USDA energy kcal
  "qe_us_fat", # USDA fat grams
  "qe_us_fibtg", # USDA total dietary fiber grams
  "qe_us_procnt", # USDA protein grams
  "hli_score", # healthy lifestyle index
  "c_death_o", # originating cause of death
  "death_status", # dead = 1; alive = 0
  "age_exit_death", # age at exit followuo or death
  "age_exit_cancer_1st", # age at first cancer diagnosis/end of followup
  "ageevent", # age at first cancer diagnosis/end of followup
  "case_t2d", # type 2 diabetes case
  "cvd_t2d_coh" # cvd or t2d event 
)

columns_cancer <- c(
  "indevent", # incident event
  "cancer_1st_sit_tumo", # site of the tumour
  "cancer_1st_mor_tumo", # morphology of the tumour
  "cncr_mal_anyc", # any malignant cancer
  "cncr_mal_clrt", # colorectal cancer
  "cncr_mal_clrt_rectum", # rectal colorectal cancer
  "cncr_mal_clrt_colon" # colon colorectal cancer
)

columns_plate <- c(
  "plateid",
  "scannerid", 
  "plateposition",
  "slideid",
  "subarray",
  "samplenotes",
  "assaynotes",
  "rowcheck"
)

data <- data %>%
  select(any_of(c(columns_covariates, columns_cancer, columns_plate)),
         contains("seq"))

# convert sex: Female=1 and Male=2 ====
data <- data %>%
  mutate(sex = ifelse(sex == "Female", 1, 2))

# subset for cases and non-cases ====
data_cancer <- subset(data, cncr_mal_clrt == "Incident")
data_cancer <- droplevels(data_cancer)
data_control <- subset(data, cncr_mal_anyc == "Non-case")
data_control <- droplevels(data_control)
data_control <- subset(data_control, cancer_1st_sit_tumo == "") # remove all individuals coded with a cancer site
data <- bind_rows(data_cancer, data_control)
data <- droplevels(data)

# subset for subtypes 
data_cancer_colon <- subset(data_cancer, cncr_mal_clrt_colon == "Incident")
data_colon <- bind_rows(data_cancer_colon, data_control)
data_cancer_rectum <- subset(data_cancer, cncr_mal_clrt_rectum == "Incident")
data_rectum <- bind_rows(data_cancer_rectum, data_control)
data_cancer_eo <- subset(data_cancer, ageevent < 55)
data_control_eo <- subset(data_control, age < 55)
data_control_eo <- data_control_eo %>% # change age of exit to 55 for controls
  mutate(age_exit_cancer_1st = ifelse(age_exit_cancer_1st > 55, 55, age_exit_cancer_1st))
data_eo <- bind_rows(data_cancer_eo, data_control_eo)

# save ====
## combined ====
write.table(data, "analysis/001_phenofile/complete/combined/overall/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/combined/overall/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(data_colon, "analysis/001_phenofile/complete/combined/colon/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_colon %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/combined/colon/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(data_rectum, "analysis/001_phenofile/complete/combined/rectum/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_rectum %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/combined/rectum/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(data_eo, "analysis/001_phenofile/complete/combined/early_onset/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_eo %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/combined/early_onset/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## female ====
data_female <- subset(data, sex == 1)
write.table(data_female, "analysis/001_phenofile/complete/female/overall/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_female %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/female/overall/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(data_female, "analysis/001_phenofile/complete/female/colon/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_female %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/female/colon/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(data_female, "analysis/001_phenofile/complete/female/rectum/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_female %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/female/rectum/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(data_female, "analysis/001_phenofile/complete/female/early_onset/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_female %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/female/early_onset/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## male ====
data_male <- subset(data, sex == 2)
write.table(data_male, "analysis/001_phenofile/complete/male/overall/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_male %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/male/overall/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(data_male, "analysis/001_phenofile/complete/male/colon/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_male %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/male/colon/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(data_male, "analysis/001_phenofile/complete/male/rectum/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_male %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/male/rectum/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

write.table(data_male, "analysis/001_phenofile/complete/male/early_onset/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_male %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/male/early_onset/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
