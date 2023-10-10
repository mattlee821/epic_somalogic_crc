rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(haven)

# data ====
data <- readRDS("data/Data/ForAnalysis_cancer.rds")
data <- data %>%
  mutate(across(where(~ is.labelled(.)), as_factor))

## exclusions ====
data <- subset(data, rowcheck == "PASS") # remove samples that are FLAGGED = 112
data <- data[is.na(data$samplenotes), ] # remove samples with sample notes = 36
data <- data[is.na(data$assaynotes), ] # remove samples with assay notes = 113
data <- subset(data, why_del == "") # remove samples with a reason for exclusion = 179

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

# subset for colorectal cancer and non-cases ====
data_cancer <- subset(data, cncr_mal_clrt == "Incident")
data_cancer <- droplevels(data_cancer)
table(data_cancer$cancer_1st_sit_tumo)

data_control <- subset(data, cncr_mal_anyc == "Non-case")
data_control <- droplevels(data_control)
table(data_control$cancer_1st_sit_tumo)
data_control <- subset(data_control, cancer_1st_sit_tumo == "") # remove all individuals coded with a cancer site

data <- bind_rows(data_cancer, data_control)
data <- droplevels(data)

# convert sex: Female=1 and Male=2 ====
data <- data %>%
  mutate(sex = ifelse(sex == "Female", 1, 2))

# save ====
## combined ====
write.table(data, "analysis/001_phenofile/complete/combined/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/combined/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_cancer <- subset(data, cncr_mal_clrt == "Incident")
write.table(data_cancer, "analysis/001_phenofile/cancer/combined/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_cancer %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/cancer/combined/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_control <- subset(data, cncr_mal_anyc == "Non-case")
write.table(data_control, "analysis/001_phenofile/non-cancer/combined/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_control %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/non-cancer/combined/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## female ====
data_female <- subset(data, sex == 1)
write.table(data_female, "analysis/001_phenofile/complete/female/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_female %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/female/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_cancer <- subset(data_female, cncr_mal_clrt == "Incident")
write.table(data_cancer, "analysis/001_phenofile/cancer/female/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_cancer %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/cancer/female/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_control <- subset(data_female, cncr_mal_anyc == "Non-case")
write.table(data_control, "analysis/001_phenofile/non-cancer/female/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_control %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/non-cancer/female/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

## male ====
data_male <- subset(data, sex == 2)
write.table(data_male, "analysis/001_phenofile/complete/male/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_male %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/male/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_cancer <- subset(data_male, cncr_mal_clrt == "Incident")
write.table(data_cancer, "analysis/001_phenofile/cancer/male/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_cancer %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/cancer/male/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_control <- subset(data_male, cncr_mal_anyc == "Non-case")
write.table(data_control, "analysis/001_phenofile/non-cancer/male/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_control %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/non-cancer/male/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
