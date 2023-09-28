rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(haven)

# data ====
data <- readRDS("data/Data/ForAnalysis_cancer.rds")
data <- data %>%
  mutate(across(where(~ is.labelled(.)), as_factor))

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
table(data$cancer_1st_mor_tumo)
data <- data %>%
  mutate(sex = ifelse(sex == "Female", 1, 2))

# save ====
write.table(data, "analysis/001_phenofile/complete/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/complete/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_cancer <- subset(data, cncr_mal_clrt == "Incident")
write.table(data_cancer, "analysis/001_phenofile/cancer/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_cancer %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/cancer/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

data_control <- subset(data, cncr_mal_anyc == "Non-case")
write.table(data_control, "analysis/001_phenofile/non-cancer/phenofile.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
data_proteins <- data_control %>%
  select(idepic,
         contains("seq"))
write.table(data_proteins, "analysis/001_phenofile/non-cancer/proteins.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
