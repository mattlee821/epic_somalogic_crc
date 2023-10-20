rm(list=ls())
set.seed(821)

# enviroment ====
library(dplyr)
library(data.table)
library(survival)

# data ====
file_list <- list.files("analysis/002_outlier-exclusions/complete/female/colon/", pattern = "phenofile", full.names = TRUE, recursive = TRUE)

# Loop through each file
for (files in file_list) {
  cat("Processing file:", files, "\n")
  
  # make labels ====  
  path_in <- gsub("phenofile.txt", "", files)
  path_out <- gsub("002_outlier-exclusions", "003_analysis", path_in)
  path_label <- unlist(strsplit(path_in, "//"))
  path_label <- gsub("/", "-", path_label[length(path_label)])
  path_label <- sub("-*$", "", path_label)
  
  # read data ====
  data <- fread(files, header = TRUE)
  proteins <- names(data)[grepl("seq", names(data))]
  data$followup <- (data$age_exit_cancer_1st - data$age)
  median_follow_up <- median(data$age_exit_cancer_1st - data$age)
  
  ## format data ====
  data <- data %>%  mutate(cvd_t2d_coh = ifelse(cvd_t2d_coh == "No", 0, 1))
  data <- data %>% mutate(age_group = as.integer(factor(cut(age, breaks = seq(0, max(age) + 5, by = 5), right = FALSE, labels = seq(1, max(age), by = 5)))))
  data_median_above <- subset(data, followup >= median_follow_up)
  data_median_below <- subset(data, followup >= median_follow_up)
  
  # prentic weighting ====
  
  ## full data: prentice weighting ====
  CasesOutsideSubCohort  <- data %>% filter(cvd_t2d_coh == 0 & indevent == 1) %>% mutate(age = ageevent-1e-4)
  CasesInSubcohort_Ctrls <- data %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(ageevent = ageevent-1e-4, indevent = 0)
  CasesInSubcohort_Cases <- data %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(age = ageevent-1e-4, cvd_t2d_coh = 0) # so that they are not included in the computation of the sd
  ControlsSubcohort      <- data %>% filter(cvd_t2d_coh == 1 & indevent == 0)
  data_analysis <- bind_rows(CasesOutsideSubCohort, CasesInSubcohort_Ctrls, CasesInSubcohort_Cases, ControlsSubcohort) # by construction, participants who are not cases and not in the subcohort are excluded (useful - "required" even, for the analysis of site-specific cancers; for the other events, they were excluded in the construction of ForAnalysis...)
  
  ## 2 year followup data ====
  data_excl2year <- data %>% filter(ageevent > age + 2) %>% mutate(age = age + 2) # exclude diagnosis within 2 years
  CasesOutsideSubCohort  <- data_excl2year %>% filter(cvd_t2d_coh == 0 & indevent == 1) %>% mutate(age = ageevent-1e-4)
  CasesInSubcohort_Ctrls <- data_excl2year %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(ageevent = ageevent-1e-4, indevent = 0)
  CasesInSubcohort_Cases <- data_excl2year %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(age = ageevent-1e-4, cvd_t2d_coh = 0) # so that they are not included in the computation of the sd
  ControlsSubcohort      <- data_excl2year %>% filter(cvd_t2d_coh == 1 & indevent == 0)
  data_analysis_excl2year <- bind_rows(CasesOutsideSubCohort, CasesInSubcohort_Ctrls, CasesInSubcohort_Cases, ControlsSubcohort) # by construction, participants who are not cases and not in the subcohort are excluded (useful - "required" even, for the analysis of site-specific cancers; for the other events, they were excluded in the construction of ForAnalysis...)
  
  ## median follow-up time: above ====
  CasesOutsideSubCohort  <- data_median_above %>% filter(cvd_t2d_coh == 0 & indevent == 1) %>% mutate(age = ageevent-1e-4)
  CasesInSubcohort_Ctrls <- data_median_above %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(ageevent = ageevent-1e-4, indevent = 0)
  CasesInSubcohort_Cases <- data_median_above %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(age = ageevent-1e-4, cvd_t2d_coh = 0) # so that they are not included in the computation of the sd
  ControlsSubcohort      <- data_median_above %>% filter(cvd_t2d_coh == 1 & indevent == 0)
  data_analysis_median_above <- bind_rows(CasesOutsideSubCohort, CasesInSubcohort_Ctrls, CasesInSubcohort_Cases, ControlsSubcohort) # by construction, participants who are not cases and not in the subcohort are excluded (useful - "required" even, for the analysis of site-specific cancers; for the other events, they were excluded in the construction of ForAnalysis...)
  
  ### 2 year followup
  data_excl2year_median_above <- data_median_above %>% filter(ageevent > age + 2) %>% mutate(age = age + 2) # exclude diagnosis within 2 years
  CasesOutsideSubCohort  <- data_excl2year_median_above %>% filter(cvd_t2d_coh == 0 & indevent == 1) %>% mutate(age = ageevent-1e-4)
  CasesInSubcohort_Ctrls <- data_excl2year_median_above %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(ageevent = ageevent-1e-4, indevent = 0)
  CasesInSubcohort_Cases <- data_excl2year_median_above %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(age = ageevent-1e-4, cvd_t2d_coh = 0) # so that they are not included in the computation of the sd
  ControlsSubcohort      <- data_excl2year_median_above %>% filter(cvd_t2d_coh == 1 & indevent == 0)
  data_analysis_median_above_excl2year <- bind_rows(CasesOutsideSubCohort, CasesInSubcohort_Ctrls, CasesInSubcohort_Cases, ControlsSubcohort) # by construction, participants who are not cases and not in the subcohort are excluded (useful - "required" even, for the analysis of site-specific cancers; for the other events, they were excluded in the construction of ForAnalysis...)
  
  ## median follow-up time: below ====
  CasesOutsideSubCohort  <- data_median_below %>% filter(cvd_t2d_coh == 0 & indevent == 1) %>% mutate(age = ageevent-1e-4)
  CasesInSubcohort_Ctrls <- data_median_below %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(ageevent = ageevent-1e-4, indevent = 0)
  CasesInSubcohort_Cases <- data_median_below %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(age = ageevent-1e-4, cvd_t2d_coh = 0) # so that they are not included in the computation of the sd
  ControlsSubcohort      <- data_median_below %>% filter(cvd_t2d_coh == 1 & indevent == 0)
  data_analysis_median_below <- bind_rows(CasesOutsideSubCohort, CasesInSubcohort_Ctrls, CasesInSubcohort_Cases, ControlsSubcohort) # by construction, participants who are not cases and not in the subcohort are excluded (useful - "required" even, for the analysis of site-specific cancers; for the other events, they were excluded in the construction of ForAnalysis...)
  
  ### 2 year followup
  data_excl2year_median_below <- data_median_below %>% filter(ageevent > age + 2) %>% mutate(age = age + 2) # exclude diagnosis within 2 years
  CasesOutsideSubCohort  <- data_excl2year_median_below %>% filter(cvd_t2d_coh == 0 & indevent == 1) %>% mutate(age = ageevent-1e-4)
  CasesInSubcohort_Ctrls <- data_excl2year_median_below %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(ageevent = ageevent-1e-4, indevent = 0)
  CasesInSubcohort_Cases <- data_excl2year_median_below %>% filter(cvd_t2d_coh == 1 & indevent == 1) %>% mutate(age = ageevent-1e-4, cvd_t2d_coh = 0) # so that they are not included in the computation of the sd
  ControlsSubcohort      <- data_excl2year_median_below %>% filter(cvd_t2d_coh == 1 & indevent == 0)
  data_analysis_median_below_excl2year <- bind_rows(CasesOutsideSubCohort, CasesInSubcohort_Ctrls, CasesInSubcohort_Cases, ControlsSubcohort) # by construction, participants who are not cases and not in the subcohort are excluded (useful - "required" even, for the analysis of site-specific cancers; for the other events, they were excluded in the construction of ForAnalysis...)
  
  ## clean up ====
  rm(list=setdiff(ls(), c("file_list", "files", "data_analysis", "data_analysis_excl2year", "data_analysis_median_above", "data_analysis_median_above_excl2year", "data_analysis_median_below", "data_analysis_median_below_excl2year", "proteins", "path_in", "path_out", "path_label")))
  
  # analysis ====
  ## make empty data frames for results ====
  table_model1 <- data.frame()
  table_model1_excl2year <- data.frame()
  table_model1_median_above <- data.frame()
  table_model1_median_above_excl2year <- data.frame()
  table_model1_median_below <- data.frame()
  table_model1_median_below_excl2year <- data.frame()
  table_model2 <- data.frame()
  table_model2_excl2year <- data.frame()
  table_model2_median_above <- data.frame()
  table_model2_median_above_excl2year <- data.frame()
  table_model2_median_below <- data.frame()
  table_model2_median_below_excl2year <- data.frame()
  table_model3 <- data.frame()
  table_model3_excl2year <- data.frame()
  table_model3_median_above <- data.frame()
  table_model3_median_above_excl2year <- data.frame()
  table_model3_median_below <- data.frame()
  table_model3_median_below_excl2year <- data.frame()
  table_model4 <- data.frame()
  table_model4_excl2year <- data.frame()
  table_model4_median_above <- data.frame()
  table_model4_median_above_excl2year <- data.frame()
  table_model4_median_below <- data.frame()
  table_model4_median_below_excl2year <- data.frame()
  table_model5 <- data.frame()
  table_model5_excl2year <- data.frame()
  table_model5_median_above <- data.frame()
  table_model5_median_above_excl2year <- data.frame()
  table_model5_median_below <- data.frame()
  table_model5_median_below_excl2year <- data.frame()
  
  ## analysis ====
  ksq <- 0
  for (exposure in proteins)
  {
    ksq <- ksq +1
    if(ksq%%100 ==1){cat("+++++++++", ksq, "\n")}
    
    ## model 1 full data: no adjustment ====
    model1 <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                      cluster(idepic), 
                    data = data_analysis)  
    
    ## model 1 excl2year: no adjustment ====
    model1_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                cluster(idepic), 
                              data = data_analysis_excl2year)
    
    ## model 1 full data median_above: no adjustment ====
    model1_median_above <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                   cluster(idepic), 
                                 data = data_analysis_median_above)  
    
    ## model 1 excl2year median_above: no adjustment ====
    model1_median_above_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                             cluster(idepic), 
                                           data = data_analysis_median_above_excl2year)
    
    ## model 1 full data median_below: no adjustment ====
    model1_median_below <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                   cluster(idepic), 
                                 data = data_analysis_median_below)  
    
    ## model 1 excl2year median_below: no adjustment ====
    model1_median_below_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                             cluster(idepic), 
                                           data = data_analysis_median_below_excl2year)
    
    ## model 2 full data: center, 5 year age group ====
    model2 <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                      strata(center, age_group) +
                      cluster(idepic), 
                    data = data_analysis)  
    
    ## model 2 excl2year: center, 5 year age group ====
    model2_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                strata(center, age_group) +
                                cluster(idepic), 
                              data = data_analysis_excl2year)  
    
    ## model 2 full data median_above: center, 5 year age group ====
    model2_median_above <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                   strata(center, age_group) +
                                   cluster(idepic), 
                                 data = data_analysis_median_above)  
    
    ## model 2 excl2year median_above: center, 5 year age group ====
    model2_median_above_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                             strata(center, age_group) +
                                             cluster(idepic), 
                                           data = data_analysis_median_above_excl2year)  
    
    ## model 2 full data median_below: center, 5 year age group ====
    model2_median_below <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                   strata(center, age_group) +
                                   cluster(idepic), 
                                 data = data_analysis_median_below)  
    
    ## model 2 excl2year median_below: center, 5 year age group ====
    model2_median_below_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                             strata(center, age_group) +
                                             cluster(idepic), 
                                           data = data_analysis_median_below_excl2year)  
    
    ## model 3 full data: center, 5 year age group, fasting, bmi, alcohol, smoking ====
    model3 <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                      fasting_c + bmi_adj + alc_re + smoke_stat + 
                      strata(center, age_group) +
                      cluster(idepic), 
                    data = data_analysis)  
    
    ## model 3 excl2year: center, 5 year age group, fasting, bmi, alcohol, smoking ====
    model3_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                fasting_c + bmi_adj + alc_re + smoke_stat + 
                                strata(center, age_group) +
                                cluster(idepic), 
                              data = data_analysis_excl2year) 
    
    ## model 3 full data median_above: center, 5 year age group, fasting, bmi, alcohol, smoking ====
    model3_median_above <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                   fasting_c + bmi_adj + alc_re + smoke_stat + 
                                   strata(center, age_group) +
                                   cluster(idepic), 
                                 data = data_analysis_median_above)  
    
    ## model 3 excl2year median_above: center, 5 year age group, fasting, bmi, alcohol, smoking ====
    model3_median_above_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                             fasting_c + bmi_adj + alc_re + smoke_stat + 
                                             strata(center, age_group) +
                                             cluster(idepic), 
                                           data = data_analysis_median_above_excl2year) 
    
    ## model 3 full data median_below: center, 5 year age group, fasting, bmi, alcohol, smoking ====
    model3_median_below <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                   fasting_c + bmi_adj + alc_re + smoke_stat + 
                                   strata(center, age_group) +
                                   cluster(idepic), 
                                 data = data_analysis_median_below)  
    
    ## model 3 excl2year median_below: center, 5 year age group, fasting, bmi, alcohol, smoking ====
    model3_median_below_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                             fasting_c + bmi_adj + alc_re + smoke_stat + 
                                             strata(center, age_group) +
                                             cluster(idepic), 
                                           data = data_analysis_median_below_excl2year) 
    
    ## model 4 full data: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education ====
    model4 <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                      fasting_c + bmi_adj + alc_re + smoke_stat + pa_index + l_school + 
                      strata(center, age_group) +
                      cluster(idepic), 
                    data = data_analysis)  
    
    ## model 4 excl2year: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education ====
    model4_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                fasting_c + bmi_adj + alc_re + smoke_stat + pa_index + l_school + 
                                strata(center, age_group) +
                                cluster(idepic), 
                              data = data_analysis_excl2year) 
    
    ## model 4 full data median_above: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education ====
    model4_median_above <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                   fasting_c + bmi_adj + alc_re + smoke_stat + pa_index + l_school + 
                                   strata(center, age_group) +
                                   cluster(idepic), 
                                 data = data_analysis_median_above)  
    
    ## model 4 excl2year median_above: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education ====
    model4_median_above_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                             fasting_c + bmi_adj + alc_re + smoke_stat + pa_index + l_school + 
                                             strata(center, age_group) +
                                             cluster(idepic), 
                                           data = data_analysis_median_above_excl2year) 
    
    ## model 4 full data median_below: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education ====
    model4_median_below <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                   fasting_c + bmi_adj + alc_re + smoke_stat + pa_index + l_school + 
                                   strata(center, age_group) +
                                   cluster(idepic), 
                                 data = data_analysis_median_below)  
    
    ## model 4 excl2year median_below: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education ====
    model4_median_below_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                             fasting_c + bmi_adj + alc_re + smoke_stat + pa_index + l_school + 
                                             strata(center, age_group) +
                                             cluster(idepic), 
                                           data = data_analysis_median_below_excl2year) 
    
    ## model 5 full data: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre ====
    model5 <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                      fasting_c + bmi_adj + alc_re + smoke_stat + pa_index + l_school + qe_us_energy + qe_us_fat + qe_us_fibtg +
                      cluster(idepic), 
                    data = data_analysis)  
    
    ## model 5 excl2year: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre ====
    model5_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                fasting_c + 
                                cluster(idepic), 
                              data = data_analysis_excl2year) 
    
    ## model 5 full data median_above: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre, ====
    model5_median_above <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                   fasting_c + 
                                   cluster(idepic), 
                                 data = data_analysis_median_above)  
    
    ## model 5 excl2year median_above: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre ====
    model5_median_above_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                             fasting_c + 
                                             cluster(idepic), 
                                           data = data_analysis_median_above_excl2year) 
    
    ## model 5 full data median_below: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre ====
    model5_median_below <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                   fasting_c + 
                                   cluster(idepic), 
                                 data = data_analysis_median_below)  
    
    ## model 5 excl2year median_below: center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre ====
    model5_median_below_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                             fasting_c + 
                                             cluster(idepic), 
                                           data = data_analysis_median_below_excl2year) 
    
    ## make temp results tables ====
    ### model 1 ====
    table_temp_model1 <- data.frame(
      exposure = exposure,
      n = model1[["n"]],
      nevent = model1[["nevent"]],
      coef = summary(model1)$coef[1,1],
      coef_exp = summary(model1)$conf.int[1,1],
      se = summary(model1)$coef[1,3],
      ci_lower_exp = summary(model1)$conf.int[1,3],
      ci_upper_exp = summary(model1)$conf.int[1,4],
      se_robust = summary(model1)$coef[1,4],
      pval = summary(model1)$coef[1,6],
      degrees_freedom = summary(model1)$logtest[["df"]],
      concordance = summary(model1)$concordance[["C"]],
      concordance_se = summary(model1)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model1)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model1)$logtest[["pvalue"]],
      wald_test = summary(model1)$waldtest[["test"]],
      wald_test_pval = summary(model1)$waldtest[["pvalue"]],
      score_test = summary(model1)$sctest[["test"]],
      score_test_pval = summary(model1)$sctest[["pvalue"]],
      robust_score_test = summary(model1)$robscore[["test"]],
      robust_score_test_pval = summary(model1)$robscore[["pvalue"]],
      exclusion_missing = length(model1[["na.action"]])
    )
    
    table_temp_model1_excl2year <- data.frame(
      exposure = exposure,
      n = model1_excl2year[["n"]],
      nevent = model1_excl2year[["nevent"]],
      coef = summary(model1_excl2year)$coef[1,1],
      coef_exp = summary(model1_excl2year)$conf.int[1,1],
      se = summary(model1_excl2year)$coef[1,3],
      ci_lower_exp = summary(model1_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model1_excl2year)$conf.int[1,4],
      se_robust = summary(model1_excl2year)$coef[1,4],
      pval = summary(model1_excl2year)$coef[1,6],
      degrees_freedom = summary(model1_excl2year)$logtest[["df"]],
      concordance = summary(model1_excl2year)$concordance[["C"]],
      concordance_se = summary(model1_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model1_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model1_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model1_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model1_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model1_excl2year)$sctest[["test"]],
      score_test_pval = summary(model1_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model1_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model1_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model1_excl2year[["na.action"]])
    )
    
    table_temp_model1_median_above <- data.frame(
      exposure = exposure,
      n = model1_median_above[["n"]],
      nevent = model1_median_above[["nevent"]],
      coef = summary(model1_median_above)$coef[1,1],
      coef_exp = summary(model1_median_above)$conf.int[1,1],
      se = summary(model1_median_above)$coef[1,3],
      ci_lower_exp = summary(model1_median_above)$conf.int[1,3],
      ci_upper_exp = summary(model1_median_above)$conf.int[1,4],
      se_robust = summary(model1_median_above)$coef[1,4],
      pval = summary(model1_median_above)$coef[1,6],
      degrees_freedom = summary(model1_median_above)$logtest[["df"]],
      concordance = summary(model1_median_above)$concordance[["C"]],
      concordance_se = summary(model1_median_above)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model1_median_above)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model1_median_above)$logtest[["pvalue"]],
      wald_test = summary(model1_median_above)$waldtest[["test"]],
      wald_test_pval = summary(model1_median_above)$waldtest[["pvalue"]],
      score_test = summary(model1_median_above)$sctest[["test"]],
      score_test_pval = summary(model1_median_above)$sctest[["pvalue"]],
      robust_score_test = summary(model1_median_above)$robscore[["test"]],
      robust_score_test_pval = summary(model1_median_above)$robscore[["pvalue"]],
      exclusion_missing = length(model1_median_above[["na.action"]])
    )
    
    table_temp_model1_median_above_excl2year <- data.frame(
      exposure = exposure,
      n = model1_median_above_excl2year[["n"]],
      nevent = model1_median_above_excl2year[["nevent"]],
      coef = summary(model1_median_above_excl2year)$coef[1,1],
      coef_exp = summary(model1_median_above_excl2year)$conf.int[1,1],
      se = summary(model1_median_above_excl2year)$coef[1,3],
      ci_lower_exp = summary(model1_median_above_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model1_median_above_excl2year)$conf.int[1,4],
      se_robust = summary(model1_median_above_excl2year)$coef[1,4],
      pval = summary(model1_median_above_excl2year)$coef[1,6],
      degrees_freedom = summary(model1_median_above_excl2year)$logtest[["df"]],
      concordance = summary(model1_median_above_excl2year)$concordance[["C"]],
      concordance_se = summary(model1_median_above_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model1_median_above_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model1_median_above_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model1_median_above_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model1_median_above_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model1_median_above_excl2year)$sctest[["test"]],
      score_test_pval = summary(model1_median_above_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model1_median_above_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model1_median_above_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model1_median_above_excl2year[["na.action"]])
    )
    
    table_temp_model1_median_below <- data.frame(
      exposure = exposure,
      n = model1_median_below[["n"]],
      nevent = model1_median_below[["nevent"]],
      coef = summary(model1_median_below)$coef[1,1],
      coef_exp = summary(model1_median_below)$conf.int[1,1],
      se = summary(model1_median_below)$coef[1,3],
      ci_lower_exp = summary(model1_median_below)$conf.int[1,3],
      ci_upper_exp = summary(model1_median_below)$conf.int[1,4],
      se_robust = summary(model1_median_below)$coef[1,4],
      pval = summary(model1_median_below)$coef[1,6],
      degrees_freedom = summary(model1_median_below)$logtest[["df"]],
      concordance = summary(model1_median_below)$concordance[["C"]],
      concordance_se = summary(model1_median_below)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model1_median_below)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model1_median_below)$logtest[["pvalue"]],
      wald_test = summary(model1_median_below)$waldtest[["test"]],
      wald_test_pval = summary(model1_median_below)$waldtest[["pvalue"]],
      score_test = summary(model1_median_below)$sctest[["test"]],
      score_test_pval = summary(model1_median_below)$sctest[["pvalue"]],
      robust_score_test = summary(model1_median_below)$robscore[["test"]],
      robust_score_test_pval = summary(model1_median_below)$robscore[["pvalue"]],
      exclusion_missing = length(model1_median_below[["na.action"]])
    )
    
    table_temp_model1_median_below_excl2year <- data.frame(
      exposure = exposure,
      n = model1_median_below_excl2year[["n"]],
      nevent = model1_median_below_excl2year[["nevent"]],
      coef = summary(model1_median_below_excl2year)$coef[1,1],
      coef_exp = summary(model1_median_below_excl2year)$conf.int[1,1],
      se = summary(model1_median_below_excl2year)$coef[1,3],
      ci_lower_exp = summary(model1_median_below_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model1_median_below_excl2year)$conf.int[1,4],
      se_robust = summary(model1_median_below_excl2year)$coef[1,4],
      pval = summary(model1_median_below_excl2year)$coef[1,6],
      degrees_freedom = summary(model1_median_below_excl2year)$logtest[["df"]],
      concordance = summary(model1_median_below_excl2year)$concordance[["C"]],
      concordance_se = summary(model1_median_below_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model1_median_below_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model1_median_below_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model1_median_below_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model1_median_below_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model1_median_below_excl2year)$sctest[["test"]],
      score_test_pval = summary(model1_median_below_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model1_median_below_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model1_median_below_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model1_median_below_excl2year[["na.action"]])
    )
    
    ### model 2 ====
    table_temp_model2 <- data.frame(
      exposure = exposure,
      n = model2[["n"]],
      nevent = model2[["nevent"]],
      coef = summary(model2)$coef[1,1],
      coef_exp = summary(model2)$conf.int[1,1],
      se = summary(model2)$coef[1,3],
      ci_lower_exp = summary(model2)$conf.int[1,3],
      ci_upper_exp = summary(model2)$conf.int[1,4],
      se_robust = summary(model2)$coef[1,4],
      pval = summary(model2)$coef[1,6],
      degrees_freedom = summary(model2)$logtest[["df"]],
      concordance = summary(model2)$concordance[["C"]],
      concordance_se = summary(model2)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model2)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model2)$logtest[["pvalue"]],
      wald_test = summary(model2)$waldtest[["test"]],
      wald_test_pval = summary(model2)$waldtest[["pvalue"]],
      score_test = summary(model2)$sctest[["test"]],
      score_test_pval = summary(model2)$sctest[["pvalue"]],
      robust_score_test = summary(model2)$robscore[["test"]],
      robust_score_test_pval = summary(model2)$robscore[["pvalue"]],
      exclusion_missing = length(model2[["na.action"]])
    )
    
    table_temp_model2_excl2year <- data.frame(
      exposure = exposure,
      n = model2_excl2year[["n"]],
      nevent = model2_excl2year[["nevent"]],
      coef = summary(model2_excl2year)$coef[1,1],
      coef_exp = summary(model2_excl2year)$conf.int[1,1],
      se = summary(model2_excl2year)$coef[1,3],
      ci_lower_exp = summary(model2_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model2_excl2year)$conf.int[1,4],
      se_robust = summary(model2_excl2year)$coef[1,4],
      pval = summary(model2_excl2year)$coef[1,6],
      degrees_freedom = summary(model2_excl2year)$logtest[["df"]],
      concordance = summary(model2_excl2year)$concordance[["C"]],
      concordance_se = summary(model2_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model2_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model2_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model2_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model2_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model2_excl2year)$sctest[["test"]],
      score_test_pval = summary(model2_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model2_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model2_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model2_excl2year[["na.action"]])
    )
    
    table_temp_model2_median_above <- data.frame(
      exposure = exposure,
      n = model2_median_above[["n"]],
      nevent = model2_median_above[["nevent"]],
      coef = summary(model2_median_above)$coef[1,1],
      coef_exp = summary(model2_median_above)$conf.int[1,1],
      se = summary(model2_median_above)$coef[1,3],
      ci_lower_exp = summary(model2_median_above)$conf.int[1,3],
      ci_upper_exp = summary(model2_median_above)$conf.int[1,4],
      se_robust = summary(model2_median_above)$coef[1,4],
      pval = summary(model2_median_above)$coef[1,6],
      degrees_freedom = summary(model2_median_above)$logtest[["df"]],
      concordance = summary(model2_median_above)$concordance[["C"]],
      concordance_se = summary(model2_median_above)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model2_median_above)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model2_median_above)$logtest[["pvalue"]],
      wald_test = summary(model2_median_above)$waldtest[["test"]],
      wald_test_pval = summary(model2_median_above)$waldtest[["pvalue"]],
      score_test = summary(model2_median_above)$sctest[["test"]],
      score_test_pval = summary(model2_median_above)$sctest[["pvalue"]],
      robust_score_test = summary(model2_median_above)$robscore[["test"]],
      robust_score_test_pval = summary(model2_median_above)$robscore[["pvalue"]],
      exclusion_missing = length(model2_median_above[["na.action"]])
    )
    
    table_temp_model2_median_above_excl2year <- data.frame(
      exposure = exposure,
      n = model2_median_above_excl2year[["n"]],
      nevent = model2_median_above_excl2year[["nevent"]],
      coef = summary(model2_median_above_excl2year)$coef[1,1],
      coef_exp = summary(model2_median_above_excl2year)$conf.int[1,1],
      se = summary(model2_median_above_excl2year)$coef[1,3],
      ci_lower_exp = summary(model2_median_above_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model2_median_above_excl2year)$conf.int[1,4],
      se_robust = summary(model2_median_above_excl2year)$coef[1,4],
      pval = summary(model2_median_above_excl2year)$coef[1,6],
      degrees_freedom = summary(model2_median_above_excl2year)$logtest[["df"]],
      concordance = summary(model2_median_above_excl2year)$concordance[["C"]],
      concordance_se = summary(model2_median_above_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model2_median_above_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model2_median_above_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model2_median_above_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model2_median_above_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model2_median_above_excl2year)$sctest[["test"]],
      score_test_pval = summary(model2_median_above_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model2_median_above_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model2_median_above_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model2_median_above_excl2year[["na.action"]])
    )
    
    table_temp_model2_median_below <- data.frame(
      exposure = exposure,
      n = model2_median_below[["n"]],
      nevent = model2_median_below[["nevent"]],
      coef = summary(model2_median_below)$coef[1,1],
      coef_exp = summary(model2_median_below)$conf.int[1,1],
      se = summary(model2_median_below)$coef[1,3],
      ci_lower_exp = summary(model2_median_below)$conf.int[1,3],
      ci_upper_exp = summary(model2_median_below)$conf.int[1,4],
      se_robust = summary(model2_median_below)$coef[1,4],
      pval = summary(model2_median_below)$coef[1,6],
      degrees_freedom = summary(model2_median_below)$logtest[["df"]],
      concordance = summary(model2_median_below)$concordance[["C"]],
      concordance_se = summary(model2_median_below)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model2_median_below)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model2_median_below)$logtest[["pvalue"]],
      wald_test = summary(model2_median_below)$waldtest[["test"]],
      wald_test_pval = summary(model2_median_below)$waldtest[["pvalue"]],
      score_test = summary(model2_median_below)$sctest[["test"]],
      score_test_pval = summary(model2_median_below)$sctest[["pvalue"]],
      robust_score_test = summary(model2_median_below)$robscore[["test"]],
      robust_score_test_pval = summary(model2_median_below)$robscore[["pvalue"]],
      exclusion_missing = length(model2_median_below[["na.action"]])
    )
    
    table_temp_model2_median_below_excl2year <- data.frame(
      exposure = exposure,
      n = model2_median_below_excl2year[["n"]],
      nevent = model2_median_below_excl2year[["nevent"]],
      coef = summary(model2_median_below_excl2year)$coef[1,1],
      coef_exp = summary(model2_median_below_excl2year)$conf.int[1,1],
      se = summary(model2_median_below_excl2year)$coef[1,3],
      ci_lower_exp = summary(model2_median_below_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model2_median_below_excl2year)$conf.int[1,4],
      se_robust = summary(model2_median_below_excl2year)$coef[1,4],
      pval = summary(model2_median_below_excl2year)$coef[1,6],
      degrees_freedom = summary(model2_median_below_excl2year)$logtest[["df"]],
      concordance = summary(model2_median_below_excl2year)$concordance[["C"]],
      concordance_se = summary(model2_median_below_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model2_median_below_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model2_median_below_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model2_median_below_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model2_median_below_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model2_median_below_excl2year)$sctest[["test"]],
      score_test_pval = summary(model2_median_below_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model2_median_below_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model2_median_below_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model2_median_below_excl2year[["na.action"]])
    )
    ### model 3 ====
    table_temp_model3 <- data.frame(
      exposure = exposure,
      n = model3[["n"]],
      nevent = model3[["nevent"]],
      coef = summary(model3)$coef[1,1],
      coef_exp = summary(model3)$conf.int[1,1],
      se = summary(model3)$coef[1,3],
      ci_lower_exp = summary(model3)$conf.int[1,3],
      ci_upper_exp = summary(model3)$conf.int[1,4],
      se_robust = summary(model3)$coef[1,4],
      pval = summary(model3)$coef[1,6],
      degrees_freedom = summary(model3)$logtest[["df"]],
      concordance = summary(model3)$concordance[["C"]],
      concordance_se = summary(model3)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model3)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model3)$logtest[["pvalue"]],
      wald_test = summary(model3)$waldtest[["test"]],
      wald_test_pval = summary(model3)$waldtest[["pvalue"]],
      score_test = summary(model3)$sctest[["test"]],
      score_test_pval = summary(model3)$sctest[["pvalue"]],
      robust_score_test = summary(model3)$robscore[["test"]],
      robust_score_test_pval = summary(model3)$robscore[["pvalue"]],
      exclusion_missing = length(model3[["na.action"]])
    )
    
    table_temp_model3_excl2year <- data.frame(
      exposure = exposure,
      n = model3_excl2year[["n"]],
      nevent = model3_excl2year[["nevent"]],
      coef = summary(model3_excl2year)$coef[1,1],
      coef_exp = summary(model3_excl2year)$conf.int[1,1],
      se = summary(model3_excl2year)$coef[1,3],
      ci_lower_exp = summary(model3_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model3_excl2year)$conf.int[1,4],
      se_robust = summary(model3_excl2year)$coef[1,4],
      pval = summary(model3_excl2year)$coef[1,6],
      degrees_freedom = summary(model3_excl2year)$logtest[["df"]],
      concordance = summary(model3_excl2year)$concordance[["C"]],
      concordance_se = summary(model3_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model3_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model3_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model3_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model3_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model3_excl2year)$sctest[["test"]],
      score_test_pval = summary(model3_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model3_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model3_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model3_excl2year[["na.action"]])
    )
    
    table_temp_model3_median_above <- data.frame(
      exposure = exposure,
      n = model3_median_above[["n"]],
      nevent = model3_median_above[["nevent"]],
      coef = summary(model3_median_above)$coef[1,1],
      coef_exp = summary(model3_median_above)$conf.int[1,1],
      se = summary(model3_median_above)$coef[1,3],
      ci_lower_exp = summary(model3_median_above)$conf.int[1,3],
      ci_upper_exp = summary(model3_median_above)$conf.int[1,4],
      se_robust = summary(model3_median_above)$coef[1,4],
      pval = summary(model3_median_above)$coef[1,6],
      degrees_freedom = summary(model3_median_above)$logtest[["df"]],
      concordance = summary(model3_median_above)$concordance[["C"]],
      concordance_se = summary(model3_median_above)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model3_median_above)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model3_median_above)$logtest[["pvalue"]],
      wald_test = summary(model3_median_above)$waldtest[["test"]],
      wald_test_pval = summary(model3_median_above)$waldtest[["pvalue"]],
      score_test = summary(model3_median_above)$sctest[["test"]],
      score_test_pval = summary(model3_median_above)$sctest[["pvalue"]],
      robust_score_test = summary(model3_median_above)$robscore[["test"]],
      robust_score_test_pval = summary(model3_median_above)$robscore[["pvalue"]],
      exclusion_missing = length(model3_median_above[["na.action"]])
    )
    
    table_temp_model3_median_above_excl2year <- data.frame(
      exposure = exposure,
      n = model3_median_above_excl2year[["n"]],
      nevent = model3_median_above_excl2year[["nevent"]],
      coef = summary(model3_median_above_excl2year)$coef[1,1],
      coef_exp = summary(model3_median_above_excl2year)$conf.int[1,1],
      se = summary(model3_median_above_excl2year)$coef[1,3],
      ci_lower_exp = summary(model3_median_above_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model3_median_above_excl2year)$conf.int[1,4],
      se_robust = summary(model3_median_above_excl2year)$coef[1,4],
      pval = summary(model3_median_above_excl2year)$coef[1,6],
      degrees_freedom = summary(model3_median_above_excl2year)$logtest[["df"]],
      concordance = summary(model3_median_above_excl2year)$concordance[["C"]],
      concordance_se = summary(model3_median_above_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model3_median_above_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model3_median_above_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model3_median_above_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model3_median_above_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model3_median_above_excl2year)$sctest[["test"]],
      score_test_pval = summary(model3_median_above_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model3_median_above_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model3_median_above_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model3_median_above_excl2year[["na.action"]])
    )
    
    table_temp_model3_median_below <- data.frame(
      exposure = exposure,
      n = model3_median_below[["n"]],
      nevent = model3_median_below[["nevent"]],
      coef = summary(model3_median_below)$coef[1,1],
      coef_exp = summary(model3_median_below)$conf.int[1,1],
      se = summary(model3_median_below)$coef[1,3],
      ci_lower_exp = summary(model3_median_below)$conf.int[1,3],
      ci_upper_exp = summary(model3_median_below)$conf.int[1,4],
      se_robust = summary(model3_median_below)$coef[1,4],
      pval = summary(model3_median_below)$coef[1,6],
      degrees_freedom = summary(model3_median_below)$logtest[["df"]],
      concordance = summary(model3_median_below)$concordance[["C"]],
      concordance_se = summary(model3_median_below)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model3_median_below)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model3_median_below)$logtest[["pvalue"]],
      wald_test = summary(model3_median_below)$waldtest[["test"]],
      wald_test_pval = summary(model3_median_below)$waldtest[["pvalue"]],
      score_test = summary(model3_median_below)$sctest[["test"]],
      score_test_pval = summary(model3_median_below)$sctest[["pvalue"]],
      robust_score_test = summary(model3_median_below)$robscore[["test"]],
      robust_score_test_pval = summary(model3_median_below)$robscore[["pvalue"]],
      exclusion_missing = length(model3_median_below[["na.action"]])
    )
    
    table_temp_model3_median_below_excl2year <- data.frame(
      exposure = exposure,
      n = model3_median_below_excl2year[["n"]],
      nevent = model3_median_below_excl2year[["nevent"]],
      coef = summary(model3_median_below_excl2year)$coef[1,1],
      coef_exp = summary(model3_median_below_excl2year)$conf.int[1,1],
      se = summary(model3_median_below_excl2year)$coef[1,3],
      ci_lower_exp = summary(model3_median_below_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model3_median_below_excl2year)$conf.int[1,4],
      se_robust = summary(model3_median_below_excl2year)$coef[1,4],
      pval = summary(model3_median_below_excl2year)$coef[1,6],
      degrees_freedom = summary(model3_median_below_excl2year)$logtest[["df"]],
      concordance = summary(model3_median_below_excl2year)$concordance[["C"]],
      concordance_se = summary(model3_median_below_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model3_median_below_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model3_median_below_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model3_median_below_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model3_median_below_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model3_median_below_excl2year)$sctest[["test"]],
      score_test_pval = summary(model3_median_below_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model3_median_below_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model3_median_below_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model3_median_below_excl2year[["na.action"]])
    )
    
    ### model 4 ====
    table_temp_model4 <- data.frame(
      exposure = exposure,
      n = model4[["n"]],
      nevent = model4[["nevent"]],
      coef = summary(model4)$coef[1,1],
      coef_exp = summary(model4)$conf.int[1,1],
      se = summary(model4)$coef[1,3],
      ci_lower_exp = summary(model4)$conf.int[1,3],
      ci_upper_exp = summary(model4)$conf.int[1,4],
      se_robust = summary(model4)$coef[1,4],
      pval = summary(model4)$coef[1,6],
      degrees_freedom = summary(model4)$logtest[["df"]],
      concordance = summary(model4)$concordance[["C"]],
      concordance_se = summary(model4)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model4)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model4)$logtest[["pvalue"]],
      wald_test = summary(model4)$waldtest[["test"]],
      wald_test_pval = summary(model4)$waldtest[["pvalue"]],
      score_test = summary(model4)$sctest[["test"]],
      score_test_pval = summary(model4)$sctest[["pvalue"]],
      robust_score_test = summary(model4)$robscore[["test"]],
      robust_score_test_pval = summary(model4)$robscore[["pvalue"]],
      exclusion_missing = length(model4[["na.action"]])
    )
    
    table_temp_model4_excl2year <- data.frame(
      exposure = exposure,
      n = model4_excl2year[["n"]],
      nevent = model4_excl2year[["nevent"]],
      coef = summary(model4_excl2year)$coef[1,1],
      coef_exp = summary(model4_excl2year)$conf.int[1,1],
      se = summary(model4_excl2year)$coef[1,3],
      ci_lower_exp = summary(model4_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model4_excl2year)$conf.int[1,4],
      se_robust = summary(model4_excl2year)$coef[1,4],
      pval = summary(model4_excl2year)$coef[1,6],
      degrees_freedom = summary(model4_excl2year)$logtest[["df"]],
      concordance = summary(model4_excl2year)$concordance[["C"]],
      concordance_se = summary(model4_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model4_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model4_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model4_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model4_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model4_excl2year)$sctest[["test"]],
      score_test_pval = summary(model4_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model4_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model4_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model4_excl2year[["na.action"]])
    )
    
    table_temp_model4_median_above <- data.frame(
      exposure = exposure,
      n = model4_median_above[["n"]],
      nevent = model4_median_above[["nevent"]],
      coef = summary(model4_median_above)$coef[1,1],
      coef_exp = summary(model4_median_above)$conf.int[1,1],
      se = summary(model4_median_above)$coef[1,3],
      ci_lower_exp = summary(model4_median_above)$conf.int[1,3],
      ci_upper_exp = summary(model4_median_above)$conf.int[1,4],
      se_robust = summary(model4_median_above)$coef[1,4],
      pval = summary(model4_median_above)$coef[1,6],
      degrees_freedom = summary(model4_median_above)$logtest[["df"]],
      concordance = summary(model4_median_above)$concordance[["C"]],
      concordance_se = summary(model4_median_above)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model4_median_above)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model4_median_above)$logtest[["pvalue"]],
      wald_test = summary(model4_median_above)$waldtest[["test"]],
      wald_test_pval = summary(model4_median_above)$waldtest[["pvalue"]],
      score_test = summary(model4_median_above)$sctest[["test"]],
      score_test_pval = summary(model4_median_above)$sctest[["pvalue"]],
      robust_score_test = summary(model4_median_above)$robscore[["test"]],
      robust_score_test_pval = summary(model4_median_above)$robscore[["pvalue"]],
      exclusion_missing = length(model4_median_above[["na.action"]])
    )
    
    table_temp_model4_median_above_excl2year <- data.frame(
      exposure = exposure,
      n = model4_median_above_excl2year[["n"]],
      nevent = model4_median_above_excl2year[["nevent"]],
      coef = summary(model4_median_above_excl2year)$coef[1,1],
      coef_exp = summary(model4_median_above_excl2year)$conf.int[1,1],
      se = summary(model4_median_above_excl2year)$coef[1,3],
      ci_lower_exp = summary(model4_median_above_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model4_median_above_excl2year)$conf.int[1,4],
      se_robust = summary(model4_median_above_excl2year)$coef[1,4],
      pval = summary(model4_median_above_excl2year)$coef[1,6],
      degrees_freedom = summary(model4_median_above_excl2year)$logtest[["df"]],
      concordance = summary(model4_median_above_excl2year)$concordance[["C"]],
      concordance_se = summary(model4_median_above_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model4_median_above_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model4_median_above_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model4_median_above_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model4_median_above_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model4_median_above_excl2year)$sctest[["test"]],
      score_test_pval = summary(model4_median_above_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model4_median_above_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model4_median_above_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model4_median_above_excl2year[["na.action"]])
    )
    
    table_temp_model4_median_below <- data.frame(
      exposure = exposure,
      n = model4_median_below[["n"]],
      nevent = model4_median_below[["nevent"]],
      coef = summary(model4_median_below)$coef[1,1],
      coef_exp = summary(model4_median_below)$conf.int[1,1],
      se = summary(model4_median_below)$coef[1,3],
      ci_lower_exp = summary(model4_median_below)$conf.int[1,3],
      ci_upper_exp = summary(model4_median_below)$conf.int[1,4],
      se_robust = summary(model4_median_below)$coef[1,4],
      pval = summary(model4_median_below)$coef[1,6],
      degrees_freedom = summary(model4_median_below)$logtest[["df"]],
      concordance = summary(model4_median_below)$concordance[["C"]],
      concordance_se = summary(model4_median_below)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model4_median_below)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model4_median_below)$logtest[["pvalue"]],
      wald_test = summary(model4_median_below)$waldtest[["test"]],
      wald_test_pval = summary(model4_median_below)$waldtest[["pvalue"]],
      score_test = summary(model4_median_below)$sctest[["test"]],
      score_test_pval = summary(model4_median_below)$sctest[["pvalue"]],
      robust_score_test = summary(model4_median_below)$robscore[["test"]],
      robust_score_test_pval = summary(model4_median_below)$robscore[["pvalue"]],
      exclusion_missing = length(model4_median_below[["na.action"]])
    )
    
    table_temp_model4_median_below_excl2year <- data.frame(
      exposure = exposure,
      n = model4_median_below_excl2year[["n"]],
      nevent = model4_median_below_excl2year[["nevent"]],
      coef = summary(model4_median_below_excl2year)$coef[1,1],
      coef_exp = summary(model4_median_below_excl2year)$conf.int[1,1],
      se = summary(model4_median_below_excl2year)$coef[1,3],
      ci_lower_exp = summary(model4_median_below_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model4_median_below_excl2year)$conf.int[1,4],
      se_robust = summary(model4_median_below_excl2year)$coef[1,4],
      pval = summary(model4_median_below_excl2year)$coef[1,6],
      degrees_freedom = summary(model4_median_below_excl2year)$logtest[["df"]],
      concordance = summary(model4_median_below_excl2year)$concordance[["C"]],
      concordance_se = summary(model4_median_below_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model4_median_below_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model4_median_below_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model4_median_below_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model4_median_below_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model4_median_below_excl2year)$sctest[["test"]],
      score_test_pval = summary(model4_median_below_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model4_median_below_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model4_median_below_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model4_median_below_excl2year[["na.action"]])
    )
    
    ### model 5 ====
    table_temp_model5 <- data.frame(
      exposure = exposure,
      n = model5[["n"]],
      nevent = model5[["nevent"]],
      coef = summary(model5)$coef[1,1],
      coef_exp = summary(model5)$conf.int[1,1],
      se = summary(model5)$coef[1,3],
      ci_lower_exp = summary(model5)$conf.int[1,3],
      ci_upper_exp = summary(model5)$conf.int[1,4],
      se_robust = summary(model5)$coef[1,4],
      pval = summary(model5)$coef[1,6],
      degrees_freedom = summary(model5)$logtest[["df"]],
      concordance = summary(model5)$concordance[["C"]],
      concordance_se = summary(model5)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model5)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model5)$logtest[["pvalue"]],
      wald_test = summary(model5)$waldtest[["test"]],
      wald_test_pval = summary(model5)$waldtest[["pvalue"]],
      score_test = summary(model5)$sctest[["test"]],
      score_test_pval = summary(model5)$sctest[["pvalue"]],
      robust_score_test = summary(model5)$robscore[["test"]],
      robust_score_test_pval = summary(model5)$robscore[["pvalue"]],
      exclusion_missing = length(model5[["na.action"]])
    )
    
    table_temp_model5_excl2year <- data.frame(
      exposure = exposure,
      n = model5_excl2year[["n"]],
      nevent = model5_excl2year[["nevent"]],
      coef = summary(model5_excl2year)$coef[1,1],
      coef_exp = summary(model5_excl2year)$conf.int[1,1],
      se = summary(model5_excl2year)$coef[1,3],
      ci_lower_exp = summary(model5_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model5_excl2year)$conf.int[1,4],
      se_robust = summary(model5_excl2year)$coef[1,4],
      pval = summary(model5_excl2year)$coef[1,6],
      degrees_freedom = summary(model5_excl2year)$logtest[["df"]],
      concordance = summary(model5_excl2year)$concordance[["C"]],
      concordance_se = summary(model5_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model5_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model5_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model5_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model5_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model5_excl2year)$sctest[["test"]],
      score_test_pval = summary(model5_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model5_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model5_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model5_excl2year[["na.action"]])
    )
    
    table_temp_model5_median_above <- data.frame(
      exposure = exposure,
      n = model5_median_above[["n"]],
      nevent = model5_median_above[["nevent"]],
      coef = summary(model5_median_above)$coef[1,1],
      coef_exp = summary(model5_median_above)$conf.int[1,1],
      se = summary(model5_median_above)$coef[1,3],
      ci_lower_exp = summary(model5_median_above)$conf.int[1,3],
      ci_upper_exp = summary(model5_median_above)$conf.int[1,4],
      se_robust = summary(model5_median_above)$coef[1,4],
      pval = summary(model5_median_above)$coef[1,6],
      degrees_freedom = summary(model5_median_above)$logtest[["df"]],
      concordance = summary(model5_median_above)$concordance[["C"]],
      concordance_se = summary(model5_median_above)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model5_median_above)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model5_median_above)$logtest[["pvalue"]],
      wald_test = summary(model5_median_above)$waldtest[["test"]],
      wald_test_pval = summary(model5_median_above)$waldtest[["pvalue"]],
      score_test = summary(model5_median_above)$sctest[["test"]],
      score_test_pval = summary(model5_median_above)$sctest[["pvalue"]],
      robust_score_test = summary(model5_median_above)$robscore[["test"]],
      robust_score_test_pval = summary(model5_median_above)$robscore[["pvalue"]],
      exclusion_missing = length(model5_median_above[["na.action"]])
    )
    
    table_temp_model5_median_above_excl2year <- data.frame(
      exposure = exposure,
      n = model5_median_above_excl2year[["n"]],
      nevent = model5_median_above_excl2year[["nevent"]],
      coef = summary(model5_median_above_excl2year)$coef[1,1],
      coef_exp = summary(model5_median_above_excl2year)$conf.int[1,1],
      se = summary(model5_median_above_excl2year)$coef[1,3],
      ci_lower_exp = summary(model5_median_above_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model5_median_above_excl2year)$conf.int[1,4],
      se_robust = summary(model5_median_above_excl2year)$coef[1,4],
      pval = summary(model5_median_above_excl2year)$coef[1,6],
      degrees_freedom = summary(model5_median_above_excl2year)$logtest[["df"]],
      concordance = summary(model5_median_above_excl2year)$concordance[["C"]],
      concordance_se = summary(model5_median_above_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model5_median_above_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model5_median_above_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model5_median_above_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model5_median_above_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model5_median_above_excl2year)$sctest[["test"]],
      score_test_pval = summary(model5_median_above_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model5_median_above_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model5_median_above_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model5_median_above_excl2year[["na.action"]])
    )
    
    table_temp_model5_median_below <- data.frame(
      exposure = exposure,
      n = model5_median_below[["n"]],
      nevent = model5_median_below[["nevent"]],
      coef = summary(model5_median_below)$coef[1,1],
      coef_exp = summary(model5_median_below)$conf.int[1,1],
      se = summary(model5_median_below)$coef[1,3],
      ci_lower_exp = summary(model5_median_below)$conf.int[1,3],
      ci_upper_exp = summary(model5_median_below)$conf.int[1,4],
      se_robust = summary(model5_median_below)$coef[1,4],
      pval = summary(model5_median_below)$coef[1,6],
      degrees_freedom = summary(model5_median_below)$logtest[["df"]],
      concordance = summary(model5_median_below)$concordance[["C"]],
      concordance_se = summary(model5_median_below)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model5_median_below)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model5_median_below)$logtest[["pvalue"]],
      wald_test = summary(model5_median_below)$waldtest[["test"]],
      wald_test_pval = summary(model5_median_below)$waldtest[["pvalue"]],
      score_test = summary(model5_median_below)$sctest[["test"]],
      score_test_pval = summary(model5_median_below)$sctest[["pvalue"]],
      robust_score_test = summary(model5_median_below)$robscore[["test"]],
      robust_score_test_pval = summary(model5_median_below)$robscore[["pvalue"]],
      exclusion_missing = length(model5_median_below[["na.action"]])
    )
    
    table_temp_model5_median_below_excl2year <- data.frame(
      exposure = exposure,
      n = model5_median_below_excl2year[["n"]],
      nevent = model5_median_below_excl2year[["nevent"]],
      coef = summary(model5_median_below_excl2year)$coef[1,1],
      coef_exp = summary(model5_median_below_excl2year)$conf.int[1,1],
      se = summary(model5_median_below_excl2year)$coef[1,3],
      ci_lower_exp = summary(model5_median_below_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(model5_median_below_excl2year)$conf.int[1,4],
      se_robust = summary(model5_median_below_excl2year)$coef[1,4],
      pval = summary(model5_median_below_excl2year)$coef[1,6],
      degrees_freedom = summary(model5_median_below_excl2year)$logtest[["df"]],
      concordance = summary(model5_median_below_excl2year)$concordance[["C"]],
      concordance_se = summary(model5_median_below_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(model5_median_below_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(model5_median_below_excl2year)$logtest[["pvalue"]],
      wald_test = summary(model5_median_below_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(model5_median_below_excl2year)$waldtest[["pvalue"]],
      score_test = summary(model5_median_below_excl2year)$sctest[["test"]],
      score_test_pval = summary(model5_median_below_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(model5_median_below_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(model5_median_below_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(model5_median_below_excl2year[["na.action"]])
    )
    
    ## make final tables ====
    table_model1 <- bind_rows(table_model1, table_temp_model1)
    table_model2 <- bind_rows(table_model2, table_temp_model2)
    table_model3 <- bind_rows(table_model3, table_temp_model3)
    table_model4 <- bind_rows(table_model4, table_temp_model4)
    table_model5 <- bind_rows(table_model5, table_temp_model5)
    
    table_model1_excl2year <- bind_rows(table_model1_excl2year, table_temp_model1_excl2year)
    table_model2_excl2year <- bind_rows(table_model2_excl2year, table_temp_model2_excl2year)
    table_model3_excl2year <- bind_rows(table_model3_excl2year, table_temp_model3_excl2year)
    table_model4_excl2year <- bind_rows(table_model4_excl2year, table_temp_model4_excl2year)
    table_model5_excl2year <- bind_rows(table_model5_excl2year, table_temp_model5_excl2year)
    
    ## multiple testing
    table_model1 <- table_model1 %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model1 <- table_model1 %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model2 <- table_model2 %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model2 <- table_model2 %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model3 <- table_model3 %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model3 <- table_model3 %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model4 <- table_model4 %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model4 <- table_model4 %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model5 <- table_model5 %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model5 <- table_model5 %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    
    table_model1_excl2year <- table_model1_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model1_excl2year <- table_model1_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model2_excl2year <- table_model2_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model2_excl2year <- table_model2_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model3_excl2year <- table_model3_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model3_excl2year <- table_model3_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model4_excl2year <- table_model4_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model4_excl2year <- table_model4_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model5_excl2year <- table_model5_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model5_excl2year <- table_model5_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    
    ## add analysis column
    table_model1$analysis <- "complete"
    table_model2$analysis <- "complete"
    table_model3$analysis <- "complete"
    table_model4$analysis <- "complete"
    table_model5$analysis <- "complete"
    
    table_model1$followup <- "complete"
    table_model2$followup <- "complete"
    table_model3$followup <- "complete"
    table_model4$followup <- "complete"
    table_model5$followup <- "complete"
    
    table_model1$model <- "1"
    table_model2$model <- "2"
    table_model3$model <- "3"
    table_model4$model <- "4"
    table_model5$model <- "5"
    
    table_model1$covariates <- ""
    table_model2$covariates <- "center, 5 year age group"
    table_model3$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking"
    table_model4$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education"
    table_model5$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre"
    
    table_model1_excl2year$analysis <- "excl2year"
    table_model2_excl2year$analysis <- "excl2year"
    table_model3_excl2year$analysis <- "excl2year"
    table_model4_excl2year$analysis <- "excl2year"
    table_model5_excl2year$analysis <- "excl2year"
    
    table_model1_excl2year$followup <- "complete"
    table_model2_excl2year$followup <- "complete"
    table_model3_excl2year$followup <- "complete"
    table_model4_excl2year$followup <- "complete"
    table_model5_excl2year$followup <- "complete"
    
    table_model1_excl2year$model <- "1"
    table_model2_excl2year$model <- "2"
    table_model3_excl2year$model <- "3"
    table_model4_excl2year$model <- "4"
    table_model5_excl2year$model <- "5"
    
    table_model1_excl2year$covariates <- ""
    table_model2_excl2year$covariates <- "center, 5 year age group"
    table_model3_excl2year$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking"
    table_model4_excl2year$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education"
    table_model5_excl2year$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre"
    
    ## median above ====
    table_model1_median_above <- bind_rows(table_model1_median_above, table_temp_model1_median_above)
    table_model2_median_above <- bind_rows(table_model2_median_above, table_temp_model2_median_above)
    table_model3_median_above <- bind_rows(table_model3_median_above, table_temp_model3_median_above)
    table_model4_median_above <- bind_rows(table_model4_median_above, table_temp_model4_median_above)
    table_model5_median_above <- bind_rows(table_model5_median_above, table_temp_model5_median_above)
    
    table_model1_median_above_excl2year <- bind_rows(table_model1_median_above_excl2year, table_temp_model1_median_above_excl2year)
    table_model2_median_above_excl2year <- bind_rows(table_model2_median_above_excl2year, table_temp_model2_median_above_excl2year)
    table_model3_median_above_excl2year <- bind_rows(table_model3_median_above_excl2year, table_temp_model3_median_above_excl2year)
    table_model4_median_above_excl2year <- bind_rows(table_model4_median_above_excl2year, table_temp_model4_median_above_excl2year)
    table_model5_median_above_excl2year <- bind_rows(table_model5_median_above_excl2year, table_temp_model5_median_above_excl2year)
    
    ## multiple testing
    table_model1_median_above <- table_model1_median_above %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model1_median_above <- table_model1_median_above %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model2_median_above <- table_model2_median_above %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model2_median_above <- table_model2_median_above %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model3_median_above <- table_model3_median_above %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model3_median_above <- table_model3_median_above %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model4_median_above <- table_model4_median_above %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model4_median_above <- table_model4_median_above %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model5_median_above <- table_model5_median_above %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model5_median_above <- table_model5_median_above %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    
    table_model1_median_above_excl2year <- table_model1_median_above_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model1_median_above_excl2year <- table_model1_median_above_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model2_median_above_excl2year <- table_model2_median_above_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model2_median_above_excl2year <- table_model2_median_above_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model3_median_above_excl2year <- table_model3_median_above_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model3_median_above_excl2year <- table_model3_median_above_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model4_median_above_excl2year <- table_model4_median_above_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model4_median_above_excl2year <- table_model4_median_above_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model5_median_above_excl2year <- table_model5_median_above_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model5_median_above_excl2year <- table_model5_median_above_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    
    ## add analysis column
    table_model1_median_above$analysis <- "complete"
    table_model2_median_above$analysis <- "complete"
    table_model3_median_above$analysis <- "complete"
    table_model4_median_above$analysis <- "complete"
    table_model5_median_above$analysis <- "complete"
    
    table_model1_median_above$followup <- "above"
    table_model2_median_above$followup <- "above"
    table_model3_median_above$followup <- "above"
    table_model4_median_above$followup <- "above"
    table_model5_median_above$followup <- "above"
    
    table_model1_median_above$model <- "1"
    table_model2_median_above$model <- "2"
    table_model3_median_above$model <- "3"
    table_model4_median_above$model <- "4"
    table_model5_median_above$model <- "5"
    
    table_model1_median_above$covariates <- ""
    table_model2_median_above$covariates <- "center, 5 year age group"
    table_model3_median_above$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking"
    table_model4_median_above$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education"
    table_model5_median_above$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre"
    
    table_model1_median_above_excl2year$analysis <- "excl2year"
    table_model2_median_above_excl2year$analysis <- "excl2year"
    table_model3_median_above_excl2year$analysis <- "excl2year"
    table_model4_median_above_excl2year$analysis <- "excl2year"
    table_model5_median_above_excl2year$analysis <- "excl2year"
    
    table_model1_median_above_excl2year$followup <- "above"
    table_model2_median_above_excl2year$followup <- "above"
    table_model3_median_above_excl2year$followup <- "above"
    table_model4_median_above_excl2year$followup <- "above"
    table_model5_median_above_excl2year$followup <- "above"
    
    table_model1_median_above_excl2year$model <- "1"
    table_model2_median_above_excl2year$model <- "2"
    table_model3_median_above_excl2year$model <- "3"
    table_model4_median_above_excl2year$model <- "4"
    table_model5_median_above_excl2year$model <- "5"
    
    table_model1_median_above_excl2year$covariates <- ""
    table_model2_median_above_excl2year$covariates <- "center, 5 year age group"
    table_model3_median_above_excl2year$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking"
    table_model4_median_above_excl2year$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education"
    table_model5_median_above_excl2year$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre"
    
    
    ## median below ====
    table_model1_median_below <- bind_rows(table_model1_median_below, table_temp_model1_median_below)
    table_model2_median_below <- bind_rows(table_model2_median_below, table_temp_model2_median_below)
    table_model3_median_below <- bind_rows(table_model3_median_below, table_temp_model3_median_below)
    table_model4_median_below <- bind_rows(table_model4_median_below, table_temp_model4_median_below)
    table_model5_median_below <- bind_rows(table_model5_median_below, table_temp_model5_median_below)
    
    table_model1_median_below_excl2year <- bind_rows(table_model1_median_below_excl2year, table_temp_model1_median_below_excl2year)
    table_model2_median_below_excl2year <- bind_rows(table_model2_median_below_excl2year, table_temp_model2_median_below_excl2year)
    table_model3_median_below_excl2year <- bind_rows(table_model3_median_below_excl2year, table_temp_model3_median_below_excl2year)
    table_model4_median_below_excl2year <- bind_rows(table_model4_median_below_excl2year, table_temp_model4_median_below_excl2year)
    table_model5_median_below_excl2year <- bind_rows(table_model5_median_below_excl2year, table_temp_model5_median_below_excl2year)
    
    ## multiple testing
    table_model1_median_below <- table_model1_median_below %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model1_median_below <- table_model1_median_below %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model2_median_below <- table_model2_median_below %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model2_median_below <- table_model2_median_below %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model3_median_below <- table_model3_median_below %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model3_median_below <- table_model3_median_below %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model4_median_below <- table_model4_median_below %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model4_median_below <- table_model4_median_below %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model5_median_below <- table_model5_median_below %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model5_median_below <- table_model5_median_below %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    
    table_model1_median_below_excl2year <- table_model1_median_below_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model1_median_below_excl2year <- table_model1_median_below_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model2_median_below_excl2year <- table_model2_median_below_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model2_median_below_excl2year <- table_model2_median_below_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model3_median_below_excl2year <- table_model3_median_below_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model3_median_below_excl2year <- table_model3_median_below_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model4_median_below_excl2year <- table_model4_median_below_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model4_median_below_excl2year <- table_model4_median_below_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    table_model5_median_below_excl2year <- table_model5_median_below_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_model5_median_below_excl2year <- table_model5_median_below_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    
    ## add analysis column
    table_model1_median_below$analysis <- "complete"
    table_model2_median_below$analysis <- "complete"
    table_model3_median_below$analysis <- "complete"
    table_model4_median_below$analysis <- "complete"
    table_model5_median_below$analysis <- "complete"
    
    table_model1_median_below$followup <- "below"
    table_model2_median_below$followup <- "below"
    table_model3_median_below$followup <- "below"
    table_model4_median_below$followup <- "below"
    table_model5_median_below$followup <- "below"
    
    table_model1_median_below$model <- "1"
    table_model2_median_below$model <- "2"
    table_model3_median_below$model <- "3"
    table_model4_median_below$model <- "4"
    table_model5_median_below$model <- "5"
    
    table_model1_median_below$covariates <- ""
    table_model2_median_below$covariates <- "center, 5 year age group"
    table_model3_median_below$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking"
    table_model4_median_below$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education"
    table_model5_median_below$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre"
    
    table_model1_median_below_excl2year$analysis <- "excl2year"
    table_model2_median_below_excl2year$analysis <- "excl2year"
    table_model3_median_below_excl2year$analysis <- "excl2year"
    table_model4_median_below_excl2year$analysis <- "excl2year"
    table_model5_median_below_excl2year$analysis <- "excl2year"
    
    table_model1_median_below_excl2year$followup <- "below"
    table_model2_median_below_excl2year$followup <- "below"
    table_model3_median_below_excl2year$followup <- "below"
    table_model4_median_below_excl2year$followup <- "below"
    table_model5_median_below_excl2year$followup <- "below"
    
    table_model1_median_below_excl2year$model <- "1"
    table_model2_median_below_excl2year$model <- "2"
    table_model3_median_below_excl2year$model <- "3"
    table_model4_median_below_excl2year$model <- "4"
    table_model5_median_below_excl2year$model <- "5"
    
    table_model1_median_below_excl2year$covariates <- ""
    table_model2_median_below_excl2year$covariates <- "center, 5 year age group"
    table_model3_median_below_excl2year$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking"
    table_model4_median_below_excl2year$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education"
    table_model5_median_below_excl2year$covariates <- "center, 5 year age group, fasting, bmi, alcohol, smoking, physical activity, education, total energy, fat, fibre"
  }
  
  ## save ====
  table <- bind_rows(
    table_model1, table_model2, table_model3, table_model4, table_model5,
    table_model1_excl2year, table_model2_excl2year, table_model3_excl2year, table_model4_excl2year, table_model5_excl2year,
    table_model1_median_above, table_model2_median_above, table_model3_median_above, table_model4_median_above, table_model5_median_above,
    table_model1_median_above_excl2year, table_model2_median_above_excl2year, table_model3_median_above_excl2year, table_model4_median_above_excl2year, table_model5_median_above_excl2year,
    table_model1_median_below, table_model2_median_below, table_model3_median_below, table_model4_median_below, table_model5_median_below,
    table_model1_median_below_excl2year, table_model2_median_below_excl2year, table_model3_median_below_excl2year, table_model4_median_below_excl2year, table_model5_median_below_excl2year
  )
  table$sex <- "female"
  write.table(table, paste0(path_out, "results.txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}
