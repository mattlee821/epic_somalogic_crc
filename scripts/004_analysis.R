rm(list=ls())
set.seed(821)

# enviroment ====
library(dplyr)
library(data.table)

# data ====
data <- fread("analysis/003_data-for-analysis/complete/data.txt")
data_analysis <- data

data_analysis <- data_analysis %>%  mutate(cvd_t2d_coh = ifelse(cvd_t2d_coh == "No", 0, 1))
data_analysis <- data_analysis %>% mutate(age_group = as.integer(factor(cut(age, breaks = seq(0, max(age) + 5, by = 5), right = FALSE, labels = seq(1, max(age), by = 5)))))
data_analysis_excl2year <- data_analysis %>% filter(ageevent > age + 2) %>% mutate(age = age + 2)
table(data_analysis$indevent)
table(data_analysis_excl2year$indevent)
proteins <- names(data_analysis)[grepl("seq", names(data_analysis))]

# analysis ====
## make empty data frames for results
table_results <- data.frame()
table_results_noadj <- data.frame()
table_results_excl2year <- data.frame()
table_results_noadj_excl2year <- data.frame()

## percentage complete
  ksq <- 0
  for (exposure in proteins)
  {
    ksq <- ksq +1
    if(ksq%%100 ==1){cat("+++++++++", ksq, "\n")}
    
    ## model 1: adjustment
    results <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                   alc_re + pa_index + smoke_stat + bmi_adj + l_school + height_adj + menopause +
                   strata(sex, center, age_group) +
                   cluster(idepic), #, age_cat
                 data = data_analysis)
    
    ## model 2: no adjustment
    results_noadj <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                         cluster(idepic), #, age_cat
                       data = data_analysis)
    
    ## model3: model 1 adjustment and excluding 2 year followup
    results_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                                 alc_re + pa_index + smoke_stat + bmi_adj + l_school + height_adj + menopause +
                                 strata(sex, center, age_group) +
                                 cluster(idepic), #, age_cat
                               data = data_analysis_excl2year)
    
    ## model4: no adjustment and excluding 2 year followup
    results_noadj_excl2year <- coxph(Surv(age, ageevent, indevent) ~  get(exposure) +
                             cluster(idepic), #, age_cat
                           data = data_analysis_excl2year)
    
    ## make temp results tables    
    table_temp_results <- data.frame(
      exposure = exposure,
      n = results[["n"]],
      nevent = results[["nevent"]],
      coef = summary(results)$coef[1,1],
      coef_exp = summary(results)$conf.int[1,1],
      se = summary(results)$coef[1,3],
      ci_lower_exp = summary(results)$conf.int[1,3],
      ci_upper_exp = summary(results)$conf.int[1,4],
      se_robust = summary(results)$coef[1,4],
      pval = summary(results)$coef[1,6],
      degrees_freedom = summary(results)$logtest[["df"]],
      concordance = summary(results)$concordance[["C"]],
      concordance_se = summary(results)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(results)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(results)$logtest[["pvalue"]],
      wald_test = summary(results)$waldtest[["test"]],
      wald_test_pval = summary(results)$waldtest[["pvalue"]],
      score_test = summary(results)$sctest[["test"]],
      score_test_pval = summary(results)$sctest[["pvalue"]],
      robust_score_test = summary(results)$robscore[["test"]],
      robust_score_test_pval = summary(results)$robscore[["pvalue"]],
      exclusion_missing = length(results[["na.action"]])
    )
    
    table_temp_results_noadj <- data.frame(
      exposure = exposure,
      n = results_noadj[["n"]],
      nevent = results_noadj[["nevent"]],
      coef = summary(results_noadj)$coef[1,1],
      coef_exp = summary(results_noadj)$conf.int[1,1],
      se = summary(results_noadj)$coef[1,3],
      ci_lower_exp = summary(results_noadj)$conf.int[1,3],
      ci_upper_exp = summary(results_noadj)$conf.int[1,4],
      se_robust = summary(results_noadj)$coef[1,4],
      pval = summary(results_noadj)$coef[1,6],
      degrees_freedom = summary(results_noadj)$logtest[["df"]],
      concordance = summary(results_noadj)$concordance[["C"]],
      concordance_se = summary(results_noadj)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(results_noadj)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(results_noadj)$logtest[["pvalue"]],
      wald_test = summary(results_noadj)$waldtest[["test"]],
      wald_test_pval = summary(results_noadj)$waldtest[["pvalue"]],
      score_test = summary(results_noadj)$sctest[["test"]],
      score_test_pval = summary(results_noadj)$sctest[["pvalue"]],
      robust_score_test = summary(results_noadj)$robscore[["test"]],
      robust_score_test_pval = summary(results_noadj)$robscore[["pvalue"]],
      exclusion_missing = length(results_noadj[["na.action"]])
    )
    
    table_temp_results_excl2year <- data.frame(
      exposure = exposure,
      n = results_excl2year[["n"]],
      nevent = results_excl2year[["nevent"]],
      coef = summary(results_excl2year)$coef[1,1],
      coef_exp = summary(results_excl2year)$conf.int[1,1],
      se = summary(results_excl2year)$coef[1,3],
      ci_lower_exp = summary(results_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(results_excl2year)$conf.int[1,4],
      se_robust = summary(results_excl2year)$coef[1,4],
      pval = summary(results_excl2year)$coef[1,6],
      degrees_freedom = summary(results_excl2year)$logtest[["df"]],
      concordance = summary(results_excl2year)$concordance[["C"]],
      concordance_se = summary(results_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(results_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(results_excl2year)$logtest[["pvalue"]],
      wald_test = summary(results_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(results_excl2year)$waldtest[["pvalue"]],
      score_test = summary(results_excl2year)$sctest[["test"]],
      score_test_pval = summary(results_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(results_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(results_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(results_excl2year[["na.action"]])
    )
    
    table_temp_results_noadj_excl2year <- data.frame(
      exposure = exposure,
      n = results_noadj_excl2year[["n"]],
      nevent = results_noadj_excl2year[["nevent"]],
      coef = summary(results_noadj_excl2year)$coef[1,1],
      coef_exp = summary(results_noadj_excl2year)$conf.int[1,1],
      se = summary(results_noadj_excl2year)$coef[1,3],
      ci_lower_exp = summary(results_noadj_excl2year)$conf.int[1,3],
      ci_upper_exp = summary(results_noadj_excl2year)$conf.int[1,4],
      se_robust = summary(results_noadj_excl2year)$coef[1,4],
      pval = summary(results_noadj_excl2year)$coef[1,6],
      degrees_freedom = summary(results_noadj_excl2year)$logtest[["df"]],
      concordance = summary(results_noadj_excl2year)$concordance[["C"]],
      concordance_se = summary(results_noadj_excl2year)$concordance[["se(C)"]],
      likelihood_ratio_test = summary(results_noadj_excl2year)$logtest[["test"]],
      likelihood_ratio_test_pval = summary(results_noadj_excl2year)$logtest[["pvalue"]],
      wald_test = summary(results_noadj_excl2year)$waldtest[["test"]],
      wald_test_pval = summary(results_noadj_excl2year)$waldtest[["pvalue"]],
      score_test = summary(results_noadj_excl2year)$sctest[["test"]],
      score_test_pval = summary(results_noadj_excl2year)$sctest[["pvalue"]],
      robust_score_test = summary(results_noadj_excl2year)$robscore[["test"]],
      robust_score_test_pval = summary(results_noadj_excl2year)$robscore[["pvalue"]],
      exclusion_missing = length(results_noadj_excl2year[["na.action"]])
    )
    
    ## make final tables
    table_results <- bind_rows(table_results, table_temp_results)
    table_results_noadj <- bind_rows(table_results_noadj, table_temp_results_noadj)
    table_results_excl2year <- bind_rows(table_results_excl2year, table_temp_results_excl2year)
    table_results_noadj_excl2year <- bind_rows(table_results_noadj_excl2year, table_temp_results_noadj_excl2year)
    
    ## multiple testing
    table_results <- table_results %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_results <- table_results %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
  
    table_results_noadj <- table_results_noadj %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_results_noadj <- table_results_noadj %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    
    table_results_excl2year <- table_results_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_results_excl2year <- table_results_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    
    table_results_noadj_excl2year <- table_results_noadj_excl2year %>% mutate(FDR_BH = p.adjust(pval, method = "BH"))
    table_results_noadj_excl2year <- table_results_noadj_excl2year %>% mutate(FDR_bonferroni = p.adjust(pval, method = "bonferroni"))
    
    ## add analysis column
    table_results$analysis <- "adjustment"
    table_results_noadj$analysis <- "no-adjustment"
    table_results_excl2year$analysis <- "adjustment and excl2year"
    table_results_noadj_excl2year$analysis <- "no-adjustment and excl2year"
  }

## save ====
write.table(table_results, "analysis/004_analysis/results.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(table_results_noadj, "analysis/004_analysis/results_noadj.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(table_results_excl2year, "analysis/004_analysis/results_excl2year.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(table_results_noadj_excl2year, "analysis/004_analysis/results_noadj_excl2year.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

  