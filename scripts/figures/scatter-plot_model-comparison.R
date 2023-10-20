rm(list=ls())
set.seed(821)

# enviroment ====
library(dplyr)
library(data.table)
library(tidyr)
library(functions)
library(wesanderson)
library(GGally)
palette_discrete <- palette()

# data ====
data <- fread("analysis/003_analysis/complete/results_formatted.txt")
## rearrange 
data <- select(data,
               "Target", "SeqId", "cancer",
               "sex", "analysis", "followup", "model", "covariates",
               "n", "nevent",
               "coef_exp", "se", "ci_lower_exp", "ci_upper_exp", "pval", "FDR_BH", "FDR_bonferroni")

# plot_data ====
plot_data <- data %>%
  mutate(ID = paste(Target, SeqId, cancer, sex, analysis, followup, sep = ";")) %>%
  select(ID, model, coef_exp) %>%
  pivot_wider(names_from = model, values_from = coef_exp)
plot_data <- plot_data %>%
  separate(ID, into = c("Target", "SeqId", "cancer", "sex", "analysis", "followup"), sep = ";")
plot_data$cancer <- factor(plot_data$cancer, levels = c("overall", "colon", "rectum", "early_onset"))

## split by sex ====
plot_data_combined <- subset(plot_data, sex == "combined")
plot_data_female <- subset(plot_data, sex == "female")
plot_data_male <- subset(plot_data, sex == "male")

# plot ====
plot_combined <- ggpairs(plot_data_combined, 
                upper = list(continuous = "cor"),
                diag = list(continuous = "blankDiag"),
                lower = list(continuous = "points"),
                columns = 7:11, 
                mapping = aes(color = cancer),
                legend = 6,
                title = "Model comparison: sex-combined") +
  scale_color_manual(values = palette_discrete) +
  theme(strip.background = element_rect(fill = "grey")) +
  theme(legend.position = "bottom") 

plot_female <- ggpairs(plot_data_female, 
                         upper = list(continuous = "cor"),
                         diag = list(continuous = "blankDiag"),
                         lower = list(continuous = "points"),
                         columns = 7:11, 
                         mapping = aes(color = cancer),
                         legend = 6,
                         title = "Model comparison: female") +
  scale_color_manual(values = palette_discrete) +
  theme(strip.background = element_rect(fill = "grey")) +
  theme(legend.position = "bottom") 

plot_male <- ggpairs(plot_data_male, 
                       upper = list(continuous = "cor"),
                       diag = list(continuous = "blankDiag"),
                       lower = list(continuous = "points"),
                       columns = 7:11, 
                       mapping = aes(color = cancer),
                       legend = 6,
                       title = "Model comparison: male") +
  scale_color_manual(values = palette_discrete) +
  theme(strip.background = element_rect(fill = "grey")) +
  theme(legend.position = "bottom") 

## save ====
tiff("analysis/tables/scatter-plot_model-comparison-combined.tiff", width = 800, height = 800, units = "px")
print(plot_combined)
dev.off()

tiff("analysis/tables/scatter-plot_model-comparison-female.tiff", width = 800, height = 800, units = "px")
print(plot_female)
dev.off()

tiff("analysis/tables/scatter-plot_model-comparison-male.tiff", width = 800, height = 800, units = "px")
print(plot_male)
dev.off()




