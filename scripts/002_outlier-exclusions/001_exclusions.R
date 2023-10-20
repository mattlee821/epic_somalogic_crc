rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(data.table)
# remotes::install_github("privefl/bigutilsr")
library(bigutilsr)
library(ggplot2)
library(magrittr)
library(cowplot)
library(functions)
library(wesanderson)
library(GGally)
palette_discrete <- palette()

# lowerfun <- function(data,mapping){
#   ggplot(data = data, mapping = mapping)+
#     geom_point(size = 1, alpha = 0.8)+
#     scale_x_continuous(limits = c(-20,20))+
#     scale_y_continuous(limits = c(-20,20))
# }

# data ====
file_list <- list.files("analysis/001_phenofile/complete/", pattern = "proteins", full.names = TRUE, recursive = TRUE)
table_samples_exclude <- data.frame()
table_features_exclude <- data.frame()

# Loop through each file
for (files in file_list) {
  cat("Processing file:", files, "\n")
  
  path_in <- gsub("proteins.txt", "", files)
  path_out <- gsub("001_phenofile", "002_outlier-exclusions", path_in)
  path_label <- unlist(strsplit(path_in, "//"))
  path_label <- gsub("/", "-", path_label[length(path_label)])
  path_label <- sub("-*$", "", path_label)
  
  # Read data ====
  data <- read.table(files, header = TRUE)
  rownames(data) <- data[, 1]
  data <- data[, -1]  
  data_samples <- data
  data_features <- as.data.frame(t(data_samples))

# PCA ====
pca_samples <- prcomp(x = data_samples, rank. = 10)
data_samples_analysis <- as.data.frame(pca_samples$x)

pca_features <- prcomp(x = data_features, rank. = 10)
data_features_analysis <- as.data.frame(pca_features$x)

# samples ====
llof <- LOF(data_samples_analysis) # compute distances using local outlier factor
outliers_samples <- which(llof > tukey_mc_up(llof)) # identify outlier threshold using tukeys rule
length(outliers_samples)
id_samples_exclude <- rownames(data_samples_analysis[outliers_samples, ])
# plot PC1/2 with lof highlighting
data_samples_analysis$llof <- llof
plot <- ggpairs(data_samples_analysis,
        upper = list(continuous = "blank"),
        diag = list(continuous = "blankDiag"),
        lower = list(continuous = "points"),
        columns = 1:10,
        mapping = aes(color = llof),
        legend = 11) +
  scale_color_viridis_c(direction = -1, option = "F") +
  labs(color = paste("Local outlier factor:", round(tukey_mc_up(llof), 2))) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "grey")) +
  theme(legend.position = "bottom") +
  theme(panel.spacing = unit(1, "lines"))
tiff(paste0(path_out, "samples.tiff"), width = 800, height = 700, units = "px")
print(plot)
dev.off()

# features ====
llof <- LOF(data_features_analysis) # compute distances using local outlier factor
outliers_features <- which(llof > tukey_mc_up(llof)) # identify outlier threshold using tukeys rule
length(outliers_features)
id_features_exclude <- rownames(data_features_analysis[outliers_features, ])
# plot PC1/2 with lof highlighting
data_features_analysis$llof <- llof
plot <- ggpairs(data_features_analysis,
        upper = list(continuous = "blank"),
        diag = list(continuous = "blankDiag"),
        lower = list(continuous = "points"),
        columns = 1:10,
        mapping = aes(color = llof),
        legend = 11) +
  scale_color_viridis_c(direction = -1, option = "F") +
  labs(color = paste("Local outlier factor:", round(tukey_mc_up(llof), 2))) +
  theme_minimal() +
  theme(strip.background = element_rect(fill = "grey")) +
  theme(legend.position = "bottom") +
  theme(panel.spacing = unit(1, "lines"))
tiff(paste0(path_out, "features.tiff"), width = 800, height = 700, units = "px")
print(plot)
dev.off()

# exclude ====
if(length(outliers_samples) > 0) {data <- data[-outliers_samples, ]}
if(length(outliers_features) > 0) {data <- data[, -outliers_features]}
id_features <- names(data)
data$idepic <- rownames(data)
data <- data[, c("idepic", setdiff(names(data), "idepic"))]
id_samples <- data$idepic

# save ====
data <- fread(paste0(path_in, "phenofile.txt"), header = T)
## exclude samples
data <- data[!(data$idepic %in% id_samples_exclude), ]
## exclude features
data <- select(data, -any_of(id_features_exclude))
## save 
write.table(data, paste0(path_out, "phenofile.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# make table of exclusions ====
id_samples_exclude <- as.data.frame(id_samples_exclude)
id_features_exclude <- as.data.frame(id_features_exclude)

# Check if there are observations in either data frame
if (nrow(id_samples_exclude) > 0 || nrow(id_features_exclude) > 0) {
  # Combine table if either has observations
  if (nrow(id_samples_exclude) > 0) {
    id_samples_exclude$ID <- path_label
    table_samples_exclude <- bind_rows(table_samples_exclude, id_samples_exclude)
  }
  if (nrow(id_features_exclude) > 0) {
    id_features_exclude$ID <- path_label
    table_features_exclude <- bind_rows(table_features_exclude, id_features_exclude)
  }
}
}

# save exclusion table
write.table(table_samples_exclude, "analysis/002_outlier-exclusions/exclusions_samples.txt",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(table_features_exclude, "analysis/002_outlier-exclusions/exclusions_features.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


