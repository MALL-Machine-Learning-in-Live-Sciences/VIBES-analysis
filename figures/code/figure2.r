rm(list = ls())
library(ggplot2)
library(corrplot)
library(phyloseq)
experiment <-  "Ravel"

# Figure 2a
# ===============
# Loada Data
data <- read.delim2(file = paste0("01_get_valencias/data/", experiment, ".csv"),
                    sep = ",")
# 1.Prepare the Data
data <- data %>%
  select(-one_of("sampleID", "read_count"))
# 2.Make corr matrix
corr_matrix <- cor(data)
svg(filename = "figures/plots/Figure2/fig2a.svg",
    width = 3.543, 
    height = 3.543,
    pointsize = 8)
corrplot(corr_matrix, order = "hclust", method = "color",
         type = 'upper', tl.pos='n')
dev.off()

# Figure 2b
# ===============
# Loada Data
data <- readRDS(file = paste0("00_preprocess_cohorts/data/SpeciesIntersect/",
                              experiment, "_Species_pseq_22.rds"))
# 1.Prepare the Data
data <- as.data.frame(data@otu_table)
# 2.Make corr matrix
corr_matrix <- cor(data)
svg(filename = "figures/plots/Figure2/fig2b.svg",
    width = 3.543, 
    height = 3.543,
    pointsize = 8)
corrplot(corr_matrix, order = "hclust", method = "color",
         type = 'upper', tl.pos='n')
dev.off()
