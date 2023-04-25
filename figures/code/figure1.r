# Complex Heatmap
packgs <- c("phyloseq", "ComplexHeatmap", "microViz")
lapply(packgs, require, character.only = TRUE)
setwd("~/git/BV_Microbiome/")
rank = "Species"
nfeat = "22"
cl <- list.files(path = "02_cluster/data/C4/", pattern = paste(rank,nfeat,sep = "_"))
vl <- list.files(path = "01_get_valencias/res/")
pseqs_list <- list()

for (i in seq_along(cl)) {
  # load pseq and valencias
  pseq <- readRDS(file = paste0("02_cluster/data/C4/", cl[i]))
  val <- read.csv2(file = paste0("01_get_valencias/res/", vl[i]),
                   header = TRUE, sep =",")
  # check same order
  val <- data.frame(val[,-1], row.names = val[,1])
  identical(rownames(pseq@sam_data), rownames(val))
  # add valencia 
  pseq@sam_data$CST <- val$CST
  # order pseq by target
  df <- pseq@sam_data
  ordered_df <- df[order(df$cluster,decreasing = FALSE), ]
  pseq_sorted <- ps_reorder(ps = pseq, sample_order = rownames(ordered_df))
  pseqs_list[[i]] <- pseq_sorted
}
names(pseqs_list) <- sapply(strsplit(cl, "_"), "[", 1)
# get NS categories for prjna2085
pseqs_list$PRJNA2085@sam_data$Nugent_Score = replace(pseqs_list$PRJNA2085@sam_data$Nugent_Score,
                                                     pseqs_list$PRJNA2085@sam_data$Nugent_Score > 6, "high")
pseqs_list$PRJNA2085@sam_data$Nugent_Score = replace(pseqs_list$PRJNA2085@sam_data$Nugent_Score,
                                                    pseqs_list$PRJNA2085@sam_data$Nugent_Score > 3 & pseqs_list$PRJNA2085@sam_data$Nugent_Score < 7, "intermediate")
pseqs_list$PRJNA2085@sam_data$Nugent_Score = replace(pseqs_list$PRJNA2085@sam_data$Nugent_Score, pseqs_list$PRJNA2085@sam_data$Nugent_Score < 4, "low")

# prepare mat for plotting
mat1 <- t(as.matrix(as.data.frame(pseqs_list$Ravel@otu_table)))
mat2 <- t(as.matrix(as.data.frame(pseqs_list$Sriniv@otu_table)))
mat3 <- t(as.matrix(as.data.frame(pseqs_list$PRJNA2085@otu_table)))
mat4 <- t(as.matrix(as.data.frame(pseqs_list$PRJNA7977@otu_table)))
common_min = min(c(mat1, mat2, mat3, mat4))
common_max = max(c(mat1, mat2, mat3, mat4))
col_fun = circlize::colorRamp2(c(common_min,
                                 ((common_max+common_min)/2),
                                 common_max),transparency = 0.2,
                               c("cornflowerblue","white", "brown3"))
  
# prepare complex heatmap
# Ravel
library(viridisLite)
ann <- data.frame(pseqs_list$Ravel@sam_data$pH, pseqs_list$Ravel@sam_data$target,
                  pseqs_list$Ravel@sam_data$cluster, pseqs_list$Ravel@sam_data$CST)
colnames(ann) <- c('ph', 'Nugent', 'Cluster', 'CST')
colours <- list(
  'ph' = c("4" = "#440154FF","4.3" = "#481B6DFF", "4.4" = "#46337EFF",
            "4.5" = "#3F4889FF", "4.7" = "#365C8DFF", "5" = "#2E6E8EFF",
            "5.3" = "#277F8EFF", "5.4" = "#21908CFF", "5.5" = "#1FA187FF",
            "5.6" = "#2DB27DFF", "5.7" = "#4AC16DFF", "5.8" = "#71CF57FF",
            "6" = "#9FDA3AFF","6.1" = "#CFE11CFF","7" = "#FDE725FF"),
  'Nugent' = c('high' = 'red', 'intermediate' = 'yellow', 'low'= 'green'),
  'Cluster' = c('IDD' = "rosybrown1" ,'N' = 'palegreen1',
                'D' = 'palevioletred1', 'IDN' = 'paleturquoise1'),
  'CST' = c('I'= 'cornflowerblue', 'II' = 'tan4', 'III' = 'yellowgreen',
            'IV-A' = 'darkgray', 'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred',
            'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              'ph' = list(nrow = 1),
                              'Nugent' = list(nrow = 1),
                              'Cluster' = list(nrow = 1),
                              'CST' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            gap = unit(0.5, 'mm'))

h1 <- Heatmap(mat1, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_position = "leftcenter-rot"),
              column_title = "Ravel",
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)


# Sriniv
ann <- data.frame(pseqs_list$Sriniv@sam_data$pH, pseqs_list$Sriniv@sam_data$target,
                  pseqs_list$Sriniv@sam_data$cluster, pseqs_list$Sriniv@sam_data$CST)
colnames(ann) <- c('ph', 'Nugent', 'Cluster', 'CST')
colours <- list(
  'ph' = c("4" = "#440154FF","4.3" = "#481B6DFF", "4.4" = "#46337EFF",
           "4.5" = "#3F4889FF", "4.7" = "#365C8DFF", "5" = "#2E6E8EFF",
           "5.3" = "#277F8EFF", "5.4" = "#21908CFF", "5.5" = "#1FA187FF",
           "5.6" = "#2DB27DFF", "5.7" = "#4AC16DFF", "5.8" = "#71CF57FF",
           "6" = "#9FDA3AFF","6.1" = "#CFE11CFF","7" = "#FDE725FF", " " = "black"),
  'Nugent' = c('high' = 'red', 'intermediate' = 'yellow', 'low'= 'green'),
  'Cluster' = c('IDD' = "rosybrown1" ,'N' = 'palegreen1',
                'D' = 'palevioletred1', 'IDN' = 'paleturquoise1'),
  'CST' = c('I'= 'cornflowerblue', 'II' = 'tan4', 'III' = 'yellowgreen',
            'IV-A' = 'darkgray', 'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred',
            'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              'ph' = list(nrow = 1),
                              'Nugent' = list(nrow = 1),
                              'Cluster' = list(nrow = 1),
                              'CST' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))

h2 <- Heatmap(mat2, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_position = "leftcenter-rot"),
              column_title = "Sriniv",
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

# PR2085
ann <- data.frame(pseqs_list$PRJNA2085@sam_data$pH, pseqs_list$PRJNA2085@sam_data$Nugent_Score,
                  pseqs_list$PRJNA2085@sam_data$cluster, pseqs_list$PRJNA2085@sam_data$CST)
colnames(ann) <- c('ph', 'Nugent', 'Cluster', 'CST')
colours <- list(
  'ph' = c("4" = "#440154FF","4.3" = "#481B6DFF", "4.4" = "#46337EFF",
           "4.5" = "#3F4889FF", "4.7" = "#365C8DFF", "5" = "#2E6E8EFF",
           "5.3" = "#277F8EFF", "5.4" = "#21908CFF", "5.5" = "#1FA187FF",
           "5.6" = "#2DB27DFF", "5.7" = "#4AC16DFF", "5.8" = "#71CF57FF",
           "6" = "#9FDA3AFF","6.1" = "#CFE11CFF","7" = "#FDE725FF", "NA" = "black"),
  'Nugent' = c('high' = 'red', 'intermediate' = 'yellow', 'low'= 'green'),
  'Cluster' = c('IDD' = "rosybrown1" ,'N' = 'palegreen1',
                'D' = 'palevioletred1', 'IDN' = 'paleturquoise1'),
  'CST' = c('I'= 'cornflowerblue', 'II' = 'tan4', 'III' = 'yellowgreen',
            'IV-A' = 'darkgray', 'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred',
            'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              'ph' = list(nrow = 1),
                              'Nugent' = list(nrow = 1),
                              'Cluster' = list(nrow = 1),
                              'CST' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))

h3 <- Heatmap(mat3, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_position = "leftcenter-rot"),
              column_title = "PRJNA2085",
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)


# PR7977
ann <- data.frame(pseqs_list$PRJNA7977@sam_data$cluster, pseqs_list$PRJNA7977@sam_data$CST)
colnames(ann) <- c('Cluster', 'CST')
colours <- list(
  'Cluster' = c('IDD' = "rosybrown1" ,'N' = 'palegreen1',
                'D' = 'palevioletred1', 'IDN' = 'paleturquoise1'),
  'CST' = c('I'= 'cornflowerblue', 'II' = 'tan4', 'III' = 'yellowgreen',
            'IV-A' = 'darkgray', 'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred',
            'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              'Cluster' = list(nrow = 1),
                              'CST' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))

h4 <- Heatmap(mat4, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_position = "leftcenter-rot"),
              column_title = "PRJNA7977",
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

hlist <- h1 + h2 + h3 + h4 
h <- draw(object = hlist, heatmap_legend_side = "right", 
          annotation_legend_side = "bottom")

# METANANLYSIS??