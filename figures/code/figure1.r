# Complex Heatmap
# Load packages
packgs <- c("phyloseq", "ComplexHeatmap", "microViz", "viridis")
lapply(packgs, require, character.only = TRUE)
# Load data
pseqs_list <- readRDS(file = "figures/data/f1_data.rds")
# declare variables
# Annotation text size
ans = 7
# Rows text size
rns = 7
# Title size
ts = 9
# Legend title size
lts = 8
# legend text size
wls = 7

# prepare data2plot
# prepare mat for plotting
mat1 <- t(as.matrix(as.data.frame(pseqs_list$Ravel@otu_table)))
mat2 <- t(as.matrix(as.data.frame(pseqs_list$Sriniv@otu_table)))
mat3 <- t(as.matrix(as.data.frame(pseqs_list$PRJNA2085@otu_table)))
mat4 <- t(as.matrix(as.data.frame(pseqs_list$PRJNA7977@otu_table)))
mat5 <- t(as.matrix(as.data.frame(pseqs_list$PRJNA3020@otu_table)))
common_min = min(c(mat1, mat2, mat3, mat4, mat5))
common_max = max(c(mat1, mat2, mat3, mat4, mat5))
col_fun = circlize::colorRamp2(c(common_min,
                                 ((common_max+common_min)/2),
                                 common_max),
                               c("cornflowerblue","white", "brown3"))

# prepare pH for plotting
ph1 <- t(as.matrix(as.data.frame(as.numeric(pseqs_list$Ravel@sam_data$pH))))
ph2 <- t(as.matrix(as.data.frame(as.numeric(pseqs_list$Sriniv@sam_data$pH))))
ph3 <- t(as.matrix(as.data.frame(as.numeric(pseqs_list$PRJNA2085@sam_data$pH))))
ph_min = min(c(ph1, ph2, ph3), na.rm=T)
ph_max = max(c(ph1, ph2, ph3), na.rm=T)
ph_fun = circlize::colorRamp2(c(ph_min,
                                ((ph_max+ph_min)/2),
                                ph_max),transparency = 0.2,
                              viridis(3))

# prepare complex heatmap
# Ravel
library(viridisLite)
ann <- data.frame(as.numeric(pseqs_list$Ravel@sam_data$pH),
                  pseqs_list$Ravel@sam_data$target,
                  pseqs_list$Ravel@sam_data$cluster,
                  pseqs_list$Ravel@sam_data$CST)
colnames(ann) <- c('ph', 'Nugent', 'Cluster', 'CST')
colours <- list(
  'ph' = ph_fun,
  'Nugent' = c('high' = 'red', 'intermediate' = 'yellow', 'low'= 'green'),
  'Cluster' = c('VCS-III' = "#A6CEE3" ,'VCS-I' = '#33A02C',
                'VCS-IV' = '#1F78B4', 'VCS-II' = '#B2DF8A'),
  'CST' = c('I'= 'cornflowerblue', 'II' = 'tan4', 'III' = 'yellowgreen',
            'IV-A' = 'darkgray', 'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred',
            'V' = 'tan2')
)

colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            simple_anno_size = unit(0.3, "cm"),
                            annotation_name_gp = gpar(fontsize = ans),
                            annotation_legend_param = list(
                              'ph' = list(nrow = 1,
                                          legend_direction = "horizontal",
                                          title_gp = gpar(fontsize = lts), 
                                          labels_gp = gpar(fontsize = wls)),
                              'Nugent' =list(nrow = 1,
                                             legend_direction = "horizontal",
                                             title_gp = gpar(fontsize = lts), 
                                             labels_gp = gpar(fontsize = wls),
                                             grid_height = unit(2, "mm"), 
                                             grid_width = unit(2, "mm")),
                              'Cluster' = list(nrow = 1,
                                               legend_direction = "horizontal",
                                               title_gp = gpar(fontsize = lts), 
                                               labels_gp = gpar(fontsize = wls),
                                               grid_height = unit(2, "mm"), 
                                               grid_width = unit(2, "mm")),
                              'CST' = list(nrow = 1,
                                           legend_direction = "horizontal",
                                           title_gp = gpar(fontsize = lts), 
                                           labels_gp = gpar(fontsize = wls),
                                           grid_height = unit(2, "mm"), 
                                           grid_width = unit(2, "mm"))),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            gap = unit(0.5, 'mm'))

h1 <- Heatmap(mat1, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_gp = gpar(fontsize = lts),
                                          labels_gp = gpar(fontsize = wls),
                                          title_position = "leftcenter-rot"),
              column_title = "Discovery",
              column_title_gp = grid::gpar(fontsize = ts),
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = rns),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

# draw(object = h1, heatmap_legend_side = "right", annotation_legend_side = "bottom")
# Sriniv
ann <- data.frame(as.numeric(pseqs_list$Sriniv@sam_data$pH),
                  pseqs_list$Sriniv@sam_data$target,
                  pseqs_list$Sriniv@sam_data$cluster,
                  pseqs_list$Sriniv@sam_data$CST)
colnames(ann) <- c('ph', 'Nugent', 'Cluster', 'CST')
colours <- list(
  'ph' = ph_fun,
  'Nugent' = c('high' = 'red', 'intermediate' = 'yellow', 'low'= 'green'),
  'Cluster' = c('VCS-III' = "#A6CEE3" ,'VCS-I' = '#33A02C',
                'VCS-IV' = '#1F78B4', 'VCS-II' = '#B2DF8A'),
  'CST' = c('I'= 'cornflowerblue', 'II' = 'tan4', 'III' = 'yellowgreen',
            'IV-A' = 'darkgray', 'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred',
            'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            simple_anno_size = unit(0.3, "cm"),
                            annotation_name_gp = gpar(fontsize = ans),
                            annotation_legend_param = list(
                              'ph' = list(nrow = 1,
                                          legend_direction = "horizontal",
                                          title_gp = gpar(fontsize = lts), 
                                          labels_gp = gpar(fontsize = wls)),
                              'Nugent' =list(nrow = 1,
                                             legend_direction = "horizontal",
                                             title_gp = gpar(fontsize = lts), 
                                             labels_gp = gpar(fontsize = wls),
                                             grid_height = unit(2, "mm"), 
                                             grid_width = unit(2, "mm")),
                              'Cluster' = list(nrow = 1,
                                               legend_direction = "horizontal",
                                               title_gp = gpar(fontsize = lts), 
                                               labels_gp = gpar(fontsize = wls),
                                               grid_height = unit(2, "mm"), 
                                               grid_width = unit(2, "mm")),
                              'CST' = list(nrow = 1,
                                           legend_direction = "horizontal",
                                           title_gp = gpar(fontsize = lts), 
                                           labels_gp = gpar(fontsize = wls),
                                           grid_height = unit(2, "mm"), 
                                           grid_width = unit(2, "mm"))),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))

h2 <- Heatmap(mat2, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_gp = gpar(fontsize = lts),
                                          labels_gp = gpar(fontsize = wls),
                                          title_position = "leftcenter-rot"),
              column_title = "Val. 1",
              column_title_gp = grid::gpar(fontsize = ts),
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = rns),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)
h2
# PR2085
ann <- data.frame(as.numeric(pseqs_list$PRJNA2085@sam_data$pH), pseqs_list$PRJNA2085@sam_data$Nugent_Score,
                  pseqs_list$PRJNA2085@sam_data$cluster, pseqs_list$PRJNA2085@sam_data$CST)
colnames(ann) <- c('ph', 'Nugent', 'Cluster', 'CST')
colours <- list(
  'ph' = ph_fun,
  'Nugent' = c('high' = 'red', 'intermediate' = 'yellow', 'low'= 'green'),
  'Cluster' = c('VCS-III' = "#A6CEE3" ,'VCS-I' = '#33A02C',
                'VCS-IV' = '#1F78B4', 'VCS-II' = '#B2DF8A'),
  'CST' = c('I'= 'cornflowerblue', 'II' = 'tan4', 'III' = 'yellowgreen',
            'IV-A' = 'darkgray', 'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred',
            'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            simple_anno_size = unit(0.3, "cm"),
                            annotation_name_gp = gpar(fontsize = ans),
                            annotation_legend_param = list(
                              'ph' = list(nrow = 1,
                                          legend_direction = "horizontal",
                                          title_gp = gpar(fontsize = lts), 
                                          labels_gp = gpar(fontsize = wls)),
                              'Nugent' =list(nrow = 1,
                                             legend_direction = "horizontal",
                                             title_gp = gpar(fontsize = lts), 
                                             labels_gp = gpar(fontsize = wls),
                                             grid_height = unit(2, "mm"), 
                                             grid_width = unit(2, "mm")),
                              'Cluster' = list(nrow = 1,
                                               legend_direction = "horizontal",
                                               title_gp = gpar(fontsize = lts), 
                                               labels_gp = gpar(fontsize = wls),
                                               grid_height = unit(2, "mm"), 
                                               grid_width = unit(2, "mm")),
                              'CST' = list(nrow = 1,
                                           legend_direction = "horizontal",
                                           title_gp = gpar(fontsize = lts), 
                                           labels_gp = gpar(fontsize = wls),
                                           grid_height = unit(2, "mm"), 
                                           grid_width = unit(2, "mm"))),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))

h3 <- Heatmap(mat3, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_gp = gpar(fontsize = lts),
                                          labels_gp = gpar(fontsize = wls),
                                          title_position = "leftcenter-rot"),
              column_title = "Val. 2",
              column_title_gp = grid::gpar(fontsize = ts),
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = rns),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)


# PR7977
ann <- data.frame(pseqs_list$PRJNA7977@sam_data$cluster, pseqs_list$PRJNA7977@sam_data$CST)
colnames(ann) <- c('Cluster', 'CST')
colours <- list(
  'Cluster' = c('VCS-III' = "#A6CEE3" ,'VCS-I' = '#33A02C',
                'VCS-IV' = '#1F78B4', 'VCS-II' = '#B2DF8A'),
  'CST' = c('I'= 'cornflowerblue', 'II' = 'tan4', 'III' = 'yellowgreen',
            'IV-A' = 'darkgray', 'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred',
            'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            simple_anno_size = unit(0.3, "cm"),
                            annotation_name_gp = gpar(fontsize = ans),
                            annotation_legend_param = list(
                              'Cluster' = list(nrow = 1,
                                               legend_direction = "horizontal",
                                               title_gp = gpar(fontsize = lts), 
                                               labels_gp = gpar(fontsize = wls),
                                               grid_height = unit(2, "mm"), 
                                               grid_width = unit(2, "mm")),
                              'CST' = list(nrow = 1,
                                           legend_direction = "horizontal",
                                           title_gp = gpar(fontsize = lts), 
                                           labels_gp = gpar(fontsize = wls),
                                           grid_height = unit(2, "mm"), 
                                           grid_width = unit(2, "mm"))),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))

h4 <- Heatmap(mat4, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_gp = gpar(fontsize = lts),
                                          labels_gp = gpar(fontsize = wls),
                                          title_position = "leftcenter-rot"),
              column_title = "Val. 3",
              column_title_gp = grid::gpar(fontsize = ts),
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = rns),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

# PR3020
ann <- data.frame(pseqs_list$PRJNA3020@sam_data$cluster, pseqs_list$PRJNA3020@sam_data$CST)
colnames(ann) <- c('Cluster', 'CST')
colours <- list(
  'Cluster' = c('VCS-III' = "#A6CEE3" ,'VCS-I' = '#33A02C',
                'VCS-IV' = '#1F78B4', 'VCS-II' = '#B2DF8A'),
  'CST' = c('I'= 'cornflowerblue', 'III' = 'yellowgreen',
            'IV-A' = 'darkgray', 'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred',
            'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            simple_anno_size = unit(0.3, "cm"),
                            annotation_name_gp = gpar(fontsize = ans),
                            annotation_legend_param = list(
                              'Cluster' = list(nrow = 1,
                                               legend_direction = "horizontal",
                                               title_gp = gpar(fontsize = lts), 
                                               labels_gp = gpar(fontsize = wls),
                                               grid_height = unit(2, "mm"), 
                                               grid_width = unit(2, "mm")),
                              'CST' = list(nrow = 1,
                                           legend_direction = "horizontal",
                                           title_gp = gpar(fontsize = lts), 
                                           labels_gp = gpar(fontsize = wls),
                                           grid_height = unit(2, "mm"), 
                                           grid_width = unit(2, "mm"))),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))

h5 <- Heatmap(mat5, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_gp = gpar(fontsize = lts),
                                          labels_gp = gpar(fontsize = wls),
                                          title_position = "leftcenter-rot"),
              column_title = "Val. 4",
              column_title_gp = grid::gpar(fontsize = ts),
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = rns),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

hlist <- h1 + h2 + h3 + h4 + h5

svg(
  filename = "figures/plots/Figure1/fig1.svg",
  width = 7.08661, 
  height = 3.93701,
  )

draw(object = hlist, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
