packgs <- c("phyloseq", "ComplexHeatmap", "microViz")
lapply(packgs, require, character.only = TRUE)
setwd("~/git/BV_Microbiome/")
cl <- readRDS("02_cluster/data/C4/PRJNA3020_Cluster_Species_22_pseq.rds")
vl <- read.csv2(file  = "01_get_valencias/res/PRJNA3020.csv", sep = ",", header = TRUE)
# check same order
vl <- data.frame(vl[,-1], row.names = vl[,1])
identical(rownames(cl@sam_data), rownames(vl))
# add valencia 
cl@sam_data$CST <- vl$CST
# order pseq by target
df <- cl@sam_data
ordered_df <- df[order(df$CST, decreasing = FALSE),]
pseq_sorted <- ps_reorder(ps = cl, sample_order = rownames(ordered_df))

#D0
meta <- data.frame(sample_data(pseq_sorted))
d0 <- meta[grep("D0", meta$sample_alias), ]
ps0 <- subset_samples(pseq_sorted,(sample_names(pseq_sorted) %in% rownames(d0)))
#D7
d7 <- meta[grep("D7", meta$sample_alias), ]
ps7 <- subset_samples(pseq_sorted,(sample_names(pseq_sorted) %in% rownames(d7)))
#D30
d30 <- meta[grep("D30", meta$sample_alias), ]
ps30 <- subset_samples(pseq_sorted,(sample_names(pseq_sorted) %in% rownames(d30)))

# prepare mat for plotting
mat1 <- t(as.matrix(as.data.frame(ps0@otu_table)))
mat2 <- t(as.matrix(as.data.frame(ps7@otu_table)))
mat3 <- t(as.matrix(as.data.frame(ps30@otu_table)))
common_min = min(c(mat1, mat2, mat3))
common_max = max(c(mat1, mat2, mat3))
col_fun = circlize::colorRamp2(c(common_min,
                                 ((common_max+common_min)/2),
                                 common_max),transparency = 0.2,
                               c("cornflowerblue","white", "brown3"))
# prepare complex heatmap
# D0
ann <- data.frame(ps0@sam_data$status, ps0@sam_data$cluster,
                  ps0@sam_data$CST)

colnames(ann) <- c('Status', 'Cluster', 'CST')
colours <- list(
  'Status' = c('Cured' = 'green', 'Failed' = 'red'),
  'Cluster' = c('IDD' = "rosybrown1", 'D' = 'palevioletred1',
                'IDN' = 'paleturquoise1'),
  'CST' = c('III' = 'yellowgreen', 'IV-A' = 'darkgray',
            'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              'Status' = list(nrow = 1),
                              'Cluster' = list(nrow = 1),
                              'CST' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            gap = unit(0.5, 'mm'))

h0 <- Heatmap(mat1, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_position = "leftcenter-rot"),
              column_title = "Pre-treatment",
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

#D7
ann <- data.frame(ps7@sam_data$status, ps7@sam_data$cluster,
                  ps7@sam_data$CST)

colnames(ann) <- c('Status', 'Cluster', 'CST')
colours <- list(
  'Status' = c('Cured' = 'green', 'Failed' = 'red'),
  'Cluster' = c('IDD' = "rosybrown1", 'D' = 'palevioletred1',
                'IDN' = 'paleturquoise1', 'N' = 'palegreen1'),
  'CST' = c('I'= 'cornflowerblue', 'III' = 'yellowgreen','IV-B' = 'seagreen',
            'IV-C' = 'mediumvioletred', 'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              'Status' = list(nrow = 1),
                              'Cluster' = list(nrow = 1),
                              'CST' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))

h7 <- Heatmap(mat2, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_position = "leftcenter-rot"),
              column_title = "After one week",
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

#D30
ann <- data.frame(ps30@sam_data$status, ps30@sam_data$cluster,
                  ps30@sam_data$CST)

colnames(ann) <- c('Status', 'Cluster', 'CST')
colours <- list(
  'Status' = c('Cured' = 'green', 'Failed' = 'red'),
  'Cluster' = c('IDD' = "rosybrown1", 'D' = 'palevioletred1',
                'IDN' = 'paleturquoise1', 'N' = 'palegreen1'),
  'CST' = c('I'= 'cornflowerblue', 'III' = 'yellowgreen','IV-A' = 'darkgray',
            'IV-B' = 'seagreen','IV-C' = 'mediumvioletred')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              'Status' = list(nrow = 1),
                              'Cluster' = list(nrow = 1),
                              'CST' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))

h30 <- Heatmap(mat3, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_position = "leftcenter-rot"),
              column_title = "After one month",
              row_title = "Species",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

hlist <- h0 + h7 + h30
h <- draw(object = hlist, heatmap_legend_side = "right", 
          annotation_legend_side = "bottom")
