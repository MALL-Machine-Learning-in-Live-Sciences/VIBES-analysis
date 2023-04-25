packgs <- c("phyloseq", "ComplexHeatmap", "microViz")
lapply(packgs, require, character.only = TRUE)
setwd("~/git/BV_Microbiome/")
cl <- readRDS("02_cluster/data/C4/PRJNA2085_Cluster_Species_22_pseq.rds")
vl <- read.csv2(file  = "01_get_valencias/res/PRJNA2085.csv", sep = ",", header = TRUE)
# check same order
vl <- data.frame(vl[,-1], row.names = vl[,1])
identical(rownames(cl@sam_data), rownames(vl))
# add valencia 
cl@sam_data$CST <- vl$CST
# order pseq by target
df <- cl@sam_data
ordered_df <- df[order(df$cluster, decreasing = FALSE),]
pseq_sorted <- ps_reorder(ps = cl, sample_order = rownames(ordered_df))

# prepare mat for plotting
mat <- t(as.matrix(as.data.frame(pseq_sorted@otu_table)))
common_min = min(mat)
common_max = max(mat)
col_fun = circlize::colorRamp2(c(common_min,
                                 ((common_max+common_min)/2),
                                 common_max),transparency = 0.2,
                               c("cornflowerblue","white", "brown3"))
# prepare complex heatmap
ann <- data.frame(pseq_sorted@sam_data$SBV, pseq_sorted@sam_data$ABV,
                  pseq_sorted@sam_data$cluster, pseq_sorted@sam_data$CST)

colnames(ann) <- c('SBV','ABV', 'Cluster', 'CST')
colours <- list(
  'SBV'= c('0' = 'green', '1' = 'red'),
  'ABV'= c('0' = 'green', '1' = 'red'),
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
                              'SBV' = list(nrow = 1),
                              'ABV' = list(nrow = 1),
                              'Cluster' = list(nrow = 1),
                              'CST' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "left",
                            gap = unit(0.5, 'mm'))

h1 <- Heatmap(mat, name = "CLR Species Abundances",
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