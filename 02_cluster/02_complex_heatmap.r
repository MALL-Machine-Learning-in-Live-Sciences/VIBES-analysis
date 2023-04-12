# Complex Heatmap
packgs <- c("phyloseq", "ComplexHeatmap")
lapply(packgs, require, character.only = TRUE)
setwd("~/git/BV_Microbiome/")
ranks <- c(3, 4)
rav_pseqs <- list()
sri_pseqs <- list()
prj_pseqs <- list()
prj79_pseqs <- list()
prj20_pseqs <- list()
# Bucle for load pseqs for clustrs 3 and 4
for (i in seq_along(ranks)) {
  rav_pseqs[[i]] <- readRDS(paste0("extdata/SpeciesIntersect/Cluster/C",ranks[i],"/",
                             "Ravel_Cluster_Species_22_pseq.rds"))
  sri_pseqs[[i]] <- readRDS(paste0("extdata/SpeciesIntersect/Cluster/C",ranks[i],"/",
                             "Sriniv_Cluster_Species_22_pseq.rds"))
  prj_pseqs[[i]] <- readRDS(paste0("extdata/SpeciesIntersect/Cluster/C",ranks[i],"/",
                             "PRJNA3020_Cluster_Species_22_pseq.rds"))
  prj79_pseqs[[i]] <- readRDS(paste0("extdata/SpeciesIntersect/Cluster/C",ranks[i],"/",
                                   "PRJNA7977D0_Cluster_Species_22_pseq.rds"))
  prj20_pseqs[[i]] <- readRDS(paste0("extdata/SpeciesIntersect/Cluster/C",ranks[i],"/",
                                     "PR2085D0_Cluster_Species_22_pseq.rds"))
}
# Rename lists
names(rav_pseqs) <- c("C3", "C4")
names(sri_pseqs) <- c("C3", "C4")
names(prj_pseqs) <- c("C3", "C4")
names(prj79_pseqs) <- c("C3", "C4")
names(prj20_pseqs) <- c("C3", "C4")
# Adding on metadata 4 clusters
# ravel
identical(rownames(rav_pseqs[["C3"]]@sam_data), rownames(rav_pseqs[["C4"]]@sam_data))
names(rav_pseqs[["C3"]]@sam_data)[names(rav_pseqs[["C3"]]@sam_data) == "cluster"] <- "Cluster_3"
rav_pseqs[["C3"]]@sam_data$Cluster_4 <- rav_pseqs[["C4"]]@sam_data$cluster
# sri
identical(rownames(sri_pseqs[["C3"]]@sam_data), rownames(sri_pseqs[["C4"]]@sam_data))
names(sri_pseqs[["C3"]]@sam_data)[names(sri_pseqs[["C3"]]@sam_data) == "cluster"] <- "Cluster_3"
sri_pseqs[["C3"]]@sam_data$Cluster_4 <- sri_pseqs[["C4"]]@sam_data$cluster
#prj
identical(rownames(prj_pseqs[["C3"]]@sam_data), rownames(prj_pseqs[["C4"]]@sam_data))
names(prj_pseqs[["C3"]]@sam_data)[names(prj_pseqs[["C3"]]@sam_data) == "cluster"] <- "Cluster_3"
prj_pseqs[["C3"]]@sam_data$Cluster_4 <- prj_pseqs[["C4"]]@sam_data$cluster
#prj79
identical(rownames(prj79_pseqs[["C3"]]@sam_data), rownames(prj79_pseqs[["C4"]]@sam_data))
names(prj79_pseqs[["C3"]]@sam_data)[names(prj79_pseqs[["C3"]]@sam_data) == "cluster"] <- "Cluster_3"
prj79_pseqs[["C3"]]@sam_data$Cluster_4 <- prj79_pseqs[["C4"]]@sam_data$cluster
#prj20
identical(rownames(prj20_pseqs[["C3"]]@sam_data), rownames(prj20_pseqs[["C4"]]@sam_data))
names(prj20_pseqs[["C3"]]@sam_data)[names(prj20_pseqs[["C3"]]@sam_data) == "cluster"] <- "Cluster_3"
prj20_pseqs[["C3"]]@sam_data$Cluster_4 <- prj20_pseqs[["C4"]]@sam_data$cluster

#Aglomerate al pseqs in a list
pseqs <- list(rav_pseqs$C3, sri_pseqs$C3, prj_pseqs$C3, prj79_pseqs$C3, prj20_pseqs$C3)
names(pseqs) <- c("ravel", "srinivasan", "prjna3020", "prjna7977", "prjna2085")

# Ordenar por clusters 
require(microViz)
#ravel
dfr <- pseqs$ravel@sam_data
ordered_dfr <- dfr[order(dfr$Cluster_4), ]
rav_sorted <- ps_reorder(ps = pseqs[["ravel"]], sample_order = rownames(ordered_dfr))
#sriniv
dfs <- pseqs$srinivasan@sam_data
ordered_dfs <- dfs[order(dfs$Cluster_4), ]
sri_sorted <- ps_reorder(ps = pseqs[["srinivasan"]], sample_order = rownames(ordered_dfs))
#prjna3020
dfp <- pseqs$prjna3020@sam_data
ordered_dfp <- dfp[order(dfp$Cluster_4), ]
prj_sorted <- ps_reorder(ps = pseqs[["prjna3020"]], sample_order = rownames(ordered_dfp))
#prjna7977
dfp79 <- pseqs$prjna7977@sam_data
ordered_dfp79 <- dfp79[order(dfp79$Cluster_4), ]
prj79_sorted <- ps_reorder(ps = pseqs[["prjna7977"]], sample_order = rownames(ordered_dfp79))
#prjna2080
dfp20 <- pseqs$prjna2080@sam_data
ordered_dfp20<- dfp20[order(dfp20$Cluster_4), ]
prj20_sorted <- ps_reorder(ps = pseqs[["prjna2085"]], sample_order = rownames(ordered_dfp20))


# ComplexHeatmap
#Species
mat1 <- t(as.matrix(as.data.frame(rav_sorted@otu_table)))
mat2 <- t(as.matrix(as.data.frame(sri_sorted@otu_table)))
mat3 <- t(as.matrix(as.data.frame(prj_sorted@otu_table)))
mat4 <- t(as.matrix(as.data.frame(prj79_sorted@otu_table)))
mat5 <- t(as.matrix(as.data.frame(prj20_sorted@otu_table)))
common_min = min(c(mat1, mat2, mat3, mat4, mat5))
common_max = max(c(mat1, mat2, mat3, mat4, mat5))
col_fun = circlize::colorRamp2(c(common_min,((common_max+common_min)/2), common_max),transparency = 0.2, c("cornflowerblue","white", "brown3"))

#Ravel
ann <- data.frame(rav_sorted@sam_data$Cluster_4, rav_sorted@sam_data$Cluster_3)
colnames(ann) <- c('4 Cluster', '3 Cluster')
colours <- list('4 Cluster' = c('IDD' = "rosybrown1" ,'N' = 'palegreen1',
                                'D' = 'palevioletred1', 'IDN' = 'paleturquoise1'),
                '3 Cluster' = c('D' = "red1" ,'N' = 'royalblue1',
                                'ID' = 'orchid1'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              '4 Cluster' = list(nrow = 1),
                              '3 Cluster' = list(nrow = 1)),
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
ann <- data.frame(sri_sorted@sam_data$Cluster_4, sri_sorted@sam_data$Cluster_3)
colnames(ann) <- c('4 Cluster', '3 Cluster')
colours <- list('4 Cluster' = c('IDD' = "rosybrown1" ,'N' = 'palegreen1',
                                'D' = 'palevioletred1', 'IDN' = 'paleturquoise1'),
                '3 Cluster' = c('D' = "red1" ,'N' = 'royalblue1',
                                'ID' = 'orchid1'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_legend_param = list(
                              '4 Cluster' = list(nrow = 1),
                              '3 Cluster' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "right",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))

h2 <- Heatmap(mat2,name = "CLR Species Abundances",
              column_title = "Sriniv",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

# PRJNA3020
ann <- data.frame(prj_sorted@sam_data$Cluster_4, prj_sorted@sam_data$Cluster_3)
colnames(ann) <- c('4 Cluster', '3 Cluster')
colours <- list('4 Cluster' = c('IDD' = "rosybrown1" ,'N' = 'palegreen1',
                                'D' = 'palevioletred1', 'IDN' = 'paleturquoise1'),
                '3 Cluster' = c('D' = "red1" ,'N' = 'royalblue1',
                                'ID' = 'orchid1'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              '4 Cluster' = list(nrow = 1),
                              '3 Cluster' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "right",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))
h3 <- Heatmap(mat3,name = "CLR Species Abundances",
              column_title = "P3020",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

# PRJNA7977
ann <- data.frame(prj79_sorted@sam_data$Cluster_4, prj79_sorted@sam_data$Cluster_3)
colnames(ann) <- c('4 Cluster', '3 Cluster')
colours <- list('4 Cluster' = c('IDD' = "rosybrown1" ,'N' = 'palegreen1',
                                'D' = 'palevioletred1', 'IDN' = 'paleturquoise1'),
                '3 Cluster' = c('D' = "red1" ,'N' = 'royalblue1',
                                'ID' = 'orchid1'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              '4 Cluster' = list(nrow = 1),
                              '3 Cluster' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "right",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))
h4 <- Heatmap(mat4,name = "CLR Species Abundances",
              column_title = "P7977",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)
# PRJNA2080
ann <- data.frame(prj20_sorted@sam_data$Cluster_4, prj20_sorted@sam_data$Cluster_3)
colnames(ann) <- c('4 Cluster', '3 Cluster')
colours <- list('4 Cluster' = c('IDD' = "rosybrown1" ,'N' = 'palegreen1',
                                'D' = 'palevioletred1', 'IDN' = 'paleturquoise1'),
                '3 Cluster' = c('D' = "red1" ,'N' = 'royalblue1',
                                'ID' = 'orchid1'))
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 9),
                            annotation_legend_param = list(
                              '4 Cluster' = list(nrow = 1),
                              '3 Cluster' = list(nrow = 1)),
                            height = unit(4.5, 'mm'),
                            annotation_name_side = "right",
                            show_annotation_name = FALSE,
                            gap = unit(0.5, 'mm'))
h5 <- Heatmap(mat5,name = "CLR Species Abundances",
              column_title = "P2085",
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)


hlist <- h1 + h2 + h3 + h4 +h5
h <- draw(object = hlist, heatmap_legend_side = "right", 
          annotation_legend_side = "bottom")


# Genera
# Ravel
ravel_intersect <- readRDS(file = "extdata/SpeciesIntersect/Ravel_Species_pseq_22.rds")
ravel_intersect <- tax_glom(physeq = ravel_intersect, taxrank = "Genus")
ravel_clr <-  microbiome::transform(ravel_intersect, transform =  'clr')
t_names <- data.frame(tax_table(ravel_clr))$Genus
taxa_names(ravel_clr) <- t_names
ravel_clr <- ps_reorder(ps = ravel_clr, sample_order = rownames(ordered_dfr))
#Sriniv
sriniv_intersect <- readRDS(file = "extdata/SpeciesIntersect/Sriniv_Species_pseq_22.rds")
sriniv_intersect <- tax_glom(physeq = sriniv_intersect, taxrank = "Genus")
sriniv_clr <- microbiome::transform(sriniv_intersect, transform =  'clr')
taxa_names(sriniv_clr) <- t_names
sriniv_clr <- ps_reorder(ps = sriniv_clr, sample_order = rownames(ordered_dfs))
#PR3020
prj3020_intersect <- readRDS(file = "extdata/SpeciesIntersect/PRJNA3020_Species_pseq_22.rds")
prj3020_intersect <- tax_glom(physeq = prj3020_intersect, taxrank = "Genus")
prj3020_clr <- microbiome::transform(prj3020_intersect, transform =  'clr')
taxa_names(prj3020_clr) <- t_names
prj3020_clr <- ps_reorder(ps = prj3020_clr, sample_order = rownames(ordered_dfp))
#PR7977
prj7977_intersect <- readRDS(file = "extdata/SpeciesIntersect/PRJNA7977D0_Species_pseq_22.rds")
prj7977_intersect <- tax_glom(physeq = prj7977_intersect, taxrank = "Genus")
prj7977_clr <- microbiome::transform(prj7977_intersect, transform =  'clr')
taxa_names(prj7977_clr) <- t_names
prj7977_clr <- ps_reorder(ps = prj7977_clr, sample_order = rownames(ordered_dfp79))
#PR72085
prj2085_intersect <- readRDS(file = "extdata/SpeciesIntersect/PRJNA2085D0_Species_pseq_22.rds")
prj2085_intersect <- tax_glom(physeq = prj2085_intersect, taxrank = "Genus")
prj2085_clr <- microbiome::transform(prj2085_intersect, transform =  'clr')
taxa_names(prj2085_clr) <- t_names
prj2085_clr <- ps_reorder(ps = prj2085_clr, sample_order = rownames(ordered_dfp20))

mat1 <- t(as.matrix(as.data.frame(ravel_clr@otu_table)))
mat2 <- t(as.matrix(as.data.frame(sriniv_clr@otu_table)))
mat3 <- t(as.matrix(as.data.frame(prj3020_clr@otu_table)))
mat4 <- t(as.matrix(as.data.frame(prj7977_clr@otu_table)))
mat5 <- t(as.matrix(as.data.frame(prj2085_clr@otu_table)))

common_min = min(c(mat1, mat2, mat3, mat4, mat5))
common_max = max(c(mat1, mat2, mat3, mat4, mat5))
col_fun = circlize::colorRamp2(c(common_min,((common_max+common_min)/2), common_max),transparency = 0.4, c("darkblue","white", "darkorange2"))

h1 <- Heatmap(mat1, name = "CLR Genera Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_position = "leftcenter-rot"),
              row_title = "Genera",
              col = col_fun,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

h2 <- Heatmap(mat2,name = "CLR Genera Abundances",
              col = col_fun,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

h3 <- Heatmap(mat3,name = "CLR Genera Abundances",
              col = col_fun,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

h4 <- Heatmap(mat4,name = "CLR Genera Abundances",
              col = col_fun,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

h5 <- Heatmap(mat5,name = "CLR Genera Abundances",
              col = col_fun,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 8),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

hlist <- h1 + h2 + h3 + h4 + h5
h <- draw(object = hlist, heatmap_legend_side = "right", 
          annotation_legend_side = "bottom")
