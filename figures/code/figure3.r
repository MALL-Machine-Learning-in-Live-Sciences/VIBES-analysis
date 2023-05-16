require(pheatmap)
library(cowplot)
library(magrittr)
require(MOFA2)
require(ggplot2)
library(RColorBrewer)
require(viridis)
require(microViz)
require(ggpubr)

# Load objects
# =====
# MEFISTO
sm <- MOFA2::load_model("04_treatment/res/PRJNA302078.rds")
feat <- readRDS("04_treatment/data/PRJNA3020_featdata.rds")
identical(feat$feature, features_metadata(sm)$feature)
feat$feature <- feat$Species
features_metadata(sm) <- feat


# Figure 3a
# ==============
p1 <- plot_variance_explained(
  sm, 
  plot_total = T,
  x = "group", 
  split_by = "view")
  
p1 <- p1[[1]] +
  theme( legend.title = element_text(size = 8),
         legend.text = element_text(size = 8),
         strip.text.x = element_text(size = 8),
         axis.text.y = element_text(size = 8),
         axis.title.y = element_blank(),
         axis.title.x = element_text(size = 8),
         axis.text.x = element_blank()) + 
    scale_fill_gradient(
      low="white",
      high=viridis(1)) +
  xlab("samples")
  

ggsave(
  p1,
  filename = "fig3a.svg",
  device = "svg",
  path = "figures/plots/Figure3/",
  width = 65, 
  height = 35, 
  units = "mm")


# Figure 3b
# ===============
p2 <- plot_factors_vs_cov(
  factors = c(1, 2),
  sm,
  dot_size = 3,
  color_by = "status",
  alpha = .2, 
  scale = T) +
  theme(legend.position = "none",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 9)) +
  scale_fill_manual(values = c(viridis(5)[1], viridis(5)[5])) +
  scale_y_continuous(position = "right") +
  stat_summary(aes(col = color_by),
               geom = "line",
               size = 2,
               fun = "mean",
               show.legend = TRUE) +
  scale_colour_manual(values = c(viridis(5)[1], viridis(5)[5]))

ggsave(
  p2,
  filename = "fig3b.svg",
  device = "svg",
  path = "figures/plots/Figure3/",
  width = 100, 
  height = 45, 
  units = "mm")



# Figure 3c
# =========
rownames(sm@expectations[["W"]][["microbiome"]]) <- feat$Species
p3 <- plot_top_weights(
  sm,
  factors = c(1, 2),
  scale = TRUE) + 
  coord_flip() +
  theme(axis.title.x = element_blank(),
    strip.background=element_rect(fill="white"),
    strip.text.x = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8))

ggsave(
  p3,
  filename = "fig3c.svg",
  device = "svg",
  path = "figures/plots/Figure3/",
  width = 175, 
  height = 42, 
  units = "mm")


# Figure 3d
# ==============
ans = 7
rns = 8
ts = 7
pseq_sorted <- readRDS(file = "figures/data/f3_data.rds")

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
ann <- data.frame(ps0@sam_data$status, ps0@sam_data$p_cluster,
                  ps0@sam_data$CST)

colnames(ann) <- c('Status', 'Cluster', 'CST')
colours <- list(
  'Status' = c('Cured' = 'green', 'Failed' = 'red'),
  'Cluster' = c("VCS-III" = "#A6CEE3" , "VCS-I" = '#33A02C',
                "VCS-IV" = '#1F78B4', "VCS-II" = '#B2DF8A'),
  'CST' = c('III' = 'yellowgreen', 'IV-A' = 'darkgray',
            'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            simple_anno_size = unit(0.3, "cm"),
                            annotation_name_gp = gpar(fontsize = ans),
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
              column_names_gp = grid::gpar(fontsize = ts),
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = rns),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

#D7
ann <- data.frame(ps7@sam_data$status, ps7@sam_data$p_cluster,
                  ps7@sam_data$CST)

colnames(ann) <- c('Status', 'Cluster', 'CST')
colours <- list(
  'Status' = c('Cured' = 'green', 'Failed' = 'red'),
  'Cluster' = c("VCS-III" = "#A6CEE3" , "VCS-I" = '#33A02C',
                "VCS-IV" = '#1F78B4', "VCS-II" = '#B2DF8A'),
  'CST' = c('I'= 'cornflowerblue', 'III' = 'yellowgreen','IV-B' = 'seagreen',
            'IV-C' = 'mediumvioletred', 'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            simple_anno_size = unit(0.3, "cm"),
                            annotation_name_gp = gpar(fontsize = ans),
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
              column_names_gp = grid::gpar(fontsize = ts),
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = rns),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

#D30
ann <- data.frame(ps30@sam_data$status, ps30@sam_data$p_cluster,
                  ps30@sam_data$CST)

colnames(ann) <- c('Status', 'Cluster', 'CST')
colours <- list(
  'Status' = c('Cured' = 'green', 'Failed' = 'red'),
  'Cluster' = c("VCS-III" = "#A6CEE3" , "VCS-I" = '#33A02C',
                "VCS-IV" = '#1F78B4', "VCS-II" = '#B2DF8A'),
  'CST' = c('I'= 'cornflowerblue', 'III' = 'yellowgreen','IV-A' = 'darkgray',
            'IV-B' = 'seagreen','IV-C' = 'mediumvioletred')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            simple_anno_size = unit(0.3, "cm"),
                            annotation_name_gp = gpar(fontsize = ans),
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
              column_names_gp = grid::gpar(fontsize = ts),
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = rns),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

hlist <- h0 + h7 + h30

svg(filename = "figures/plots/Figure3/fig3d.svg",
    width = 6.89,
    height = 4,
    pointsize = 7)

draw(object = hlist, heatmap_legend_side = "right", 
     annotation_legend_side = "bottom")
dev.off()


# Figure 3e
# ============

rs = readRDS("04_treatment/res/summary_performances.rds")

p4 <- ggscatter(
  data = rs,
  palette = c("#440154FF", "#31688EFF", "#35B779FF", "#EB8055FF"),
  x = "value",
  y = "transformation",
  color = "algorithm",
  shape = "algorithm",
  size = 5
) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        axis.text.y = element_text(color = "black",
                                   size = 8),
        strip.text.x = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  )

ggsave(
  p4,
  filename = "fig3e.svg",
  device = "svg",
  path = "figures/plots/Figure3/",
  width = 170, 
  height = 40, 
  units = "mm")

