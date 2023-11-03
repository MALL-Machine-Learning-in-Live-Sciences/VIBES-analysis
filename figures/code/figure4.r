require(pheatmap)
require(ComplexHeatmap)
library(cowplot)
library(magrittr)
require(MOFA2)
require(ggplot2)
library(RColorBrewer)
require(viridis)
require(microViz)
require(ggpubr)
require(phyloseq)
# Load objects
# =====
# MEFISTO
sm <- MOFA2::load_model("04_treatment/res/PRJNA302078.rds")
feat <- readRDS("04_treatment/data/PRJNA3020_featdata.rds")
identical(feat$feature, features_metadata(sm)$feature)
feat$feature <- feat$Species
features_metadata(sm) <- feat


# Figure 4a
# ===============
p2 <- plot_factors_vs_cov(
  factors = c(1, 2),
  sm,
  dot_size = 1,
  color_by = "status",
  alpha = .2, 
  scale = T) +
  theme(legend.position = "none",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        strip.text.x = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 8)) +
  scale_fill_manual(values = c(viridis(5)[1], viridis(5)[5])) +
  scale_y_continuous(position = "right") +
  stat_summary(aes(col = color_by),
               geom = "line",
               size = 1,
               fun = "mean",
               show.legend = TRUE) +
  scale_colour_manual(values = c(viridis(5)[1], viridis(5)[5]))

ggsave(
  p2,
  filename = "fig4a.svg",
  device = "svg",
  path = "figures/plots/Figure4/",
  width = 60, 
  height = 45, 
  units = "mm")



# Figure 4b
# =========
rownames(sm@expectations[["W"]][["microbiome"]]) <- feat$Species
p3 <- plot_top_weights(
  sm,
  factors = c(1, 2),
  scale = TRUE) + 
  geom_point(size=1)+
  coord_flip() +
  theme(axis.title.x = element_blank(),
    strip.background=element_rect(fill="white"),
    strip.text.x = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7))

ggsave(
  p3,
  filename = "fig4b.svg",
  device = "svg",
  path = "figures/plots/Figure4/",
  width = 120, 
  height = 45, 
  units = "mm")


# Figure 4c
# ==============
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
pseq_sorted <- readRDS(file = "figures/data/f3_data.rds")

#D0
meta <- data.frame(sample_data(pseq_sorted))
d0 <- meta[grep("D0", meta$sample_alias), ]
ps0 <- phyloseq::subset_samples(pseq_sorted,(sample_names(pseq_sorted) %in% rownames(d0)))
#D7
d7 <- meta[grep("D7", meta$sample_alias), ]
ps7 <-  phyloseq::subset_samples(pseq_sorted,(sample_names(pseq_sorted) %in% rownames(d7)))
#D30
d30 <- meta[grep("D30", meta$sample_alias), ]
ps30 <-  phyloseq::subset_samples(pseq_sorted,(sample_names(pseq_sorted) %in% rownames(d30)))

# prepare mat for plotting
mat1 <- t(as.matrix(as.data.frame(ps0@otu_table)))
mat2 <- t(as.matrix(as.data.frame(ps7@otu_table)))
mat3 <- t(as.matrix(as.data.frame(ps30@otu_table)))
common_min = min(c(mat1, mat2, mat3))
common_max = max(c(mat1, mat2, mat3))
col_fun = circlize::colorRamp2(c(common_min,
                                 ((common_max+common_min)/2),
                                 common_max),
                               c("cornflowerblue","white", "brown3"))
# prepare complex heatmap
# D0
ann <- data.frame(ps0@sam_data$status, ps0@sam_data$p_cluster,
                  ps0@sam_data$CST)

colnames(ann) <- c('Status', 'Cluster', 'CST')
colours <- list(
  'Status' = c('Cured' = "#440154FF", 'Failed' = "#FDE725FF"),
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
                              'Status' =list(nrow = 1,
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

h0 <- Heatmap(mat1, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_gp = gpar(fontsize = lts),
                                          labels_gp = gpar(fontsize = wls),
                                          title_position = "leftcenter-rot"),
              column_title = "Pre-treatment",
              column_title_gp = grid::gpar(fontsize = ts),
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
  'Status' = c('Cured' = "#440154FF", 'Failed' = "#FDE725FF"),
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
                              'Status' =list(nrow = 1,
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

h7 <- Heatmap(mat2, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_gp = gpar(fontsize = lts),
                                          labels_gp = gpar(fontsize = wls),
                                          title_position = "leftcenter-rot"),
              column_title = "After one week",
              column_title_gp = grid::gpar(fontsize = ts),
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
  'Status' = c('Cured' = "#440154FF", 'Failed' = "#FDE725FF"),
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
                              'Status' =list(nrow = 1,
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

h30 <- Heatmap(mat3, name = "CLR Species Abundances",
              heatmap_legend_param = list(legend_height = unit(6, "cm"),
                                          title_gp = gpar(fontsize = lts),
                                          labels_gp = gpar(fontsize = wls),
                                          title_position = "leftcenter-rot"),
              column_title = "After one month",
              column_title_gp = grid::gpar(fontsize = ts),
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = rns),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

hlist <- h0 + h7 + h30

svg(filename = "figures/plots/Figure4/fig4c.svg",
    width = 7.08661, 
    height = 3.93701)

draw(object = hlist, heatmap_legend_side = "right", 
     annotation_legend_side = "bottom")
dev.off()


# Figure 4d
# ============

sp <- readRDS("04_treatment/res/BV_Microbiome_treatment/summary_performances.rds")
sp$data[sp$data == "valencias"] <- "VALENCIA"
sp$data[sp$data == "all"] <- "VALENCIA-VIBES"
sp$data[sp$data== "clusters"] <- "VIBES"
sp$data <- factor(sp$data, levels = c("VIBES", "VALENCIA", "VALENCIA-VIBES"))

p4 <- sp %>%
  subset(algorithm == "rf") %>%
  filter((data %in% c("VALENCIA-VIBES", "VALENCIA"))) %>%
  group_by(data, variable) %>%
  summarize(mean_value = mean(value)) %>%
  ggscatter(
    x = "mean_value",
    y = "data",
    color = "data",
    shape = "data",
    size = 2
  ) + 
  theme(legend.position = "none") +
  facet_wrap("variable", scales = "free_x") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 7),
        axis.text.y = element_text(color = "black",
                                   size = 7),
        strip.text.x = element_text(size = 7),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white"),
        legend.position = "none"
  )+ scale_color_manual(values = viridis(4)[2:3])
ggsave(
  p4,
  filename = "fig4d.svg",
  device = "svg",
  path = "figures/plots/Figure4/",
  width = 135, 
  height = 40, 
  units = "mm")


