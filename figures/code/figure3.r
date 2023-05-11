# MEFISTO
# Explore the model
# =========
require(pheatmap)
library(cowplot)
library(magrittr)
require(MOFA2)
require(ggplot2)
library(RColorBrewer)
# Load model
sm <- MOFA2::load_model("04_treatment/res/PRJNA302078.rds")

# Feature data
feat <- readRDS("04_treatment/data/PRJNA3020_featdata.rds")
identical(feat$feature, features_metadata(sm)$feature)
feat$feature <- feat$Species
features_metadata(sm) <- feat

# Variance that a factor explains in each sample
p <- plot_variance_explained(
  sm, 
  plot_total = T,
  x = "group", 
  split_by = "view",)
p[[1]] +
  theme(
  axis.text.x = element_blank()) + scale_fill_gradient(low="white", high=viridis(1)) +
  xlab("samples")

plot_factor_cor(sm)

# Factors versus time 517 x 528
# F1
plot_factors_vs_cov(factors = 1,
  sm,
  dot_size = 6,
  shape_by = "status",
  color_by = "status") + theme_light(base_size = 12) +
  
  theme(legend.position =c(0.6, 0.2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        strip.text.x = element_text(size = 18, 
                                    face = "bold",
                                    color = viridis(1)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values=c ("green", "red")) +
  scale_y_continuous(position = "right") +
  stat_summary(aes(col = color_by),
               geom = "line",
               size = 1.1, linetype = "dotdash",
               fun = "mean",
               show.legend = FALSE) +
  scale_colour_manual(values = c("green", "red"))

# F2
plot_factors_vs_cov(factors = 2,
                    sm,
                    dot_size = 6,
                    shape_by = "status",
                    color_by = "status") + theme_light(base_size = 12) +
  
  theme(legend.position =c(0.6, 0.2),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        strip.text.x = element_text(size = 18, 
                                    face = "bold",
                                    color = viridis(1)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values=c ("green", "red")) +
  scale_y_continuous(position = "right") +
  stat_summary(aes(col = color_by),
               geom = "line",
               size = 1.1, linetype = "dotdash",
               fun = "mean",
               show.legend = FALSE) +
  scale_colour_manual(values = c("green", "red"))


# Plot weigths
rownames(sm@expectations[["W"]][["microbiome"]]) <- feat$Species
# Edit for colors and size
trace(plot_top_weights,
      edit=TRUE)
# F1
plot_top_weights(sm,
                 factors = 1,
)  + theme_light(base_size = 12) +
  scale_colour_gradient(viridis(1)) +
  theme(legend.position ="bottom",
        legend.text = element_text(size = 18),
        strip.text.x = element_text(size = 18, 
                                    face = "bold",
                                    color = viridis(1)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

#F2
plot_top_weights(sm,
                 factors = 2,
)  + theme_light(base_size = 12) +
  scale_colour_gradient(viridis(1)) +
  theme(legend.position ="bottom",
        legend.text = element_text(size = 18),
        strip.text.x = element_text(size = 18, 
                                    face = "bold",
                                    color = viridis(1)),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())



#HM 1361 x 725
model <- readRDS("03_machine_learning/model/pruned_glmnet.rds")
pruned_predict <- function(pruned_glmnet, newdata){
  require(glmnet)
  # 1.Pre-processing of data for betas multiplication
  dd = dim(newdata)
  if (inherits(newdata, "sparseMatrix"))
    newdata = as(newdata, "dMatrix")
  npred = dd[[1]]
  dn = list(names(pruned_glmnet$nbeta), dimnames(pruned_glmnet$nbeta[[1]])[[2]], dimnames(newdata)[[1]])
  dp = array(0, c(pruned_glmnet$nclass, pruned_glmnet$nlambda, npred), dimnames = dn)
  # 2.Multiplication by betas
  for (i in seq(pruned_glmnet$nclass)) {
    rn <- rownames(pruned_glmnet$nbeta[[i]])
    rn <- rn[-1]
    fitk = cbind2(1, newdata[,rn]) %*% (pruned_glmnet$nbeta[[i]])
    dp[i, , ] = dp[i, , ] + t(as.matrix(fitk))
  }
  # 3.Results are extracted
  pp = exp(dp)
  psum = apply(pp, c(2, 3), sum)
  # Response
  response = aperm(pp/rep(psum, rep(pruned_glmnet$nclass, pruned_glmnet$nlambda * npred)), c(3,1, 2))
  response <- as.data.frame(response)
  colnames(response) <- substr(colnames(response), 1, nchar(colnames(response))-2)
  cp = aperm(dp, c(3, 1, 2))
  # Class
  class <- apply(cp, 3, glmnet:::glmnet_softmax)
  #Bind both
  response$p_cluster <- class[,]
  return(response)
}
d <- readRDS("03_machine_learning/data/PRJNA3020_clr.rds")
d <- subset(d,
            select = -c(Mageeibacillus_indolicus, Peptoniphilus_lacrimalis))
p <- pruned_predict(pruned_glmnet = model,
                    newdata = as.matrix(d[,1:20]))
identical(rownames(p), rownames(d))
dp <- cbind(d, p)
packgs <- c("phyloseq", "ComplexHeatmap", "microViz")
lapply(packgs, require, character.only = TRUE)
setwd("~/git/BV_Microbiome/")
cl <- readRDS("02_cluster/data/C4/PRJNA3020_Cluster_Species_22_pseq.rds")
identical(rownames(dp), rownames(cl@sam_data))
sample_data(cl) <- cbind(cl@sam_data, dp[,22:26])

vl <- read.csv2(file  = "01_get_valencias/res/PRJNA3020.csv", sep = ",", header = TRUE)
# check same order
vl <- data.frame(vl[,-1], row.names = vl[,1])
identical(rownames(cl@sam_data), rownames(vl))
# add valencia 
cl@sam_data$CST <- vl$CST
# order pseq by probbilities of belong to own class
df <- cl@sam_data

dfn <- df[df$p_cluster == "N",]
dfn <- dfn[order(dfn$N, decreasing = TRUE),]
dfidn <- df[df$p_cluster == "IDN",]
dfidn <- dfidn[order(dfidn$IDN, decreasing = TRUE),]
dfidd <- df[df$p_cluster == "IDD",]
dfidd <- dfidd[order(dfidd$IDD, decreasing = TRUE),]
dfd <- df[df$p_cluster == "D",]
dfd <- dfd[order(dfd$D, decreasing = TRUE),]
pseq_sorted <- ps_reorder(ps = cl,
                          sample_order = rownames(rbind(dfd,
                                                        dfidd,
                                                        dfidn,
                                                        dfn)))

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
  'Cluster' = c('IDD' = "#A6CEE3" ,'N' = '#33A02C',
                'D' = '#1F78B4', 'IDN' = '#B2DF8A'),
  'CST' = c('III' = 'yellowgreen', 'IV-A' = 'darkgray',
            'IV-B' = 'seagreen', 'IV-C' = 'mediumvioletred')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 11),
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
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 12),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

#D7
ann <- data.frame(ps7@sam_data$status, ps7@sam_data$p_cluster,
                  ps7@sam_data$CST)

colnames(ann) <- c('Status', 'Cluster', 'CST')
colours <- list(
  'Status' = c('Cured' = 'green', 'Failed' = 'red'),
  'Cluster' = c('IDD' = "#A6CEE3" ,'N' = '#33A02C',
                'D' = '#1F78B4', 'IDN' = '#B2DF8A'),
  'CST' = c('I'= 'cornflowerblue', 'III' = 'yellowgreen','IV-B' = 'seagreen',
            'IV-C' = 'mediumvioletred', 'V' = 'tan2')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 11),
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
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 12),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

#D30
ann <- data.frame(ps30@sam_data$status, ps30@sam_data$p_cluster,
                  ps30@sam_data$CST)

colnames(ann) <- c('Status', 'Cluster', 'CST')
colours <- list(
  'Status' = c('Cured' = 'green', 'Failed' = 'red'),
  'Cluster' = c('IDD' = "#A6CEE3" ,'N' = '#33A02C',
                'D' = '#1F78B4', 'IDN' = '#B2DF8A'),
  'CST' = c('I'= 'cornflowerblue', 'III' = 'yellowgreen','IV-A' = 'darkgray',
            'IV-B' = 'seagreen','IV-C' = 'mediumvioletred')
)
colAnn <- HeatmapAnnotation(df = ann,
                            which = 'col',
                            col = colours,
                            annotation_name_gp = gpar(fontsize = 11),
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
              col = col_fun,
              top_annotation = colAnn,
              row_names_side = "left",
              row_names_gp = gpar(fontsize = 12),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              show_column_names = FALSE)

hlist <- h0 + h7 + h30
h <- draw(object = hlist, heatmap_legend_side = "right", 
          annotation_legend_side = "bottom")

# ML results 
rs = readRDS("04_treatment/res/summary_performances.rds")
ggscatter(
  data = rs,
  palette = c("#440154FF", "#31688EFF", "#35B779FF", "#EB8055FF"),
  x = "value",
  y = "transformation",
  color = "algorithm",
  shape = "algorithm",
  size = 6
) +
  facet_wrap(~ variable, scales = "free_x") +
  theme_bw() + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.position = "bottom") +
theme(trip.background = element_rect(colour = "black",
                                     fill = "white",
                                     linewidth = 1.5,
                                     linetype="solid"))+
  theme_light(base_size = 12) +
  scale_fill_continuous(guide = guide_legend()) +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(size = 14, 
                                    face = "bold",
                                    color = viridis(1)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
