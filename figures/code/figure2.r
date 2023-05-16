#### Figure 2 ####
# ForestPlot
require("forestplot")
require(dplyr)
base_data = readRDS(file = "figures/data/f2_forest_data.rds")

base_data |>
  forestplot(labeltext = c(cohort, n, brier, kappa, bacc),
             vertices = TRUE,
             xlab = "Accuracy",
             clip = c(.5, 0.7, 1),
             xlog = TRUE,
             grid = structure(c(0.75), 
                                          gp = gpar(lty = 2, col = "#CCCCFF")),
             xticks = c(log(0.75), log(0.80), log(0.85),
                        log(0.9), log(0.95), log(1))
             ) |>
  fp_add_lines(h_1 = gpar(lwd = 2, columns = 1:6, col = "#000044"),
               h_2 = gpar(lwd = 1, columns = 1:6, col = "#440154FF"),
               h_6 = gpar(lwd = 1, columns = 1:6, col = "#440154FF"),
               h_7 = gpar(lwd = 2, columns = 1:6, col = "#000044")) |>
  fp_set_style(box = "#440154FF",
               line = "#453781FF",
               summary = "#440154FF",
               align = "lcccc",
               txt_gp = fpTxtGp(ticks = gpar(cex = 1),
                                xlab  = gpar(cex = 1,
                                             fontface = "bold"))) |> 
  fp_add_header(cohort = c("Cohort"),
                n = c("N"),
                brier = c("Brier"),
                kappa = c("Kappa"),
                bacc = c("Bal. Acc.")) |>
  fp_append_row(cohort = "Summary",
                mean  = 0.923,
                lower = 0.8885,
                upper = 0.9485,
                is.summary = TRUE) |> 
  fp_set_zebra_style("#EFEFEF")

svg(
  filename = "figures/plots/Figure2/fig2.svg",
  width = 7, 
  height = 2.5,
  pointsize = 8
)
base_data |>
       forestplot(labeltext = c(cohort, n, brier, kappa, bacc),
                  vertices = TRUE,
                  xlab = "Accuracy",
                  clip = c(.5, 0.7, 1),
                  xlog = TRUE,
                  grid = structure(c(0.75), 
                                   gp = gpar(lty = 2, col = "#CCCCFF")),
                  xticks = c(log(0.75), log(0.80), log(0.85),
                             log(0.9), log(0.95), log(1))
       ) |>
       fp_add_lines(h_1 = gpar(lwd = 2, columns = 1:6, col = "#000044"),
                    h_2 = gpar(lwd = 1, columns = 1:6, col = "#440154FF"),
                    h_6 = gpar(lwd = 1, columns = 1:6, col = "#440154FF"),
                    h_7 = gpar(lwd = 2, columns = 1:6, col = "#000044")) |>
       fp_set_style(box = "#440154FF",
                    line = "#453781FF",
                    summary = "#440154FF",
                    align = "lcccc",
                    txt_gp = fpTxtGp(ticks = gpar(cex = 1),
                                     xlab  = gpar(cex = 1,
                                                  fontface = "bold"))) |> 
       fp_add_header(cohort = c("Cohort"),
                     n = c("N"),
                     brier = c("Brier"),
                     kappa = c("Kappa"),
                     bacc = c("Bal. Acc.")) |>
       fp_append_row(cohort = "Summary",
                     mean  = 0.923,
                     lower = 0.8885,
                     upper = 0.9485,
                     is.summary = TRUE) |> 
       fp_set_zebra_style("#EFEFEF")
dev.off()

# Feature Importance
# Declare base
xx <- readRDS("03_machine_learning/model/pruned_glmnet.rds")
species <- unique(c(rownames(xx$nbeta$N),
                    rownames(xx$nbeta$IDN),
                    rownames(xx$nbeta$IDD),
                    rownames(xx$nbeta$D)))
base <- matrix(rep(x = as.numeric(0), 21))
row.names(base) <- species
colnames(base) <- "importance"

# For each class
fi_n <- data.frame(rowsum(x = rbind(base, as.matrix(abs(xx$nbeta$N))),
                          group = row.names(rbind(base, as.matrix(abs(xx$nbeta$N))))))
fi_idn <- data.frame(rowsum(x = rbind(base, as.matrix(abs(xx$nbeta$IDN))),
                            group = row.names(rbind(base, as.matrix(abs(xx$nbeta$IDN))))))
fi_idd <- data.frame(rowsum(x = rbind(base, as.matrix(abs(xx$nbeta$IDD))),
                            group = row.names(rbind(base, as.matrix(abs(xx$nbeta$IDD))))))
fi_d <- data.frame(rowsum(x = rbind(base, as.matrix(abs(xx$nbeta$D))),
                          group = row.names(rbind(base, as.matrix(abs(xx$nbeta$D))))))

d_n <- data.frame(species = rownames(fi_n),
                  importance = fi_n[,1]) [ -1,]
d_idn <- data.frame(species = rownames(fi_idn),
                    importance = fi_idn[,1]) [ -1,]
d_idd <- data.frame(species = rownames(fi_idd),
                    importance = fi_idd[,1]) [ -1,]
d_d <- data.frame(species = rownames(fi_d),
                  importance = fi_d[,1]) [ -1,]

# 2.Plot importances
require(ggplot2)
require(viridis)
plot_d_n = ggplot(data = d_n, aes(x=reorder(species, desc(species)), y = importance))+
  geom_point(size=2)+
  geom_segment( aes(xend=species, yend=0), size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        strip.text.x = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) + coord_flip() + facet_grid(. ~ "VCS-I")
  



plot_d_idn = ggplot( data = d_idn,  aes(x=reorder(species, desc(species)), y = importance))+
  geom_point(size=2)+
  geom_segment( aes(xend=species, yend=0), size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        strip.text.x = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) + coord_flip() + facet_grid(. ~ "VCS-II")


plot_d_idd = ggplot( data = d_idd, aes(x=reorder(species, desc(species)), y = importance))+
  geom_point(size=2)+
  geom_segment( aes(xend=species, yend=0), size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        strip.text.x = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) + coord_flip() + facet_grid(. ~ "VCS-III")


plot_d_d = ggplot(data = d_d, aes(x=reorder(species, desc(species)), y = importance))+
  geom_point(size=2)+
  geom_segment( aes(xend=species, yend=0), size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        strip.text.x = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) + coord_flip() + facet_grid(. ~ "VCS-IV")


require(ggpubr)
gp <- ggarrange(plotlist = list(plot_d_n,
                          plot_d_idn,
                          plot_d_idd,
                          plot_d_d),
          label.x = "importance",
          ncol = 4, nrow = 1)

ggsave(
  gp,
  filename = "fig2b.svg",
  device = "svg",
  path = "figures/plots/Figure2/",
  width = 148, 
  height = 70, 
  units = "mm")

#Plot again for extract names
plot_d_n = ggplot(data = d_n, aes(x=reorder(species, desc(species)), y = importance))+
  geom_point(size=2)+
  geom_segment( aes(xend=species, yend=0), size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8,
                                   hjust = 1,
                                   color = "black"),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        strip.text.x = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(color = "black"),
        legend.key = element_rect(fill = "transparent"), 
        strip.text = element_text(size = rel(1.2)),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.grid.major.y = element_blank(),
        strip.background=element_rect(fill="white")
  ) + coord_flip() + facet_grid(. ~ "VCS-I")

gp_names = ggarrange(plotlist = list(plot_d_n,
                          plot_d_idn,
                          plot_d_idd,
                          plot_d_d),
          ncol = 4, nrow = 1)

ggsave(
  gp_names,
  filename = "fig2b_names.svg",
  device = "svg",
  path = "figures/plots/Figure2/",
  width = 180, 
  height = 70, 
  units = "mm")
