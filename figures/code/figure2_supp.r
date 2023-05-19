# Plotting CV benchmark
# =====
library(ggplot2)
library(ggpubr)
toplot <- readRDS(file = "figures/data/f2_supp_data.rds")

gp <- ggscatter(
    data = toplot,
    palette = c("#440154FF", "#31688EFF", "#35B779FF", "#EB8055FF"),
    x = "value",
    y = "data",
    color = "algorithm",
    shape = "algorithm",
    size = 5
    ) +
  facet_grid(~ variable + transformation, scales = "free_x") + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 8),
        strip.text.x = element_text(size = 8),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(),
        legend.position = "bottom",
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
  gp,
  filename = "fig2_supp.svg",
  device = "svg",
  path = "figures/plots/Figure2/",
  width = 250, 
  height = 90, 
  units = "mm")



# Feature importance
# ===========
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
  filename = "supp_fig2_varimp.svg",
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
  filename = "supp_fig2_varimp.pdf",
  device = "pdf",
  path = "figures/plots/Figure2/",
  width = 180, 
  height = 70, 
  units = "mm")
