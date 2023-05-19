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
    size = 2
    ) +
  facet_grid(~ variable + transformation,
             scales = "free_x") + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "black",
                                   size = 7),
        axis.text.y = element_text(color = "black",
                                   size = 7),
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
  width = 220, 
  height = 70, 
  units = "mm")
