# Plotting CV benchmark
# =====
require(ggplot2)
require(ggpubr)
toplot <- readRDS("03_machine_learning/res/summary_performances.rds")
ggscatter(
    data = toplot,
    x = "value",
    y = "data",
    color = "algorithm",
    shape = "algorithm",
    size = 6
    ) +
  facet_wrap(~variable, scales = "free_x") +
  theme_bw()


# Plotting external validation benchmark
# =====
require(ggplot2)
require(ggpubr)
toplot <- readRDS("03_machine_learning/res/summary_performances_extval.rds")
ggscatter(
    data = toplot,
    x = "value",
    y = "data",
    color = "algorithm",
    shape = "algorithm",
    size = 6
    ) +
  facet_grid(~variable + transformation, scales = "free_x") +
  theme_bw()
