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
    size = 5
    ) +
  facet_wrap(~transformation + variable, scales = "free_x") +
  theme_bw() + theme(axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position = "none")


# Plotting external validation benchmark
# =====
require(ggplot2)
require(ggpubr)
toplot <- readRDS("03_machine_learning/res/summary_performances_extval.rds")

"classif.glmnet.tuned" = "glmnet"
"classif.ranger.tuned" = "rf"
"classif.svm.tuned" = "svm"
"classif.xgboost.tuned" = "xgboost"

ggscatter(
    data = toplot,
    x = "value",
    y = "data",
    color = "algorithm",
    shape = "algorithm",
    size = 5
    ) +
  facet_grid(~variable + transformation, scales = "free_x",) +
  theme_bw() + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.position = "none")
