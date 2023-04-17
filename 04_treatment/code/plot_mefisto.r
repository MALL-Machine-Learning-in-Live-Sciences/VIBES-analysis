# Explore the model
# =========
require(pheatmap)
library(cowplot)
library(magrittr)
require(MOFA2)
require(ggplot2)

# Load model
sm <- MOFA2::load_model("~/git/vaginosis-jlb/MEFISTO/res/PRJNA797778.rds")

# Variance that a factor explains in each sample
p <- plot_variance_explained(
  sm, 
  plot_total = T,
  x = "group", 
  split_by = "view")

p[[1]] + theme(
  axis.text.x = element_blank()) + 
  xlab("sample")

p[[2]] + theme(
  axis.text.x = element_blank()) + 
  xlab("sample")


# Factors versus time
plot_factors_vs_cov(
  sm, 
  color_by = "time") +
  stat_summary(aes(col = color_by), geom="line", fun = "mean")

# Plot weigths
plot_weights(sm, factors = 1)
plot_top_weights(sm, factors = 1)

plot_data_vs_cov(
  sm,
  factor = 1,
  features = 2,
  color_by = "time",
  dot_size = 1
)

kk <- MOFA2::cluster_samples(sm, k = 3)
MOFA2::plot_group_kernel(sm, factors = 1)

factors <- MOFA2::get_factors(sm)
factors$AYAC03

plot_data_heatmap(
  sm,
  factor = 2)

# Scatterplot
gg_scatter <- plot_grid(
  plot_factors(sm, color_by = "delivery") +
    theme(legend.position = "top"),
  plot_factors(sm, color_by = "status") +
    theme(legend.position = "top"),
  plot_factors(sm, color_by = "time") +
    theme(legend.position = "top"),
  nrow = 1, align = "h", axis = "tb")

gg_scatter


plot_factor_cor(sm)
get_scales(sm)
plot_factors_vs_cov(sm, color_by = "time")
plot_weights(sm, factors = 4, view = 1)
plot_top_weights(sm,
                 factors = 3,
                 view = 2)
plot_data_vs_cov(sm,
                 factor=3,
                 features = 2,
                 color_by = "time",
                 dot_size = 1)
