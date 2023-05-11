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
  split_by = "view")

p[[1]] + theme(
  axis.text.x = element_blank()) + 
  xlab("sample")

plot_factor_cor(sm)

# ML results 
rs = readRDS("04_treatment/res/summary_performances.rds")
ggscatter(
  data = rs,
  x = "value",
  y = "data",
  color = "algorithm",
  shape = "algorithm",
  size = 5
) +
  facet_wrap(~transformation + variable, scales = "free_x") +
  theme_bw() + theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     legend.position = "bottom")
