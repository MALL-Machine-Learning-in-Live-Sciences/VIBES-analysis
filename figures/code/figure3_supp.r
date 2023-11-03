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
