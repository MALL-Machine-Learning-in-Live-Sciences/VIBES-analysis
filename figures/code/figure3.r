#### Figure 3 ####
# ForestPlot
require("forestplot")
require(dplyr)
ext_val = readRDS("03_machine_learning/res/performances_extval.rds")
names(ext_val)[2] <- "mean"
names(ext_val)[3] <- "lower"
names(ext_val)[4] <- "upper"
ext_val <- ext_val %>% 
  mutate_if(is.numeric, round, digits = 3)

base_data <- tibble::tibble(ext_val)


# meter linea en 0.5, vamciar el eje X
# Meter el nuymero de samples (n)
# m3 ventilo todo menos acc en summary
# MEter ravel

base_data |>
  forestplot(labeltext = c(cohort, brier, kappa, bacc),
             vertice = TRUE,align = "c",
             xlab = "Accuracy",
             xlog=TRUE,
             xticks = c(-0.46,-0.10, -0.05, 0)) |>
  fp_add_lines(h_2 = gpar(lty = 2), 
               h_6 = gpar(lwd = 1, columns = 1:4, col = "#000044")) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_header(cohort = c("Cohort"),
                brier = c("Brier"),
                kappa = c("Kappa"),
                bacc = c("Balanced Acc.")) |>
  fp_append_row(cohort = "Summary",
                brier = 0.120,
                kappa =  0.883,
                bacc = 0.933,
                mean  = 0.923,
                lower = 0.8885,
                upper = 0.9485,
                is.summary = TRUE) |> 
  fp_set_zebra_style("#EFEFEF")
  
# Feature Importance
# 1.Extract betas and reduce abs for all class
xx <- readRDS("03_machine_learning/model/pruned_glmnet.rds")
species <- unique(c(rownames(xx$nbeta$N),
                    rownames(xx$nbeta$IDN),
                    rownames(xx$nbeta$IDD),
                    rownames(xx$nbeta$D)))
base <- rbind(abs(xx$nbeta$N),
             abs(xx$nbeta$IDN),
             abs(xx$nbeta$IDD),
             abs(xx$nbeta$D))
fi <- rowsum(as.matrix(base), row.names(base))
d_fi <- data.frame(species = rownames(fi),
                   importance = fi[,1]) [ -1,]
d_fi <- d_fi[order(d_fi$importance,decreasing = TRUE),]

# 2.Plot importances
require(ggplot2)
require(viridis)
plotFIC = ggplot(data = d_fi, aes(x=reorder(species, importance), y = importance))+
  geom_segment( aes(xend=species, yend=0,), color = viridis(3)[2]) +
  geom_point( size=2, color = viridis(1)) +
  coord_flip() +
  theme_light(base_size = 16)+
  theme( axis.text=element_text(size=9),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))
plotFIC

# Facet !!