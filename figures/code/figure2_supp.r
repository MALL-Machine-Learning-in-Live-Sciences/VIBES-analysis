# Plotting CV benchmark
# =====
require(ggplot2)
require(ggpubr)
toplot1 <- readRDS("03_machine_learning/res/summary_performances.rds")
toplot2 <- readRDS("03_machine_learning/res/summary_performances_extval.rds")
# Rename algorithms 
toplot2[toplot2 == "classif.glmnet.tuned"] <- "glmnet"
toplot2[toplot2 == "classif.ranger.tuned"] <- "rf"
toplot2[toplot2 == "classif.svm.tuned"] <- "svm"
toplot2[toplot2 == "classif.xgboost.tuned"] <- "xgboost"
# Rename variable
toplot2$variable <- as.character(toplot2$variable)
toplot2[toplot2 == "Brier Error"] <- "Brier"
toplot2[toplot2 == "Balanced Accuracy"] <- "BAccuracy"
toplot2$variable <- as.factor(toplot2$variable)
toplot <- rbind(toplot1, toplot2)
# Rename Studies
toplot[toplot == "Ravel"] <- "Train"
toplot[toplot == "Sriniv"] <- "Validation 1"
toplot[toplot == "PRJNA2085"] <- "Validation 2"
toplot[toplot == "PRJNA7977"] <- "Validation 3"
toplot[toplot == "PRJNA3020"] <- "Validation 4"
toplot$data <- factor(toplot$data,
                       levels = c("Validation 4","Validation 3","Validation 2",
                                  "Validation 1","Train"))
toplot$value <- round(toplot$value, digits = 2)
ggscatter(
    data = toplot,
    palette = c("#440154FF", "#31688EFF", "#35B779FF", "#EB8055FF"),
    x = "value",
    y = "data",
    color = "algorithm",
    shape = "algorithm",
    size = 5
    ) +
  facet_grid(~ variable + transformation, scales = "free_x") + 
  theme(trip.background = element_rect(colour = "black",
                                       fill = "white",
                                       linewidth = 1.5,
                                       linetype="solid"))+
  theme_light(base_size = 12) +
  scale_fill_continuous(guide = guide_legend()) +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        strip.text.x = element_text(size=12, 
                                    face = "bold",
                                    color=viridis(1)),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
