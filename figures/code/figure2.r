#### Figure 3 ####
# ForestPlot
require("forestplot")
require(dplyr)
ext_val = readRDS("03_machine_learning/res/performances_extval.rds")
ext_val$cohort[ext_val$cohort == 'Ravel'] <- 'Train'
ext_val$cohort[ext_val$cohort == 'Sriniv'] <- 'Validation 1'
ext_val$cohort[ext_val$cohort == 'PRJNA2085'] <- 'Validation 2'
ext_val$cohort[ext_val$cohort == 'PRJNA7977'] <- 'Validation 3'
ext_val$cohort[ext_val$cohort == 'PRJNA3020'] <- 'Validation 4'
ext_val <- ext_val[c(4,5,1,3,2),]
ext_val$n <-c(394, 220, 1657, 220, 201)
names(ext_val)[2] <- "mean"
names(ext_val)[3] <- "lower"
names(ext_val)[4] <- "upper"
ext_val <- ext_val %>% 
  mutate_if(is.numeric, round, digits = 3)

base_data <- tibble::tibble(ext_val)
base_data <- base_data[-1,]

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

#x <- unit(.22, 'npc')
#y <- unit(.15, 'npc')
#grid.text('*Accuracy summary did not take into account the Train performance', x, y, gp = gpar(fontsize=10, font = 3))


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
  geom_segment( aes(xend=species, yend=0,), color = viridis(3)[2]) +
  geom_point( size=2, color = viridis(1)) +
  coord_flip() +
  theme_light(base_size = 16) +
  ggtitle("N cluster") +
  theme( axis.text=element_text(size=15),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         legend.title=element_text(size=10), 
         legend.text=element_text(size=10),
         plot.title = element_text(hjust = 0.5))


plot_d_idn = ggplot( data = d_idn, aes(x=reorder(species, desc(species)), y = importance))+
  geom_segment( aes(xend=species, yend=0,), color = viridis(3)[2]) +
  geom_point( size=2, color = viridis(1)) +
  coord_flip() +
  theme_light(base_size = 16) +
  ggtitle("IDN cluster") +
  theme( axis.text=element_text(size=15),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         legend.title=element_text(size=10), 
         legend.text=element_text(size=10),
         plot.title = element_text(hjust = 0.5))


plot_d_idd = ggplot( data = d_idd, aes(x=reorder(species, desc(species)), y = importance))+
  geom_segment( aes(xend=species, yend=0,), color = viridis(3)[2]) +
  geom_point( size=2, color = viridis(1)) +
  coord_flip() +
  theme_light(base_size = 16) +
  ggtitle("IDD cluster") +
  theme( axis.text=element_text(size=15),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         legend.title=element_text(size=10), 
         legend.text=element_text(size=10),
         plot.title = element_text(hjust = 0.5))


plot_d_d = ggplot(data = d_d, aes(x=reorder(species, desc(species)), y = importance))+
  geom_segment( aes(xend=species, yend=0,), color = viridis(3)[2]) +
  geom_point( size=2, color = viridis(1)) +
  coord_flip() +
  theme_light(base_size = 16) +
  ggtitle("D cluster") +
  theme( axis.text=element_text(size=15),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x = element_blank(),
         legend.title=element_text(size=10), 
         legend.text=element_text(size=10),
         plot.title = element_text(hjust = 0.5))


require(ggpubr)
ggarrange(plotlist = list(plot_d_n,
                          plot_d_idn,
                          plot_d_idd,
                          plot_d_d),
          label.x = "importance",
          ncol = 4, nrow = 1)

#Plot again for extract names
plot_d_n = ggplot(data = d_n, aes(x=reorder(species, desc(species)), y = importance))+
  geom_segment( aes(xend=species, yend=0,), color = viridis(3)[2]) +
  geom_point( size=2, color = viridis(1)) +
  coord_flip() +
  theme_light(base_size = 16) +
  ggtitle("N cluster") +
  theme( axis.text=element_text(size=15),
         axis.title.y = element_blank(),
         axis.title.x = element_blank(),
         axis.ticks.x = element_blank(),
         legend.title=element_text(size=10), 
         legend.text=element_text(size=10),
         plot.title = element_text(hjust = 0.5))

ggarrange(plotlist = list(plot_d_n,
                          plot_d_idn,
                          plot_d_idd,
                          plot_d_d),
          label.x = "importance",
          ncol = 4, nrow = 1)
