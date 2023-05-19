#### Figure 2 ####
# ForestPlot
require("forestplot")
require(dplyr)
base_data = readRDS(file = "figures/data/f2_forest_data.rds")

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

svg(
  filename = "figures/plots/Figure2/fig2.svg",
  width = 7, 
  height = 2.5,
  pointsize = 8
)
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
dev.off()


