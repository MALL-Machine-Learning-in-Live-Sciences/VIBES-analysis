
#VenDIAGRAM
set.seed(20190708)
genes <- paste("gene",1:1000,sep="")
x <- list(
  A = sample(genes,300), 
  B = sample(genes,525), 
  C = sample(genes,440),
  D = sample(genes,350)
)
library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
library("ggVennDiagram")
ggVennDiagram(x, label_alpha = 0)
# MIRAR MLR COMO HACE LAS ROC
x <- list(A=1:5,B=2:7,C=3:6,D=4:9)
ggVennDiagram(x)  # 4d venn
ggVennDiagram(x[1:3])  # 3d venn
ggVennDiagram(x[1:2])  # 2d venn

library(VennDiagram)
venn.diagram(x, filename = "venn-4-dimensions.png")
# Further customization
display_venn(
  x,
  category.names = c("Set 1" , "Set 2 " , "Set 3", "Set 4"),
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
  # Numbers
  cex = .9,
  fontface = "italic",
  # Set names
  cat.cex = 1,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.dist = c(0.055, 0.055, 0.1, 0.1)
)
# Curva ROC
a = asROCRPrediction(prediccion)
p = ROCR::performance(a, "tpr", "fpr")
plot(p)



