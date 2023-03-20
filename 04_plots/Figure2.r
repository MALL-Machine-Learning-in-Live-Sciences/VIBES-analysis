#Figure 2
require(plotly)
library(mlr)
library(ggplot2)
library(viridis)
### The index plots are extracted from the ClusterInvestigation.R script.
#Genus
#Cluster Visualization
a = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Cluster/3C/FS/Ravel_3C_3.rds")
a$cluster = as.factor(a$cluster)
levels(a$cluster) <- c("C3", "C2", "C1")
a$cluster = as.character(a$cluster)

clusterG = plot_ly(x = a$Prevotella, 
                   y = a$Dialister, 
                   z = a$Sneathia,
                   type = "scatter3d", 
                   mode = "markers",colors =  c("#21908CFF","#440154FF","#FDE725FF" ),
                   color = as.factor(a$cluster)) %>%
  layout(legend = list(size = 35,orientation = "h", x = 0.3, y = 0),
         scene = list(xaxis = list(title = "Prevotella", tickfont = list(size = 12)),
                      yaxis = list(title = "Dialister", tickfont = list(size = 12)),
                      zaxis = list(title = "Sneathia", tickfont = list(size = 12))))
clusterG

#Species
#Cluster Visualization
a = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/FS/Ravel_3C_3.rds")
names(a)[2] <- 'Atopobium.vaginae'
clusterS = plot_ly(x = a$Lactobacillus.crispatus, 
                   y = a$Atopobium.vaginae, 
                   z = a$Sneathia.sanguinegens,
                   type = "scatter3d", 
                   mode = "markers",colors =  c("#21908CFF","#440154FF","#FDE725FF" ),
                   color = as.factor(a$cluster)) %>%
  layout(legend = list(size = 35,orientation = "h", x = 0.3, y = 0),
         scene = list(xaxis = list(title = "Lactobacillus.crispatus",tickfont = list(size = 12)),
                      yaxis = list(title = "Atopobium.vaginae", tickfont = list(size = 12)),
                      zaxis = list(title = "Sneathia.sanguinegens", tickfont = list(size = 12))))
clusterS


