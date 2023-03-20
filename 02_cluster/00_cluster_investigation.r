paquetes <- c("tidyverse", "cluster", "factoextra", "NbClust",
"tidyr", "viridis", "plotly")
lapply(paquetes, require, character.only = TRUE)
setwd("~/git/BV_Microbiome/")
rank <-  "Species" # Genus or "Species"
nfeat <- "24" # 24 and 14 for Species or 37 and 21 for Genus

Ravel <- readRDS(paste0("extdata/", rank, "Intersect/Ravel_",
rank, "_pseq_", nfeat, ".rds"))
Sriniv <-readRDS(paste0("extdata/", rank, "Intersect/Sriniv_",
rank, "_pseq_", nfeat, ".rds"))

# Select Target
select.target <- function(pseq, tget) {
  require(dplyr)
  require(phyloseq)
  df <- sample_data(pseq)
  colnames(df)[colnames(df) == tget] <- "target"
  sample_data(pseq) <- df
  return(pseq)
}
norm.dataset <- function(pseq) {
  # Change columns by rows too, interested in maintain fts in columns
  require(phyloseq)
  otu <- data.frame(otu_table(pseq))
  # Normalize that variables
  otu <- apply(X = otu, FUN = function(x) log2(x + 1), MARGIN = 2)
  otu_table(pseq) <- otu_table(otu, taxa_are_rows = FALSE)
  return(pseq)
}

#3C
Ravel_3C = select.target(pseq = Ravel, tget ="Nugent_score_category")
Sriniv_3C = select.target(pseq = Sriniv, tget = "nugent")

#Norm Datasets
Ravel_3C = norm.dataset(Ravel_3C)
Sriniv_3C = norm.dataset(Sriniv_3C)

#Cluster estimation
m.distancia <- get_dist(otu_table(Ravel_3C), method = "euclidean") 
fviz_dist(m.distancia,show_labels = FALSE, gradient = list(low = "blue", mid = "white", high = "red"))
#Elbow, silhouette o gap_stat  method
fviz_nbclust(otu_table(Ravel_3C), kmeans, method = "wss")
fviz_nbclust(otu_table(Ravel_3C), kmeans, method = "silhouette")
fviz_nbclust(otu_table(Ravel_3C), kmeans, method = "gap_stat")

#All inedx:
set.seed(1312)
resnumclust2<-NbClust(otu_table(Ravel_3C), distance = "euclidean", min.nc=2, max.nc=6, method = "kmeans", index = "alllong")
fviz_nbclust(resnumclust2)
dfresnumclust2 = data.frame("Clusters" = c(2,3,4,6),"Freq" = c(9,15,2,1))
dfresnumclust2$Freq = as.factor(dfresnumclust2$Freq)
dfresnumclust2$Clusters = as.factor(dfresnumclust2$Clusters)

#Plot Index
clIG = ggplot(data =dfresnumclust2, aes(x = Clusters, y = Freq, fill = Freq))+
  geom_bar(position = 'dodge', stat='identity', width=.5)+ theme_light(base_size = 13)+
  ylab("Freq. among all Indices")+ xlab("Optimal Number of Clusters")
require(ggpubr)
clIG = change_palette(clIG, palette= viridis(10)[1:4])
clIG = clIG + theme(legend.title = element_blank(),legend.position = "none",plot.title = element_text(hjust = 0.5))
clIG
#Clusters
k2 <- kmeans(otu_table(Ravel_3C), centers = 3, nstart = 25,iter.max = 100)
identical(names(k2$cluster), rownames(sam_data(Ravel_3C)))
sample_data(Ravel_3C)$cluster <- k2$cluster
#Plot Clusters
clusterG = plot_ly(x = data.frame(otu_table(Ravel_3C))$OTU_20, 
        y = data.frame(otu_table(Ravel_3C))$OTU_6, 
        z = data.frame(otu_table(Ravel_3C))$OTU_13,
        type = "scatter3d", 
        mode = "markers",colors =  c("#FDE725FF","#440154FF","#21908CFF" ),
        color = as.factor(sample_data(Ravel_3C)$cluster)) %>%
  layout(legend = list(x = 0.8, y = 0.5,size = 20),
         scene = list(xaxis = list(title = "Prevotella",tickfont = list(size = 15)),
                      yaxis = list(title = "Dialister", tickfont = list(size = 15)),
                      zaxis = list(title = "Sneathia", tickfont = list(size = 15))))
clusterG
xaxis = list(tickfont = list(size = 15))
