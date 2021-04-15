library(devtools)

library(dynamicTreeCut)
datos = train[,1:52]

# Data Scale
datos = scale(datos)
# Dissimilarity matrix
d <- dist(t(datos), method = "euclidean")
# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

#Dinamic tree
dtree = cutreeDynamic(dendro = hc1,cutHeight = 0.90, minClusterSize = 3,method = "tree")
