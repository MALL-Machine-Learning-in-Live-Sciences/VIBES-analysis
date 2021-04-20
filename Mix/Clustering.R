
library(dynamicTreeCut)
# Remove target column
datos = select(train, -(target))

#Data Scale
sdata = scale(as.matrix(datos))

# Dissimilarity matrix
d <- dist(t(sdata), method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

#Dinamic tree
dtree = cutreeDynamic(dendro = hc1, cutHeight = 25, minClusterSize = 3,method = "tree")

fviz_cluster(list(data = t(sdata), cluster = dtree))
