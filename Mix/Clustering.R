
library(dynamicTreeCut)
# Remove target column
datos = select(train, -(target))

#Data Scale
sdata = scale(x = as.matrix(datos))

# Dissimilarity matrix
d <- dist(t(sdata), method = "euclidean")

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

#Dinamic tree
dtree = cutreeDynamic(dendro = hc1, cutHeight = 24, minClusterSize = 3,method = "tree")
colnames(datos) = as.factor(dtree)
prefixes = unique(sub("\\..*", "", colnames(datos)))
data_clustered = as.data.frame(sapply(prefixes, function(x) rowSums(datos[,startsWith(colnames(datos), x)])))
data_clustered = select(data_clustered, -("0"))
data_clustered = as.data.frame(cbind(data_clustered, target = train$target))
