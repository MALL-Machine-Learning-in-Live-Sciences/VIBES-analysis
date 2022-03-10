Ravel = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Ravel_Genus_Counts.rds")
Sriniv = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Sriniv_Genus_Counts.rds")


# Select Target
select.target = function(dataset, tget){
  require(dplyr)
  matr = dataset[1:(length(names(dataset)) - 7)] # 7 cause are clinical variables
  target = dataset  %>% select(starts_with(tget))
  dataframe = data.frame(cbind(matr,target))
  names(dataframe)[length(names(dataframe))] <- "target"
  return(dataframe)
}
remove.intermediate = function(dataset){
  df = dataset
  df2 = df[!(df$target == "intermediate"),]
  return(df2)
}
norm.dataset = function(data){
  # Retain only numerics variables
  cols <- sapply(data, is.numeric) 
  
  # Normalize that variables
  data[cols] <- apply(X = data[cols], FUN = function(x) log2(x+1), MARGIN = 2) 
  return(data) 
}

#3C
Ravel_3C = select.target(dataset = Ravel, tget ="Nugent_score_category")
Sriniv_3C = select.target(dataset = Sriniv, tget = "nugent")

#Norm Datasets
Ravel_3C = norm.dataset(Ravel_3C)
Sriniv_3C = norm.dataset(Sriniv_3C)

paquetes <- c("tidyverse","cluster", "factoextra","NbClust","tidyr")
lapply(paquetes, require, character.only = TRUE)

m.distancia <- get_dist(Ravel_3C[1:23], method = "euclidean") 
fviz_dist(m.distancia,show_labels = FALSE, gradient = list(low = "blue", mid = "white", high = "red"))

#CLuster estimation
#Elbow, silhouette o gap_stat  method
fviz_nbclust(Ravel_3C[1:23], kmeans, method = "wss")
fviz_nbclust(Ravel_3C[1:23], kmeans, method = "silhouette")
fviz_nbclust(Ravel_3C[1:23], kmeans, method = "gap_stat")

#All inedx:
#the index to be calculated. This should be one of : "kl", "ch", "hartigan", "ccc", "scott",
#"marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda",
#"pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma",
#"gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw", "all" (all indices except GAP,
#Gamma, Gplus and Tau), "alllong" (all indices with Gap, Gamma, Gplus and Tau included).
set.seed(1312)
resnumclust2<-NbClust(Ravel_3C[1:23], distance = "euclidean", min.nc=2, max.nc=6, method = "kmeans", index = "alllong")
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


#calculamos los tres clústers
k2 <- kmeans(Ravel_3C[1:23], centers = 3, nstart = 25,iter.max = 100)
Ravel_3C$km_cluster <- k2$cluster
a = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Cluster/3C/FS/Ravel_3C_3.rds")
require(viridis)
require(plotly)
#Plot Clusters
clusterG = plot_ly(x = a$Prevotella, 
        y = a$Dialister, 
        z = a$Sneathia,
        type = "scatter3d", 
        mode = "markers",colors =  c("#FDE725FF","#440154FF","#21908CFF" ),
        color = as.factor(a$cluster)) %>%
  layout(legend = list(x = 0.8, y = 0.5,size = 30),
         scene = list(xaxis = list(title = "Prevotella",tickfont = list(size = 15)),
                      yaxis = list(title = "Dialister", tickfont = list(size = 15)),
                      zaxis = list(title = "Sneathia", tickfont = list(size = 15))))
clusterG

xaxis = list(tickfont = list(size = 15))
#Specie
Ravel = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Ravel_Sps_Counts.rds")
Sriniv = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Sriniv_Sps_Counts.rds")


# Select Target
select.target = function(dataset, tget){
  require(dplyr)
  matr = dataset[1:(length(names(dataset)) - 7)] # 7 cause are clinical variables
  target = dataset  %>% select(starts_with(tget))
  dataframe = data.frame(cbind(matr,target))
  names(dataframe)[length(names(dataframe))] <- "target"
  return(dataframe)
}
remove.intermediate = function(dataset){
  df = dataset
  df2 = df[!(df$target == "intermediate"),]
  return(df2)
}
norm.dataset = function(data){
  # Retain only numerics variables
  cols <- sapply(data, is.numeric) 
  
  # Normalize that variables
  data[cols] <- apply(X = data[cols], FUN = function(x) log2(x+1), MARGIN = 2) 
  return(data) 
}

#3C
Ravel_3C = select.target(dataset = Ravel, tget ="Nugent_score_category")
Sriniv_3C = select.target(dataset = Sriniv, tget = "nugent")

#Norm Datasets
Ravel_3C = norm.dataset(Ravel_3C)
Sriniv_3C = norm.dataset(Sriniv_3C)

paquetes <- c("tidyverse","cluster", "factoextra","NbClust","tidyr")
lapply(paquetes, require, character.only = TRUE)

m.distancia <- get_dist(Ravel_3C[1:23], method = "euclidean") #el método aceptado también puede ser: "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" o "kendall"
fviz_dist(m.distancia,show_labels = FALSE, gradient = list(low = "blue", mid = "white", high = "red"))


#Elbow, silhouette o gap_stat  method
fviz_nbclust(Ravel_3C[1:27], kmeans, method = "wss")
fviz_nbclust(Ravel_3C[1:27], kmeans, method = "silhouette")
fviz_nbclust(Ravel_3C[1:27], kmeans, method = "gap_stat")


#the index to be calculated. This should be one of : "kl", "ch", "hartigan", "ccc", "scott",
#"marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", "duda",
#"pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma",
#"gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw", "all" (all indices except GAP,
#Gamma, Gplus and Tau), "alllong" (all indices with Gap, Gamma, Gplus and Tau included).
set.seed(1312)
resnumclust2<-NbClust(Ravel_3C[1:27], distance = "euclidean", min.nc=2, max.nc=6, method = "kmeans", index = "alllong")
#fviz_nbclust(resnumclust2)

dfresnumclustS = data.frame("Clusters" = c(2,3,4,5,6),"Freq" = c(2,10,4,6,5))
dfresnumclustS$Freq = as.factor(dfresnumclustS$Freq)
dfresnumclustS$Clusters = as.factor(dfresnumclustS$Clusters)

#Plot Index
clIG = ggplot(data =dfresnumclustS, aes(x = Clusters, y = Freq, fill = Freq))+
  geom_bar(position = 'dodge', stat='identity', width=.5)+ theme_light(base_size = 13)+
  scale_y_discrete(position = "right")+
  ylab("Freq. among all Indices")+ xlab("Optimal Number of Clusters")

clIG = change_palette(clIG, palette= viridis(10)[1:5])
clIG = clIG + theme(legend.title = element_blank(),legend.position = "none",plot.title = element_text(hjust = 0.5))
clIG
ggsave()
