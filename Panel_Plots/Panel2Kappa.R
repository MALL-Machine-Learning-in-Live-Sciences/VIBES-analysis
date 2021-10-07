#Panel2
require(plotly)
library(mlr)
library(ggplot2)
library(viridis)

#Genus
#Cluster Visualization
a = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Cluster/3C/FS/Ravel_3C_3.rds")
clusterG = plot_ly(x = a$Prevotella, 
                   y = a$Dialister, 
                   z = a$Sneathia,
                   type = "scatter3d", 
                   mode = "markers",colors =  c("#FDE725FF","#440154FF","#21908CFF" ),
                   color = as.factor(a$cluster)) %>%
  layout(legend = list(x = 0.8, y = 0.5,size = 30),
         scene = list(xaxis = list(title = "Prevotella"),
                      yaxis = list(title = "Dialister"),
                      zaxis = list(title = "Sneathia")))
clusterG

# Validation Feat-feat
val = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Cluster/3C/VAL1-23.rds")

Kappa = list()
nomb = names(val$CM)
for (i in seq_along(val$CM)){
  Kappa[[i]] = val$CM[[i]]$overall[[2]]
}
df = data.frame("N.Feat" = nomb, "Kappa" = unlist(Kappa))
df$N.Feat = as.double(substr(df$N.Feat,10,12))
df = df[order(df$N.Feat),]


p = ggplot(df, aes(x = N.Feat, y = Kappa))+
  geom_line(color= viridis(3)[2]) + geom_point(color= viridis(1))+
  geom_segment(aes(x = 7, y = 0.92, xend = 7.8, yend = df$Kappa[[8]]+0.005),
               arrow = arrow(length = unit(0.3, "cm")),
               colour = viridis(3)[3])+
  geom_point(aes(x=8, y=df$Kappa[[8]]), shape=23,fill=viridis(3)[1],color=viridis(3)[3],size=3)+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.ticks.x = element_blank())
p


#Feature Importance Genera
kk = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Cluster/3C/BMR1-23.rds")
#Change BMr
models = getBMRModels(kk[[22]])
rf.vi = models[[1]]$classif.randomForest.tuned
iters = list()
for (i in 1:length(rf.vi)) {
  iters[[i]] = as.data.frame(getFeatureImportance(rf.vi[[i]])$res)$importance
}
names = as.data.frame(getFeatureImportance(rf.vi[[1]])$res)$variable
importance = as.numeric(Reduce("+", iters))
ff = data.frame(cbind(names,importance))
ff$importance = as.numeric(ff$importance)
ff = ff[order(ff$importance,decreasing = TRUE),]
ff
plotFIC = ggplot(ff, aes(x = reorder(names, importance), y = importance))+
  geom_segment( aes(xend=names, yend=0,), color = viridis(3)[2]) +
  geom_point( size=2, color=viridis(1)) +
  coord_flip() +
  theme_light()+
  theme( axis.text=element_text(size=9),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))
plotFIC


#Species
#Cluster Visualization
a = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/FS/Ravel_3C_3.rds")
clusterS = plot_ly(x = a$Lactobacillus.crispatus, 
                   y = a$Fannyhessea.vaginae, 
                   z = a$Sneathia.sanguinegens,
                   type = "scatter3d", 
                   mode = "markers",colors =  c("#21908CFF","#440154FF","#FDE725FF" ),
                   color = as.factor(a$cluster)) %>%
  layout(legend = list(x = 0.8, y = 0.5,size = 30),
         scene = list(xaxis = list(title = "Lactobacillus.crispatus"),
                      yaxis = list(title = "Fannyhessea.vaginae"),
                      zaxis = list(title = "Sneathia.sanguinegens")))
clusterS


# Validation Feat a Feat
val = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/VAL1-27.rds")

Kappa = list()
nomb = names(val$CM)
for (i in seq_along(val$CM)){
  Kappa[[i]] = val$CM[[i]]$overall[[2]]
}
df = data.frame("N.Feat" = nomb, "Kappa" = unlist(Kappa))
df$N.Feat = as.double(substr(df$N.Feat,10,12))
df = df[order(df$N.Feat),]


p = ggplot(df, aes(x = N.Feat, y = Kappa))+
  geom_line(color= viridis(3)[2]) + geom_point(color= viridis(1))+
  geom_segment(aes(x = 7, y = 0.98, xend = 7.8, yend = df$Kappa[[8]]+0.005),
               arrow = arrow(length = unit(0.3, "cm")),
               colour = viridis(3)[3])+
  geom_point(aes(x=8, y=df$Kappa[[8]]), shape=23,fill=viridis(3)[1],color=viridis(3)[3],size=3)+
  theme_light()+
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.ticks.x = element_blank())+ scale_y_continuous(position = 'right')
p

#FI
kk = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/BMR1-27.rds")
#Change BMr
models = getBMRModels(kk[[26]])
rf.vi = models[[1]]$classif.randomForest.tuned
iters = list()
for (i in 1:length(rf.vi)) {
  iters[[i]] = as.data.frame(getFeatureImportance(rf.vi[[i]])$res)$importance
}
names = as.data.frame(getFeatureImportance(rf.vi[[1]])$res)$variable
importance = as.numeric(Reduce("+", iters))
ff = data.frame(cbind(names,importance))
ff$importance = as.numeric(ff$importance)
ff = ff[order(ff$importance,decreasing = TRUE),]
ff
plotFIC = ggplot(ff, aes(x = reorder(names, importance), y = importance))+
  geom_segment( aes(xend=names, yend=0,), color = viridis(3)[2]) +
  geom_point( size=2, color=viridis(1)) +
  scale_x_discrete(position = 'top')+scale_y_reverse()+
  coord_flip() +
  theme_light()+
  theme( axis.text=element_text(size=9),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10)) 
plotFIC

