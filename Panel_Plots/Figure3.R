#Figure 3
require(plotly)
library(mlr)
library(ggplot2)
library(viridis)

# Validation Feat-feat
val= readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Cluster/3C/VAL1-23.rds")

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
  geom_point(aes(x=8, y=df$Kappa[[8]]), shape=23,fill=viridis(3)[1],color=viridis(3)[3],size=3)+ ylim(0,1)+
  theme_light(base_size = 13)+
  theme( plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.ticks.x = element_blank())+
  scale_y_continuous(limits = c(0.3,1),breaks = c(0.4,0.6,0.8,1))
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
ff$names[ff$names == "Fannyhessea"] <- "Atopobium"
  
plotFIC = ggplot(ff, aes(x = reorder(names, importance), y = importance))+
  geom_segment( aes(xend=names, yend=0,), color = viridis(3)[2]) +
  geom_point( size=2, color=viridis(1)) +
  coord_flip() +
  theme_light(base_size = 16)+
  theme( axis.text=element_text(size=9),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))
plotFIC


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
  theme_light(base_size = 13)+
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.ticks.x = element_blank())+ scale_y_continuous(limits = c(0.3,1),breaks = c(0.4,0.6,0.8,1),position = 'right')
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
ff$names[ff$names == "Fannyhessea.vaginae"] <- "Atopobium.vaginae"

plotFIC = ggplot(ff, aes(x = reorder(names, importance), y = importance))+
  geom_segment( aes(xend=names, yend=0,), color = viridis(3)[2]) +
  geom_point( size=2, color=viridis(1)) +
  scale_x_discrete(position = 'top')+scale_y_reverse()+
  coord_flip() +
  theme_light(base_size = 16)+
  theme( axis.text=element_text(size=9),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10)) 
plotFIC

