#Panel1
library(UpSetR)
library(ComplexHeatmap)
library(viridis)
library(ggplot2)
library(ggpubr)

#Genus
a = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Benchmarks/BMR_Ge.rds")
# Intersection Features UPSET Genus
lt = list(FCBF = a$FCBF$results$FCBF$classif.randomForest.tuned$models[[1]]$features,
          KW = a$KW$results$KW$classif.randomForest.tuned$models[[1]]$features,
          LDM = a$LDM$results$LDM$classif.randomForest.tuned$models[[1]]$features,
          RAW = a$Raw$results$Raw$classif.randomForest.tuned$models[[1]]$features)
list_to_matrix(lt)
up_plotG = upset(fromList(lt), order.by = "freq",sets.x.label = "Number of Genera",mainbar.y.label = "Number of Genera Intersected",
                 decreasing = FALSE,sets.bar.color = viridis(20)[17:20],text.scale = 1.5,
                 main.bar.color = "black",point.size = 4,line.size = 1,
                 queries = list(
                   list(query = intersects, params = list("FCBF", "LDM", "RAW"), color = viridis(20)[1],active = T),
                   list(query = intersects, params = list("RAW"), color = viridis(20)[5], active = T),
                   list(query = intersects, params = list("FCBF", "KW","LDM", "RAW"), color = viridis(20)[10],active = T),
                   list(query = intersects, params = list("KW", "LDM", "RAW"), color = viridis(20)[15],active = T),
                   list(query = intersects, params = list("RAW", "LDM"), color = viridis(20)[18],active = T)
                 )) 
up_plotG


#Pointplot External Validation Genus
l = list.files("~/git/BVMetaGenomics/data/GenusIntersect/Validation/")
remove <- c ("HIL")
l =  l [! l %in% remove]
ldata = list()
for (i in 1:length(l)) {
  ldata[[i]] = readRDS(paste("~/git/BVMetaGenomics/data/GenusIntersect/Validation/", l[[i]], sep=''))
}
names(ldata) = l
Algoritmo = rep(c(rep('GBM', 4), rep('GLMNET', 4), rep('RF', 4),rep('SVM', 4),rep('xGBOOST', 4)))
FS = rep(c("FCBF", "KW", "LDM", "RAW"),5)
Kappa = list()
for (i in 1:length(l)){
  Kappa[[i]] = ldata[[i]]$overall[[2]]
}
Kappa = unlist(Kappa)

data  <- data.frame (Algoritmo  = Algoritmo,
                     FS = FS,
                     Kappa = Kappa)

p = ggplot(data) +
  geom_point(aes(x=Algoritmo, y=Kappa,color = FS,shape=FS,color = factor(FS)), size =3) + theme_light()+
  xlab(NULL) + coord_cartesian(ylim=c(0.4,1))+ theme(legend.title=element_blank())
p = change_palette(p, palette= viridis(4))
p = p + theme(legend.position = "none" )
p


#Pointplot External Validation HIL Genus
l = list.files("~/git/BVMetaGenomics/data/GenusIntersect/Validation/HIL/")
ldata = list()
for (i in 1:length(l)) {
  ldata[[i]] = readRDS(paste("~/git/BVMetaGenomics/data/GenusIntersect/Validation/HIL/", l[[i]], sep=''))
}
names(ldata) = l
Algoritmo = rep(c( rep('RF', 4)))
FS = rep(c("FCBF", "KW", "LDM", "RAW"))

Kappa = list()
for (i in 1:length(l)){
  Kappa[[i]] = ldata[[i]][[1]]$overall[[2]]
}
Kappa = unlist(Kappa)

data  <- data.frame (Algoritmo  = Algoritmo,
                     FS = FS,
                     Kappa = Kappa)

p2 = ggplot(data) +
  geom_point(aes(x=Algoritmo, y=Kappa,color = FS,shape=FS,color = factor(FS)), size =4) + theme_light()+
  xlab(NULL) + coord_cartesian(ylim=c(0.4,1))+ theme(legend.title=element_blank())
p2 = change_palette(p2, palette= viridis(4))
p2 = p2 + theme(legend.position = "none" )
p2

# Species
a = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Benchmarks/BMR_Sp.rds")
# Intersection Features
lt = list(FCBF = a$FCBF$results$FCBF$classif.randomForest.tuned$models[[1]]$features,
          KW = a$KW$results$KW$classif.randomForest.tuned$models[[1]]$features,
          LDM = a$LDM$results$LDM$classif.randomForest.tuned$models[[1]]$features,
          RAW = a$Raw$results$Raw$classif.randomForest.tuned$models[[1]]$features)
list_to_matrix(lt)
up_plotS = upset(fromList(lt), order.by = "freq",sets.x.label = "Number of Species",mainbar.y.label = "Number of Species Intersected",
                 decreasing = FALSE,sets.bar.color = viridis(20)[17:20],text.scale = 1.5,
                 main.bar.color = "black",point.size = 4,line.size = 1,
                 queries = list(
                   list(query = intersects, params = list("FCBF", "LDM", "RAW"), color = viridis(20)[1],active = T),
                   list(query = intersects, params = list("RAW"), color = viridis(20)[5], active = T),
                   list(query = intersects, params = list("FCBF", "KW","LDM", "RAW"), color = viridis(20)[10],active = T),
                   list(query = intersects, params = list("KW", "LDM", "RAW"), color = viridis(20)[15],active = T),
                   list(query = intersects, params = list("RAW", "LDM"), color = viridis(20)[18],active = T)
                 )) 
up_plotS 


#Pointplot External Validation Species
l = list.files("~/git/BVMetaGenomics/data/SpeciesIntersect/Validation/")
remove <- c ("HIL")
l =  l [! l %in% remove]
ldata = list()
for (i in 1:length(l)) {
  ldata[[i]] = readRDS(paste("~/git/BVMetaGenomics/data/SpeciesIntersect/Validation/", l[[i]], sep=''))
}
names(ldata) = l
Algoritmo = rep(c(rep('GBM', 4), rep('GLMNET', 4), rep('RF', 4),rep('SVM', 4),rep('xGBOOST', 4)))
FS = rep(c("FCBF", "KW", "LDM", "RAW"),5)
Kappa = list()
for (i in 1:length(l)){
  Kappa[[i]] = ldata[[i]]$overall[[2]]
}
Kappa = unlist(Kappa)

data  <- data.frame (Algoritmo  = Algoritmo,
                     FS = FS,
                     Kappa = Kappa)

p3 = ggplot(data) +
  geom_point(aes(x=Algoritmo, y=Kappa,color = FS,shape=FS,color = factor(FS)), size =3) + theme_light()+
  xlab(NULL) + coord_cartesian(ylim=c(0.4,1))+ theme(legend.title=element_blank())+scale_y_continuous(position = 'right')
p3 = change_palette(p3, palette= viridis(4))
p3 = p3 + theme(legend.position = "none" )
p3

#Pointplot External Validation HIL Species
l = list.files("~/git/BVMetaGenomics/data/SpeciesIntersect/Validation/HIL/")
ldata = list()
for (i in 1:length(l)) {
  ldata[[i]] = readRDS(paste("~/git/BVMetaGenomics/data/SpeciesIntersect/Validation/HIL/", l[[i]], sep=''))
}
names(ldata) = l
Algoritmo = rep(c( rep('RF', 4)))
FS = rep(c("FCBF", "KW", "LDM", "RAW"))

Kappa = list()
for (i in 1:length(l)){
  Kappa[[i]] = ldata[[i]][[1]]$overall[[2]]
}
Kappa = unlist(Kappa)

data  <- data.frame (Algoritmo  = Algoritmo,
                     FS = FS,
                     Kappa = Kappa)

p4 = ggplot(data) +
  geom_point(aes(x=Algoritmo, y=Kappa,color = FS,shape=FS,color = factor(FS)), size =4) + theme_light()+
  xlab(NULL) + coord_cartesian(ylim=c(0.4,1))+ theme(legend.title=element_blank())+ scale_y_continuous(position = 'right')
p4 = change_palette(p4, palette= viridis(4))
# Extract the legend. Returns a gtable
leg <- get_legend(p4)
# Convert to a ggplot and print
leg <- as_ggplot(leg)
p4 = p4 + theme(legend.position = "none" )
p4
