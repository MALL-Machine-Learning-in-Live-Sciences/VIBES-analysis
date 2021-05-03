require(ggplot2)
require(viridis)
library("ggVennDiagram")
require(plotly)
require(ggpubr)
require(grid)

counts = readRDS(file = "projects/Entropy/data/benchmarks/Ravel_Genus_C_Benchmarks.rds")
ar = readRDS(file = "projects/Entropy/data/benchmarks/Ravel_Genus_AR_Benchmarks.rds")

#P1
## Fig1
# A) PLOTS COUNTS vs AR
# COUNTS
dclust = as.data.frame(counts$Bmr_Ravel_Genus_C_train_CLUST_7.rds)
dfcbf = as.data.frame(counts$Bmr_Ravel_Genus_C_train_FCBF_8.rds)
dkw = as.data.frame(counts$Bmr_Ravel_Genus_C_train_KW_12.rds)
dldm = as.data.frame(counts$Bmr_Ravel_Genus_C_train_LDM_40.rds)
datal = list(dclust,dfcbf,dkw,dldm)
nmes = c("CLUST", "FCBF", "KW","LDM")
names(datal) =nmes


for (i in seq_along(datal)){
  FSnames = c(rep(names(datal[i]),5))
  datal[[i]] = data.frame(FS =FSnames,Algorithm = c("RF", "GLMNET","xGBOOST", "SVM", "GBM"),
                      AUC = c(mean(datal[[i]]$auc[1:15]),mean(datal[[i]]$auc[16:30]),
                              mean(datal[[i]]$auc[31:45]),mean(datal[[i]]$auc[46:60]),
                              mean(datal[[i]]$auc[61:75])))
}

dataCounts = rbind(datal[[1]], datal[[2]], datal[[3]], datal[[4]])

# AR
dclust2 = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_CLUST_7.rds)
dfcbf2 = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds)
dkw2 = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_KW_12.rds)
dldm2 = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_LDM_39.rds)
datal2 = list(dclust2,dfcbf2,dkw2,dldm2)
nmes = c("CLUST", "FCBF", "KW","LDM")
names(datal2) =nmes

for (i in seq_along(datal2)){
  FSnames = c(rep(names(datal2[i]),5))
  datal2[[i]] = data.frame(FS =FSnames,Algorithm = c("RF", "GLMNET","xGBOOST", "SVM", "GBM"),
                          AUC = c(mean(datal2[[i]]$auc[1:15]),mean(datal2[[i]]$auc[16:30]),
                                  mean(datal2[[i]]$auc[31:45]),mean(datal2[[i]]$auc[46:60]),
                                  mean(datal2[[i]]$auc[61:75])))
}

dataAR = rbind(datal2[[1]], datal2[[2]], datal2[[3]], datal2[[4]]) 


# PLOTS
comp_C = ggplot(dataCounts, aes(x = AUC, y = FS, color = Algorithm, shape = Algorithm)) +
  geom_point(size = 3) +
  scale_color_manual(values = viridis(5))+
  theme_light()+
  theme( axis.text=element_text(size=10),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))+
  ggtitle("Benchmark Counts")+
  theme(plot.title = element_text(hjust = 0.5))
  
l1 <- get_legend(comp_C)
comp_C = comp_C +theme(legend.position = "none")
comp_C

comp_AR = ggplot(dataAR, aes(x = AUC, y = FS, color = Algorithm, shape = Algorithm)) +
  geom_point(size = 3) +
  scale_color_manual(values = viridis(5))+
  theme_light()+
  theme( axis.text=element_text(size=10),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))+
  ggtitle("Benchmark Relative Abundance")+
  theme(plot.title = element_text(hjust = 0.5))
comp_AR

# B) PLOTS COUNTS vs AR (FEATURES)
#Counts
#CLUST_C = readRDS("projects/Entropy/data/CLUST_C_names.rds")
KW_C = readRDS("projects/Entropy/data/train/Ravel_Genus_C_train_KW_12.rds")
FCBF_C = readRDS("projects/Entropy/data/train/Ravel_Genus_C_train_FCBF_8.rds")
LDM_C = readRDS("projects/Entropy/data/train/Ravel_Genus_C_train_LDM_40.rds")
list_C = list( names(KW_C), names(FCBF_C), names(LDM_C))
names = c("KW", "FCBF","LDM")
remove = "target"
names(list_C) = names
for (i in seq_along(list_C)){
  list_C[[i]] = list_C[[i]][!list_C[[i]] %in% remove]
}
#AR
#CLUST_AR = readRDS("projects/Entropy/data/CLUST_AR_names.rds")
KW_AR = readRDS("projects/Entropy/data/train/Ravel_Genus_AR_train_KW_12.rds")
FCBF_AR = readRDS("projects/Entropy/data/train/Ravel_Genus_AR_train_FCBF_9.rds")
LDM_AR = readRDS("projects/Entropy/data/train/Ravel_Genus_AR_train_LDM_39.rds")
list_AR = list( names(KW_AR), names(FCBF_AR), names(LDM_AR))
names(list_AR) = names
for (i in seq_along(list_AR)){
  list_AR[[i]] = list_AR[[i]][!list_AR[[i]] %in% remove]
}


# PLOTS
venn_C  = ggVennDiagram(list_C,color ="#21908CFF",lty = 1, label_alpha = 0.5, label= "count") 
venn_C  = venn_C  + scale_fill_gradient(high=viridis(1),low = "#FDE725FF")
venn_C  = venn_C  + labs(fill = "Nº Features")+ ggtitle("Features Intersection")+theme(plot.title = element_text(hjust = 0.5))
l2 = get_legend(venn_C)
venn_C = venn_C+theme(legend.position = "none")
venn_C

venn_AR  = ggVennDiagram(list_AR,color = "#21908CFF",lty = 1, label_alpha = 0.5, label= "count" ) 
venn_AR  = venn_AR  + scale_fill_gradient(high=viridis(1),low = "#FDE725FF")
venn_AR  = venn_AR  + labs(fill = "Nº Features")+ggtitle("Features Intersection")+ theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
venn_AR

# C) Comparacion de los mejores modelos de cada tipo de datos
require(plyr)
comp.models= function(dataset){
  reps = rep(c(rep('1-rep', 5), rep('2-rep', 5), rep('3-rep', 5)))
  data = cbind.data.frame(dataset, reps)
  data <- subset( data, select = -iter )
  aggsdata = ddply(data, .(reps), summarize, AUC = mean(auc), ACC = mean(acc), MMCE = mean(mmce))
  return(aggsdata)
}

CLUST.GBM =as.data.frame(counts$Bmr_Ravel_Genus_C_train_CLUST_7.rds$results$dataset$classif.gbm.tuned$measures.test)
FCBF.SVM = as.data.frame(counts$Bmr_Ravel_Genus_C_train_FCBF_8.rds$results$dataset$classif.ksvm.tuned$measures.test)
KW.GLMNET = as.data.frame(counts$Bmr_Ravel_Genus_C_train_KW_12.rds$results$dataset$classif.glmnet.tuned$measures.test) 
LDM.RF = as.data.frame(counts$Bmr_Ravel_Genus_C_train_LDM_40.rds$results$dataset$classif.randomForest.tuned$measures.test)

l = list(CLUST.GBM, FCBF.SVM, KW.GLMNET, LDM.RF)
names = c("CLUST.GBM", "FCBF.SVM", "KW.GLMNET", "LDM.RF")
names(l)= names
for (i in seq_along(l)){
  l[[i]] = comp.models(l[[i]])
}
for (i in seq_along(l)){
  Model = c(rep(names(l[i]),3))
  l[[i]] = cbind(l[[i]], Model)
}

df = rbind.data.frame(l[[1]],l[[2]], l[[3]], l[[4]])

my_comparaisons <- list(c("CLUST.GBM", "FCBF.SVM"),
                        c("CLUST.GBM", "KW.GLMNET"),
                        c("CLUST.GBM", "LDM.RF"),
                        c("FCBF.SVM","KW.GLMNET"),
                        c("FCBF.SVM","LDM.RF"),
                        c("KW.GLMNET","LDM.RF"))

mod_C =  ggplot(df, aes(x = Model, y = AUC, color = Model))+
  geom_jitter()+
  geom_boxplot(aes(fill = Model), notch = F)+
  scale_fill_manual(values = viridis(4))+
  theme_light()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text=element_text(size=10), axis.title.y = element_text(size = 10),
        axis.ticks = element_blank(),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10))+
  ggtitle("Best models Counts")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_compare_means(comparisons = my_comparaisons, bracket.size = 0.3,  size = 3)
mod_C

# Abundancias relativas
CLUST.RF =as.data.frame(ar$Bmr_Ravel_Genus_AR_train_CLUST_7.rds$results$dataset$classif.randomForest.tuned$measures.test)
FCBF.GBM = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds$results$dataset$classif.gbm.tuned$measures.test)
KW.GBM = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_KW_12.rds$results$dataset$classif.gbm.tuned$measures.test)
LDM.RF = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_LDM_39.rds$results$dataset$classif.randomForest.tuned$measures.test)

l = list(CLUST.RF, FCBF.GBM, KW.GBM, LDM.RF)
names = c("CLUST.RF", "FCBF.GBM", "KW.GBM", "LDM.RF")
names(l)= names
for (i in seq_along(l)){
  l[[i]] = comp.models(l[[i]])
}
for (i in seq_along(l)){
  Model = c(rep(names(l[i]),3))
  l[[i]] = cbind(l[[i]], Model)
}

df = rbind.data.frame(l[[1]],l[[2]], l[[3]], l[[4]])

my_comparaisons <- list(c("CLUST.RF", "FCBF.GBM"),
                        c("CLUST.RF", "KW.GBM"),
                        c("CLUST.RF", "LDM.RF"),
                        c("FCBF.GBM","KW.GBM"),
                        c("FCBF.GBM","LDM.RF"),
                        c("KW.GBM","LDM.RF"))

mod_AR =  ggplot(df, aes(x = Model, y = AUC, color = Model))+
  geom_jitter()+
  geom_boxplot(aes(fill = Model), notch = F)+
  scale_fill_manual(values = viridis(4))+
  theme_light()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text=element_text(size=10), axis.title.y = element_text(size = 10),
        axis.ticks = element_blank(),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10))+
  ggtitle("Best models Relative Abundance")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_compare_means(comparisons = my_comparaisons, bracket.size = 0.3,  size = 3)
mod_AR

#PANEL
p = ggarrange(comp_C, comp_AR, venn_C,venn_AR,mod_C,mod_AR,
              ncol = 2, nrow = 3,widths = c(4,4,1),
              labels = list("A", "", "B", "", "C"))
p = p + theme(plot.margin = unit(c(1.2, 0.5, 1.2, 0.5), "cm"))
p

#P2
#P2
require(mlr)
# A) Feature importance
#Counts
models = as.data.frame(counts$Bmr_Ravel_Genus_C_train_LDM_40.rds)
model = models[models$learner.id == "classif.randomForest.tuned",]
index = which.max(model$auc)
# Get best model
modelos = getBMRModels(counts$Bmr_Ravel_Genus_C_train_LDM_40.rds)
alg_index = which(names(modelos$dataset)=="classif.randomForest.tuned")
best_mod_C = getLearnerModel(modelos$dataset[[alg_index]][[index]])
FI_C = getFeatureImportance(best_mod_C)
View(FI_C$res)

#AR
models = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds)
model = models[models$learner.id == "classif.gbm.tuned",]
index = which.max(model$auc)

# Get best model
modelos = getBMRModels(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds)
alg_index = which(names(modelos$dataset)=="classif.gbm.tuned")
best_mod_AR = getLearnerModel(modelos$dataset[[alg_index]][[index]])
FI_AR = getFeatureImportance(best_mod_AR)
View(FI_AR$res)

## Fig2
# AUROC Curve Predict
source("git/Entropy/functions/FunctionsDataFilter.R")
source("git/Entropy/functions/FunctionsValidation.R")
source("git/Entropy/functions/FunctionsGetSplitData.R")
source("git/Entropy/functions/FunctionsML.R")
require(phyloseq)
# Load phyloseqs
Ravel = readRDS("projects/Entropy/data/Ravel_phyloseq.rds")
Ravel = phy.aglomerate(phyobject = Ravel, rank = "Rank6")
Sriniv = readRDS("projects/Entropy/data/Sriniv_Nugent_phyloseq.rds")
Sriniv = relat.abun(Sriniv)
Sriniv = phy.aglomerate(phyobject = Sriniv, rank = "Rank6")

best= readRDS(file = "projects/Entropy/data/models/RetrainBestGBMFCBF7Feat")
features = best$features

# Saco las correspondencias entre el genero y number acceso
Ravel.df = as.data.frame(tax_table(Ravel))
Ravel.df = Ravel.df[features,]
#Creo una tabla de igualdad entre number acces y generos
features = as.data.frame(cbind(rownames(Ravel.df),Ravel.df$Rank6))
# Selecciono los generos en el dataframe de Sriniv
Sriniv_sub <- subset_taxa(Sriniv, Rank6 %in% c(features[,2])) 
Sriniv.df = as.data.frame(tax_table(Sriniv_sub))

test.Sriniv = get.sriniv.test(Sriniv_sub =Sriniv_sub,Sriniv.df = Sriniv.df, features = features  )
# Make task
test_task = makeClassifTask(data = test.Sriniv, target = "target")
test_task = normalizeFeatures(
  test_task,
  method = "range",
  cols = NULL,
  range = c(0, 1),
  on.constant = "quiet")

prediccion <- predict(best, task= test_task, type = "prob")
library(caret)
table(test.Sriniv$target,prediccion$data$response)
confusionMatrix(data = prediccion$data$response, reference = as.factor(test.Sriniv$target))
a = asROCRPrediction(prediccion)
p = ROCR::performance(a, "tpr", "fpr")
plot(p)

## Fig3
# Heatmap Phyloseq



