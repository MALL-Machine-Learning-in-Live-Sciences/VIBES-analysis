require(ggplot2)
require(viridis)
library("ggVennDiagram")
require(plotly)
require(ggpubr)
require(grid)

#P1_SRINIV

S_counts = readRDS(file = "projects/Entropy/data/benchmarks/Sriniv_Amsel_Genus_C_Benchmarks.rds")
S_AR = readRDS(file = "projects/Entropy/data/benchmarks/Sriniv_Amsel_Genus_AR_Benchmarks.rds")

#P1
## Fig1
# A) PLOTS COUNTS vs AR
# COUNTS
sclust = as.data.frame(S_counts$Bmr_Sriniv_Amsel_Genus_C_train_CLUST_3.rds)
sfcbf = as.data.frame(S_counts$Bmr_Sriniv_Amsel_Genus_C_train_FCBF_4.rds)
skw = as.data.frame(S_counts$Bmr_Sriniv_Amsel_Genus_C_train_KW_9.rds)
sldm = as.data.frame(S_counts$Bmr_Sriniv_Amsel_Genus_C_train_LDM_29.rds)
satal = list(sclust,sfcbf,skw,sldm)
nmes = c("CLUST", "FCBF", "KW","LDM")
names(satal) =nmes


for (i in seq_along(satal)){
  FSnames = c(rep(names(satal[i]),5))
  satal[[i]] = data.frame(FS =FSnames,Algorithm = c("RF", "GLMNET","xGBOOST", "SVM", "GBM"),
                          AUC = c(mean(satal[[i]]$auc[1:15]),mean(satal[[i]]$auc[16:30]),
                                  mean(satal[[i]]$auc[31:45]),mean(satal[[i]]$auc[46:60]),
                                  mean(satal[[i]]$auc[61:75])))
}

Sriniv_Counts = rbind(satal[[1]], satal[[2]], satal[[3]], satal[[4]])

# AR
sclust2 = as.data.frame(S_AR$Bmr_Sriniv_Amsel_Genus_AR_train_CLUST_5.rds)
sfcbf2 = as.data.frame(S_AR$Bmr_Sriniv_Amsel_Genus_AR_train_FCBF_3.rds)
skw2 = as.data.frame(S_AR$Bmr_Sriniv_Amsel_Genus_AR_train_KW_9.rds)
sldm2 = as.data.frame(S_AR$Bmr_Sriniv_Amsel_Genus_AR_train_LDM_24.rds)
satal2 = list(sclust2,sfcbf2,skw2,sldm2)
nmes = c("CLUST", "FCBF", "KW","LDM")
names(satal2) =nmes

for (i in seq_along(satal2)){
  FSnames = c(rep(names(satal2[i]),5))
  satal2[[i]] = data.frame(FS =FSnames,Algorithm = c("RF", "GLMNET","xGBOOST", "SVM", "GBM"),
                           AUC = c(mean(satal2[[i]]$auc[1:15]),mean(satal2[[i]]$auc[16:30]),
                                   mean(satal2[[i]]$auc[31:45]),mean(satal2[[i]]$auc[46:60]),
                                   mean(satal2[[i]]$auc[61:75])))
}

Sriniv_AR = rbind(satal2[[1]], satal2[[2]], satal2[[3]], satal2[[4]]) 


# PLOTS
plot_Sriniv_C = ggplot(Sriniv_Counts, aes(x = AUC, y = FS, color = Algorithm, shape = Algorithm)) +
  geom_point(size = 3) +
  scale_color_manual(values = viridis(5))+
  theme_light()+
  theme( axis.text=element_text(size=10),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))+
  ggtitle("Benchmark Counts")+
  theme(plot.title = element_text(hjust = 0.5))

l1 <- get_legend(plot_Sriniv_C)
plot_Sriniv_C = plot_Sriniv_C +theme(legend.position = "none")
plot_Sriniv_C

plot_Sriniv_AR = ggplot(Sriniv_AR, aes(x = AUC, y = FS, color = Algorithm, shape = Algorithm)) +
  geom_point(size = 3) +
  scale_color_manual(values = viridis(5))+
  theme_light()+
  theme( axis.text=element_text(size=10),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))+
  ggtitle("Benchmark Relative Abundance")+
  theme(plot.title = element_text(hjust = 0.5))
plot_Sriniv_AR



#B) PLOTS COUNTS vs AR (FEATURES)
#Counts
SKW_C = readRDS("projects/Entropy/data/train/Sriniv_Amsel_Genus_C_train_KW_9.rds")
SFCBF_C = readRDS("projects/Entropy/data/train/Sriniv_Amsel_Genus_C_train_FCBF_4.rds")
SLDM_C = readRDS("projects/Entropy/data/train/Sriniv_Amsel_Genus_C_train_LDM_29.rds")
Slist_C = list( names(SKW_C), names(SFCBF_C), names(SLDM_C))
names = c("KW", "FCBF","LDM")
remove = "target"
names(Slist_C) = names
for (i in seq_along(Slist_C)){
  Slist_C[[i]] = Slist_C[[i]][!Slist_C[[i]] %in% remove]
}
#AR
#CLUST_AR = readRDS("projects/Entropy/data/CLUST_AR_names.rds")
SKW_AR = readRDS("projects/Entropy/data/train/Sriniv_Amsel_Genus_AR_train_KW_9.rds")
SFCBF_AR = readRDS("projects/Entropy/data/train/Sriniv_Amsel_Genus_AR_train_FCBF_3.rds")
SLDM_AR = readRDS("projects/Entropy/data/train/Sriniv_Amsel_Genus_AR_train_LDM_24.rds")
Slist_AR = list( names(SKW_AR), names(SFCBF_AR), names(SLDM_AR))
names(Slist_AR) = names
for (i in seq_along(Slist_AR)){
  Slist_AR[[i]] = Slist_AR[[i]][!Slist_AR[[i]] %in% remove]
}


# PLOTS
Svenn_C  = ggVennDiagram(Slist_C,color ="#21908CFF",lty = 1, label_alpha = 0.5, label= "count") 
Svenn_C  = Svenn_C  + scale_fill_gradient(high=viridis(1),low = "#FDE725FF")
Svenn_C  = Svenn_C  + labs(fill = "Nº Features")+ ggtitle("Features Intersection")+theme(plot.title = element_text(hjust = 0.5))
l2 = get_legend(Svenn_C)
Svenn_C = Svenn_C+theme(legend.position = "none")
Svenn_C

Svenn_AR  = ggVennDiagram(Slist_AR,color = "#21908CFF",lty = 1, label_alpha = 0.5, label= "count" ) 
Svenn_AR  = Svenn_AR  + scale_fill_gradient(high=viridis(1),low = "#FDE725FF")
Svenn_AR  = Svenn_AR  + labs(fill = "Nº Features")+ggtitle("Features Intersection")+ theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
Svenn_AR


# C) Comparacion de los mejores modelos de cada tipo de datos
require(plyr)
comp.models= function(dataset){
  reps = rep(c(rep('1-rep', 5), rep('2-rep', 5), rep('3-rep', 5)))
  data = cbind.data.frame(dataset, reps)
  data <- subset( data, select = -iter )
  aggsdata = ddply(data, .(reps), summarize, AUC = mean(auc), ACC = mean(acc), MMCE = mean(mmce))
  return(aggsdata)
}
#COUNTS
SCLUST.GLMNET =as.data.frame(S_counts$Bmr_Sriniv_Amsel_Genus_C_train_CLUST_3.rds$results$dataset$classif.glmnet.tuned$measures.test)
SFCBF.GLMNET = as.data.frame(S_counts$Bmr_Sriniv_Amsel_Genus_C_train_FCBF_4.rds$results$dataset$classif.glmnet.tuned$measures.test)
SKW.GLMNET = as.data.frame(S_counts$Bmr_Sriniv_Amsel_Genus_C_train_KW_9.rds$results$dataset$classif.glmnet.tuned$measures.test) 
SLDM.RF = as.data.frame(S_counts$Bmr_Sriniv_Amsel_Genus_C_train_LDM_29.rds$results$dataset$classif.randomForest.tuned$measures.test)

l = list(SCLUST.GLMNET, SFCBF.GLMNET, SKW.GLMNET, SLDM.RF)
names = c("CLUST.GLMNET", "FCBF.GLMNET", "KW.GLMNET", "LDM.RF")
names(l)= names
for (i in seq_along(l)){
  l[[i]] = comp.models(l[[i]])
}
for (i in seq_along(l)){
  Model = c(rep(names(l[i]),3))
  l[[i]] = cbind(l[[i]], Model)
}

Sri_C_df = rbind.data.frame(l[[1]],l[[2]], l[[3]], l[[4]])

my_comparaisons <- list(c("CLUST.GLMNET", "FCBF.GLMNET"),
                        c("CLUST.GLMNET", "KW.GLMNET"),
                        c("CLUST.GLMNET", "LDM.RF"),
                        c("FCBF.GLMNET","KW.GLMNET"),
                        c("FCBF.GLMNET","LDM.RF"),
                        c("KW.GLMNET","LDM.RF"))

Sriniv_mod_C =  ggplot(Sri_C_df, aes(x = Model, y = AUC, color = Model))+
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
Sriniv_mod_C

# Abundancias relativas
SAR.CLUST.SVM =as.data.frame(S_AR$Bmr_Sriniv_Amsel_Genus_AR_train_CLUST_5.rds$results$dataset$classif.ksvm.tuned$measures.test)
SAR.FCBF.GLMNET = as.data.frame(S_AR$Bmr_Sriniv_Amsel_Genus_AR_train_FCBF_3.rds$results$dataset$classif.glmnet.tuned$measures.test)
SAR.KW.RF = as.data.frame(S_AR$Bmr_Sriniv_Amsel_Genus_AR_train_KW_9.rds$results$dataset$classif.randomForest.tuned$measures.test)
SAR.LDM.RF = as.data.frame(S_AR$Bmr_Sriniv_Amsel_Genus_AR_train_LDM_24.rds$results$dataset$classif.randomForest.tuned$measures.test)

l = list(SAR.CLUST.SVM, SAR.FCBF.GLMNET, SAR.KW.RF, SAR.LDM.RF)
names = c("CLUST.SVM", "FCBF.GLMNET", "KW.RF", "LDM.RF")
names(l)= names
for (i in seq_along(l)){
  l[[i]] = comp.models(l[[i]])
}
for (i in seq_along(l)){
  Model = c(rep(names(l[i]),3))
  l[[i]] = cbind(l[[i]], Model)
}

Sri_AR_df = rbind.data.frame(l[[1]],l[[2]], l[[3]], l[[4]])

my_comparaisons <- list(c("CLUST.SVM", "FCBF.GLMNET"),
                        c("CLUST.SVM", "KW.RF"),
                        c("CLUST.SVM", "LDM.RF"),
                        c("FCBF.GLMNET","KW.RF"),
                        c("FCBF.GLMNET","LDM.RF"),
                        c("KW.RF","LDM.RF"))

Sriniv_mod_AR =  ggplot(Sri_AR_df, aes(x = Model, y = AUC, color = Model))+
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
Sriniv_mod_AR

p1_Sriniv = ggarrange(plot_Sriniv_C, plot_Sriniv_AR, Svenn_C,Svenn_AR,Sriniv_mod_C,Sriniv_mod_AR,
              ncol = 2, nrow = 3,widths = c(4,4,1),
              labels = list("A", "", "B", "", "C"))
#p1_Sriniv = p1_Sriniv + theme(plot.margin = unit(c(1.2, 0.5, 1.2, 0.5), "cm"))
p1_Sriniv

#P2
require(mlr)
# A) Feature importance
#Counts
models = as.data.frame(S_counts$Bmr_Sriniv_Amsel_Genus_C_train_LDM_29.rds)
model = models[models$learner.id == "classif.randomForest.tuned",]
index = which.max(model$auc)
# Get best model
modelos = getBMRModels(S_counts$Bmr_Sriniv_Amsel_Genus_C_train_LDM_29.rds)
alg_index = which(names(modelos$dataset)=="classif.randomForest.tuned")
best_mod_C = getLearnerModel(modelos$dataset[[alg_index]][[index]])
FI_C = getFeatureImportance(best_mod_C)
View(FI_C$res)

#AR
models = as.data.frame(S_AR$Bmr_Sriniv_Amsel_Genus_AR_train_LDM_24.rds)
model = models[models$learner.id == "classif.randomForest.tuned",]
index = which.max(model$auc)

# Get best model
modelos = getBMRModels(S_AR$Bmr_Sriniv_Amsel_Genus_AR_train_LDM_24.rds)
alg_index = which(names(modelos$dataset)=="classif.randomForest.tuned")
best_mod_AR = getLearnerModel(modelos$dataset[[alg_index]][[index]])
FI_AR = getFeatureImportance(best_mod_AR)
View(FI_AR$res)

#B) Validation
source("git/Entropy/functions/FunctionsGetSplitData.R")
test_C = readRDS("projects/Entropy/data/test/Sriniv_Amsel_Genus_C_test.rds")
test_AR = readRDS("projects/Entropy/data/test/Sriniv_Amsel_Genus_AR_test.rds")
# Counts
#We have to retain only the features used in the model
features = best_mod_C$features
features = c(features, "target")
test_C = subset(test_C, select= features)
test_C= norm.dataset(test_C)
# Make task
test_task = makeClassifTask(data = test_C, target = "target")
test_task = normalizeFeatures(
  test_task,
  method = "range",
  cols = NULL,
  range = c(0, 1),
  on.constant = "quiet")
prediccion_C <- predict(best_mod_C, task= test_task, type = "prob")
library(caret)
confusionMatrix(data = prediccion_C$data$response, reference = as.factor(test_C$target))
a = asROCRPrediction(prediccion_C)
p = ROCR::performance(a, "tpr", "fpr")
plot(p)
#AR
#We have to retain only the features used in the model
features = best_mod_AR$features
features = c(features, "target")
test_AR = subset(test_AR, select= features)
test_AR= norm.dataset(test_AR)
# Make task
test_task = makeClassifTask(data = test_AR, target = "target")
test_task = normalizeFeatures(
  test_task,
  method = "range",
  cols = NULL,
  range = c(0, 1),
  on.constant = "quiet")
prediccionAR <- predict(best_mod_AR, task= test_task, type = "prob")
library(caret)
confusionMatrix(data = prediccionAR$data$response, reference = as.factor(test_AR$target))
a = asROCRPrediction(prediccionAR)
p2 = ROCR::performance(a, "tpr", "fpr")
plot(p2)



