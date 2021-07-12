require(ggplot2)
require(viridis)
library("ggVennDiagram")
require(plotly)
require(ggpubr)
require(grid)
require(ggplotify)
library(easyGgplot2)
counts = readRDS(file = "projects/Entropy/data/benchmarks/Nugent/Ravel_Genus_C_Second_Benchmarks.rds")
ar = readRDS(file = "projects/Entropy/data/benchmarks/Nugent/Ravel_Genus_AR_Second_Benchmarks.rds")
two_ft = readRDS(file = "projects/Entropy/data/benchmarks/")
#P1
###### PAnel 1
# A) PLOTS COUNTS vs AR
# COUNTS
d2f = as.data.frame()
dfcbf = as.data.frame(counts$Bmr_Ravel_Genus_C_train_FCBF_8.rds)
dkw = as.data.frame(counts$Bmr_Ravel_Genus_C_train_KW_12.rds)
dldm = as.data.frame(counts$Bmr_Ravel_Genus_C_train_LDM_40.rds)

l = list(dfcbf,dkw, dldm )
fs = list("FCBF", "KW", "LDM")
lista= list()
comp.models= function(dataset, fs){
  require(plyr)
  reps = rep(c(rep('1-rep', 5), rep('2-rep', 5), rep('3-rep', 5),rep('4-rep', 5),rep('5-rep', 5),rep('6-rep', 5),rep('7-rep', 5),rep('8-rep', 5),rep('9-rep', 5),rep('10-rep', 5)),5)
  data = cbind.data.frame(dataset, reps)
  data <- subset( data, select = -iter )
  data$learner.id[data$learner.id == "classif.randomForest.tuned"] <- "RF"
  data$learner.id[data$learner.id == "classif.glmnet.tuned"] <- "GLMNET"
  data$learner.id[data$learner.id == "classif.xgboost.tuned"] <- "xGBOOST"
  data$learner.id[data$learner.id == "classif.ksvm.tuned"] <- "SVM"
  data$learner.id[data$learner.id == "classif.gbm.tuned"] <- "GBM"
  data$task.id = fs
  names(data)[names(data) == "task.id"] <- "FS"
  names(data)[names(data) == "learner.id"] <- "Algorithm"
  aggsdata = ddply(data, .(FS,Algorithm,reps ), summarize, AUC = mean(auc), ACC = mean(acc), MMCE = mean(mmce))
  return(aggsdata)
}

for (i in seq_along(l)){
  lista[[i]]= comp.models(l[[i]], fs[[i]])
}
df = rbind(lista[[1]], lista[[2]], lista[[3]])


bp <- ggplot(df, aes(x=FS, y=AUC, group=FS)) +
  geom_boxplot(aes(fill=FS)) +
  theme_light()+ theme(axis.title.y = element_blank(), axis.ticks.x = element_blank(),legend.position = "bottom",legend.title = element_blank(),
                                     axis.text.x = element_blank())+
  facet_wrap(~ Algorithm, scales ="fixed")+#coord_cartesian(ylim = c(0.5, 1))+
  #stat_compare_means(comparisons = list(c("FCBF", "KW"),c("FCBF", "LDM"), c("KW", "LDM")),method = "wilcox.test", bracket.size = 0.2,  size = 2, paired = FALSE,label = "p.signif",hide.ns = TRUE,vjust = 0.5)+
  theme(strip.background = element_blank(),strip.text.x = element_text(
    size = 8, color = "black"), axis.ticks = element_blank(),axis.title.x.bottom = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank())

bp =change_palette(bp, palette= viridis(3))
bp_counts = bp+ggtitle("Counts Benchmark")+theme(legend.position = c(0.82, 0.25),
                                          legend.text =  element_text(size=10),
                                          legend.key = element_rect(size = 6),
                                          legend.key.height = unit(1, "cm"),
                                          legend.key.width = unit(1, "cm"),
                                          plot.title = element_text(hjust = 0.5))



##### Fried Test #####                   
k_GBM = df[df$Algorithm == "GBM",]
k_GLMNET = df[df$Algorithm == "GLMNET",]
k_RF = df[df$Algorithm == "RF",]
k_SVM = df[df$Algorithm == "SVM",]
k_xGBOOST = df[df$Algorithm == "xGBOOST",]

ft_GBM = friedman.test(k_GBM$AUC, k_GBM$FS)
ft_GLMNET = friedman.test(k_GLMNET$AUC, k_GLMNET$FS, k_GLMNET$rep)
ft_RF = friedman.test(k_RF$AUC, k_RF$FS, k_RF$rep)
ft_SVM = friedman.test(k_SVM$AUC, k_SVM$FS, k_SVM$rep)
ft_xGBOOST = friedman.test(k_xGBOOST$AUC, k_xGBOOST$FS, k_xGBOOST$rep)
ft_Counts = list(ft_GBM,ft_GLMNET,ft_RF,ft_SVM,ft_xGBOOST )
##### Fried Test #####



#Relative Abundances 
#######
dfcbf = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds)
dkw = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_KW_12.rds)
dldm = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_LDM_39.rds)

l = list(dfcbf,dkw, dldm )
fs = list("FCBF", "KW", "LDM")
lista= list()

for (i in seq_along(l)){
  lista[[i]]= comp.models(l[[i]], fs[[i]])
}
df = rbind(lista[[1]], lista[[2]], lista[[3]])


bp <- ggplot(df, aes(x=FS, y=AUC, group=FS)) +
  geom_boxplot(aes(fill=FS)) +
  theme_light()+ theme(axis.title.y = element_blank(), axis.ticks.x = element_blank(),legend.position = "bottom",legend.title = element_blank(),
                       axis.text.x = element_blank())+
  facet_wrap(~ Algorithm, scales ="fixed")+ #coord_cartesian(ylim = c(0.5, 1))+
  #stat_compare_means(comparisons = list(c("FCBF", "KW"),c("FCBF", "LDM"), c("KW", "LDM")),method = "wilcox.test", bracket.size = 0.2,  size = 2, paired = FALSE,label = "p.signif",hide.ns = TRUE,vjust = 0.5)+
  theme(strip.background = element_blank(),strip.text.x = element_text(
    size = 8, color = "black"), axis.ticks = element_blank(),axis.title.x.bottom = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank())

bp =change_palette(bp, palette= viridis(3))
bp_RA = bp+ggtitle("R.A Benchmark")+theme(legend.position = c(0.82, 0.25),
                                          legend.text =  element_text(size=10),
                                          legend.key = element_rect(size = 6),
                                          legend.key.height = unit(1, "cm"),
                                          legend.key.width = unit(1, "cm"),
                                          plot.title = element_text(hjust = 0.5))
##### Fried Test #####
k_GBM = df[df$Algorithm == "GBM",]
k_GLMNET = df[df$Algorithm == "GLMNET",]
k_RF = df[df$Algorithm == "RF",]
k_SVM = df[df$Algorithm == "SVM",]
k_xGBOOST = df[df$Algorithm == "xGBOOST",]

ft_GBM = friedman.test(k_GBM$AUC, k_GBM$FS, k_GBM$rep)
ft_GLMNET = friedman.test(k_GLMNET$AUC, k_GLMNET$FS, k_GLMNET$rep)
ft_RF = friedman.test(k_RF$AUC, k_RF$FS, k_RF$rep)
ft_SVM = friedman.test(k_SVM$AUC, k_SVM$FS, k_SVM$rep)
ft_xGBOOST = friedman.test(k_xGBOOST$AUC, k_xGBOOST$FS, k_xGBOOST$rep)
ft_RA = list(ft_GBM,ft_GLMNET,ft_RF,ft_SVM,ft_xGBOOST )
##### Fried Test #####

#######
#COUNTS
GBM.KW =as.data.frame(counts$Bmr_Ravel_Genus_C_train_KW_12.rds$results$dataset$classif.gbm.tuned$measures.test)
GLMNET.FCBF = as.data.frame(counts$Bmr_Ravel_Genus_C_train_FCBF_8.rds$results$dataset$classif.glmnet.tuned$measures.test)
RF.FCBF = as.data.frame(counts$Bmr_Ravel_Genus_C_train_FCBF_8.rds$results$dataset$classif.randomForest.tuned$measures.test)
SVM.FCBF = as.data.frame(counts$Bmr_Ravel_Genus_C_train_FCBF_8.rds$results$dataset$classif.ksvm.tuned$measures.test)
xGBOOST.FCBF = as.data.frame(counts$Bmr_Ravel_Genus_C_train_FCBF_8.rds$results$dataset$classif.xgboost.tuned$measures.test)

l = list(GBM.KW, GLMNET.FCBF, RF.FCBF, SVM.FCBF, xGBOOST.FCBF)
names = c("GBM.KW", "GLMNET.FCBF", "RF.FCBF", "SVM.FCBF", "xGBOOST.FCBF")
names(l)= names

require(plyr)
comp.models2= function(dataset){
  reps = rep(c(rep('1-rep', 5), rep('2-rep', 5), rep('3-rep', 5),rep('4-rep', 5),rep('5-rep', 5),rep('6-rep', 5),rep('7-rep', 5),rep('8-rep', 5),rep('9-rep', 5),rep('10-rep', 5)))
  data = cbind.data.frame(dataset, reps)
  data <- subset( data, select = -iter )
  aggsdata = ddply(data, .(reps), summarize, AUC = mean(auc), ACC = mean(acc), MMCE = mean(mmce))
  return(aggsdata)
}

for (i in seq_along(l)){
  l[[i]] = comp.models2(l[[i]])
}

for (i in seq_along(l)){
  Model = c(rep(names(l[i]),10))
  l[[i]] = cbind(l[[i]], Model)
}

df = rbind.data.frame(l[[1]],l[[2]], l[[3]], l[[4]], l[[5]])

my_comparaisons <- list(c("GBM.KW", "GLMNET.FCBF"),
                        c("GBM.KW", "RF.FCBF"),
                        c("GBM.KW", "SVM.FCBF"),
                        c("GBM.KW","xGBOOST.FCBF"),
                        c("GLMNET.FCBF","RF.FCBF"),
                        c("GLMNET.FCBF","SVM.FCBF"),
                        c("GLMNET.FCBF","xGBOOST.FCBF"),
                        c("RF.FCBF","SVM.FCBF"),
                        c("RF.FCBF","xGBOOST.FCBF"),
                        c("SVM.FCBF","xGBOOST.FCBF"))


mod_C =  ggplot(df, aes(x = Model, y = AUC, color = Model))+
  geom_jitter()+
  geom_boxplot(aes(fill = Model), notch = F)+
  scale_fill_manual(values = viridis(5))+
  theme_light()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text=element_text(size=10), axis.title.y = element_text(size = 10),
        axis.ticks = element_blank(),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10),legend.position = "bottom")+
  ggtitle("Counts Best models")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_compare_means(comparisons = my_comparaisons, bracket.size = 0.3,  size = 3,paired = TRUE,label = "p.signif",hide.ns = TRUE,vjust = 0.5)
mod_C

ft_comparation_C = friedman.test(df$AUC, df$Model, df$reps)

#AR
GBM.FCBF = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds$results$dataset$classif.gbm.tuned$measures.test)
GLMNET.FCBF = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds$results$dataset$classif.glmnet.tuned$measures.test)
RF.FCBF = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds$results$dataset$classif.randomForest.tuned$measures.test)
SVM.FCBF = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds$results$dataset$classif.ksvm.tuned$measures.test)
xGBOOST.FCBF = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds$results$dataset$classif.xgboost.tuned$measures.test)

l = list(GBM.FCBF, GLMNET.FCBF, RF.FCBF, SVM.FCBF, xGBOOST.FCBF)
names = c("GBM.FCBF", "GLMNET.FCBF", "RF.FCBF", "SVM.FCBF", "xGBOOST.FCBF")
names(l)= names


for (i in seq_along(l)){
  l[[i]] = comp.models2(l[[i]])
}

for (i in seq_along(l)){
  Model = c(rep(names(l[i]),10))
  l[[i]] = cbind(l[[i]], Model)
}

df = rbind.data.frame(l[[1]],l[[2]], l[[3]], l[[4]], l[[5]])

my_comparaisons <- list(c("GBM.FCBF", "GLMNET.FCBF"),
                        c("GBM.FCBF", "RF.FCBF"),
                        c("GBM.FCBF", "SVM.FCBF"),
                        c("GBM.FCBF","xGBOOST.FCBF"),
                        c("GLMNET.FCBF","RF.FCBF"),
                        c("GLMNET.FCBF","SVM.FCBF"),
                        c("GLMNET.FCBF","xGBOOST.FCBF"),
                        c("RF.FCBF","SVM.FCBF"),
                        c("RF.FCBF","xGBOOST.FCBF"),
                        c("SVM.FCBF","xGBOOST.FCBF"))
mod_AR =  ggplot(df, aes(x = Model, y = AUC, color = Model))+
  geom_jitter()+
  geom_boxplot(aes(fill = Model), notch = F)+
  scale_fill_manual(values = viridis(5))+
  theme_light()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text=element_text(size=10), axis.title.y = element_text(size = 10),
        axis.ticks = element_blank(),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10),legend.position = "bottom")+
  ggtitle("RA Best models")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_compare_means(comparisons = my_comparaisons, bracket.size = 0.3,  size = 3,paired = TRUE,label = "p.signif",hide.ns = TRUE,vjust = 0.5)
mod_AR

ft_comparation_RA = friedman.test(df$AUC, df$Model, df$reps)
###########
# B) VENNPLOTS COUNTS vs AR (FEATURES)
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
venn_C  = ggVennDiagram(list_C,color =viridis(1),lty = 1, label_alpha = 0.5, label= "count") 
venn_C  = venn_C  + scale_fill_gradient(high="#21908CFF",low = "#FDE725FF")
venn_C  = venn_C  + labs(fill = "Nº Features")+ ggtitle("Counts Features Intersection")+theme(plot.title = element_text(hjust = 0.5))
l2 = get_legend(venn_C)
venn_C = venn_C+theme(legend.position = "none")
venn_C

venn_AR  = ggVennDiagram(list_AR,color =viridis(1),lty = 1, label_alpha = 0.5, label= "count" ) 
venn_AR  = venn_AR  + scale_fill_gradient(high="#21908CFF",low = "#FDE725FF")
venn_AR  = venn_AR  + labs(fill = "Nº Features")+ggtitle("R.A Features Intersection")+ theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
venn_AR

mod_C = mod_C + theme(legend.position = "bottom")
mod_AR = mod_AR + theme(legend.position = "bottom")


panel1 = ggarrange(bp_counts,bp_RA,venn_C,venn_AR,
                   ncol = 2,nrow =2,
                   labels = list("A","","B"))
#######

#### PANEL 2 #####
#A) Feature importance
#Counts
require(mlr)
source("git/Entropy/functions/FunctionsDataFilter.R")
models = as.data.frame(counts$Bmr_Ravel_Genus_C_train_FCBF_8.rds)
model = models[models$learner.id == "classif.randomForest.tuned",]
index = which.max(model$acc)
# Get best model
modelos = getBMRModels(counts$Bmr_Ravel_Genus_C_train_FCBF_8.rds)
alg_index = which(names(modelos$dataset)=="classif.randomForest.tuned")
best_mod_C = getLearnerModel(modelos$dataset[[alg_index]][[index]])
FI_C = getFeatureImportance(best_mod_C)
#View(FI_C$res)
df_FIc = as.data.frame(FI_C$res)
kk = as.vector(df_FIc$variable)
Ravel = readRDS("projects/Entropy/data/Ravel_phyloseq.rds")
Ravel = phy.aglomerate(phyobject = Ravel, rank = "Rank6")
Ravel.subC <- subset_taxa(Ravel, rownames(tax_table(Ravel)) %in% kk)
Ravel.df = as.data.frame(tax_table(Ravel.subC))
Ravel.df = Ravel.df[kk,]
Ravel.df = as.data.frame(cbind(rownames(Ravel.df),Ravel.df$Rank6))
df_FIc$variable = Ravel.df$V2
df_FIc$variable = substr(df_FIc$variable, 4,100)

plotFIC = ggplot(df_FIc, aes(x = reorder(variable, importance), y = importance))+
  geom_segment( aes(xend=variable, yend=0,), color = viridis(3)[2]) +
  geom_point( size=4, color=viridis(1)) +
  coord_flip() +
  theme_light()+
  theme( axis.text=element_text(size=10),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))+
  ggtitle("Feature Importance: Counts+FCBF+RF model")+theme(plot.title = element_text(hjust = 0.5))
plotFIC



#AR
models = as.data.frame(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds)
model = models[models$learner.id == "classif.randomForest.tuned",]
index = which.max(model$acc)

# Get best model
modelos = getBMRModels(ar$Bmr_Ravel_Genus_AR_train_FCBF_9.rds)
alg_index = which(names(modelos$dataset)=="classif.randomForest.tuned")
best_mod_AR = getLearnerModel(modelos$dataset[[alg_index]][[index]])
FI_AR = getFeatureImportance(best_mod_AR)
#View(FI_AR$res)
df_FIAR= as.data.frame(FI_AR$res)
kk2 = as.vector(df_FIAR$variable)
Ravel = readRDS("projects/Entropy/data/Ravel_phyloseq.rds")
Ravel = phy.aglomerate(phyobject = Ravel, rank = "Rank6")
Ravel.subAR <- subset_taxa(Ravel, rownames(tax_table(Ravel)) %in% kk2)
Ravel.df = as.data.frame(tax_table(Ravel.subAR))
Ravel.df = Ravel.df[kk2,]
Ravel.df = as.data.frame(cbind(rownames(Ravel.df),Ravel.df$Rank6))
df_FIAR$variable = Ravel.df$V2
df_FIAR$variable = substr(df_FIAR$variable, 4,100)



plotFIAR = ggplot(df_FIAR, aes(x = reorder(variable, importance), y = importance))+
  geom_segment( aes(xend=variable, yend=0,), color = viridis(3)[2]) +
  geom_point( size=4, color=viridis(1)) +
  coord_flip() +
  theme_light()+
  theme( axis.text=element_text(size=10),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))+
  ggtitle("Feature Importance: RA+FCBF+RF model")+theme(plot.title = element_text(hjust = 0.5))

plotFIAR


#Test Validation
source("git/Entropy/functions/FunctionsGetSplitData.R")
test_C = readRDS("projects/Entropy/data/test/Ravel_Genus_C_test.rds")
test_AR = readRDS("projects/Entropy/data/test/Ravel_Genus_AR_test.rds")
require(ggplotify)
source("git/Entropy/functions/PlotCM.R")
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
library(ggplotify)
CM_C = confusionMatrix(data = prediccion_C$data$response, reference = as.factor(test_C$target),positive = "High")

#TestsPlots
#Plot CM
table(test_task$env$data$target, prediccion_C$data$response)
plotCM_C= as.ggplot(~draw_confusion_matrix(CM_C))
plotCM_C = plotCM_C + ggtitle("Counts Confusion matrix")+theme(plot.title = element_text(hjust = 0.5))
plotCM_C


#Plot AUC_Curve
asRoc = asROCRPrediction(prediccion_C)
p_Counts = ROCR::performance(asRoc, "tpr", "fpr")
plotAUC_C = as.ggplot(~plot(p_Counts))
plotAUC_C= plotAUC_C+ ggtitle("Test validation using Counts")+theme(plot.title = element_text(hjust = 0.5))

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
CM_RA = confusionMatrix(data = prediccionAR$data$response, reference = as.factor(test_AR$target),positive = "High")

#TestsPlots
#Plot CM
plotCM_RA= as.ggplot(~draw_confusion_matrix(CM_RA))
plotCM_RA = plotCM_RA + ggtitle("RA Confusion matrix")+theme(plot.title = element_text(hjust = 0.5))
plotCM_RA

#Plot AUC_Curve
asRoc2 = asROCRPrediction(prediccionAR)
p_RA = ROCR::performance(asRoc2, "tpr", "fpr")
plotAUC_AR = as.ggplot(~plot(p_RA))
plotAUC_AR = plotAUC_AR + ggtitle("Test validation using RA")+theme(plot.title = element_text(hjust = 0.5))

#PANEL
panel2 = ggarrange(plotFIC, plotFIAR,plotCM_C,plotCM_RA,plotAUC_C ,plotAUC_AR,
                   ncol = 2, nrow = 3,widths = c(4,4,1),
                   labels = list("A", "", "B", "", "C"))
panel2



