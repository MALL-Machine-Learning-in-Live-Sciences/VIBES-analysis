require(ggplot2)
require(viridis)
library("ggVennDiagram")
require(plotly)
require(ggpubr)
require(grid)
require(ggplotify)

counts = readRDS(file = "projects/Entropy/data/benchmarks/Sriniv_Amsel_Genus_C_Second_Benchmarks.rds")
ar = readRDS(file = "projects/Entropy/data/benchmarks/Sriniv_Amsel_Genus_AR_Second_Benchmarks.rds")

#P1
#-----------
# A) PLOTS COUNTS vs AR
# COUNTS
dclust = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_CLUST_3.rds)
dfcbf = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_FCBF_4.rds)
dkw = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_KW_9.rds)
dldm = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_LDM_29.rds)
datal = list(dclust,dfcbf,dkw,dldm)
nmes = c("CLUST", "FCBF", "KW","LDM")
names(datal) =nmes


for (i in seq_along(datal)){
  FSnames = c(rep(names(datal[i]),5))
  datal[[i]] = data.frame(FS =FSnames,Algorithm = c("RF", "GLMNET","xGBOOST", "SVM", "GBM"),
                          AUC = c(mean(datal[[i]]$auc[1:50]),mean(datal[[i]]$auc[51:100]),
                                  mean(datal[[i]]$auc[101:150]),mean(datal[[i]]$auc[151:200]),
                                  mean(datal[[i]]$auc[201:250])))
}

dataCounts = rbind(datal[[1]], datal[[2]], datal[[3]], datal[[4]])

# AR
dclust2 = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_CLUST_5.rds)
dfcbf2 = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_FCBF_3.rds)
dkw2 = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_KW_9.rds)
dldm2 = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_LDM_24.rds)
datal2 = list(dclust2,dfcbf2,dkw2,dldm2)
nmes = c("CLUST", "FCBF", "KW","LDM")
names(datal2) =nmes

for (i in seq_along(datal2)){
  FSnames = c(rep(names(datal2[i]),5))
  datal2[[i]] = data.frame(FS =FSnames,Algorithm = c("RF", "GLMNET","xGBOOST", "SVM", "GBM"),
                           AUC = c(mean(datal2[[i]]$auc[1:50]),mean(datal2[[i]]$auc[51:100]),
                                   mean(datal2[[i]]$auc[101:150]),mean(datal2[[i]]$auc[151:200]),
                                   mean(datal2[[i]]$auc[201:250])))
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

#--------------
# C) Comparacion de los mejores modelos de cada tipo de datos
require(plyr)
comp.models= function(dataset){
  reps = rep(c(rep('1-rep', 5), rep('2-rep', 5), rep('3-rep', 5),rep('4-rep', 5),rep('5-rep', 5),rep('6-rep', 5),rep('7-rep', 5),rep('8-rep', 5),rep('9-rep', 5),rep('10-rep', 5)))
  data = cbind.data.frame(dataset, reps)
  data <- subset( data, select = -iter )
  aggsdata = ddply(data, .(reps), summarize, AUC = mean(auc), ACC = mean(acc), MMCE = mean(mmce))
  return(aggsdata)
}

CLUST.GLMNET =as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_CLUST_3.rds$results$dataset$classif.glmnet.tuned$measures.test)
FCBF.GLMNET = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_FCBF_4.rds$results$dataset$classif.glmnet.tuned$measures.test)
KW.RF = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_KW_9.rds$results$dataset$classif.randomForest.tuned$measures.test) 
LDM.RF = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_LDM_29.rds$results$dataset$classif.randomForest.tuned$measures.test)

l = list(CLUST.GLMNET, FCBF.GLMNET, KW.RF, LDM.RF)
names = c("CLUST.GLMNET", "FCBF.GLMNET", "KW.RF", "LDM.RF")
names(l)= names
for (i in seq_along(l)){
  l[[i]] = comp.models(l[[i]])
}
for (i in seq_along(l)){
  Model = c(rep(names(l[i]),10))
  l[[i]] = cbind(l[[i]], Model)
}

df = rbind.data.frame(l[[1]],l[[2]], l[[3]], l[[4]])

my_comparaisons <- list(c("CLUST.GLMNET", "FCBF.GLMNET"),
                        c("CLUST.GLMNET", "KW.RF"),
                        c("CLUST.GLMNET", "LDM.RF"),
                        c("FCBF.GLMNET","KW.RF"),
                        c("FCBF.GLMNET","LDM.RF"),
                        c("KW.RF","LDM.RF"))

mod_C =  ggplot(df, aes(x = Model, y = AUC, color = Model))+
  geom_jitter()+
  geom_boxplot(aes(fill = Model), notch = F)+
  scale_fill_manual(values = viridis(4))+
  theme_light()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text=element_text(size=10), axis.title.y = element_text(size = 10),
        axis.ticks = element_blank(),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10),
        legend.position = "none")+
  ggtitle("Best models Counts")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_compare_means(comparisons = my_comparaisons, bracket.size = 0.3,  size = 3)
mod_C

# Abundancias relativas
CLUST.GLMNET =as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_CLUST_5.rds$results$dataset$classif.glmnet.tuned$measures.test)
FCBF.GLMNET = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_FCBF_3.rds$results$dataset$classif.glmnet.tuned$measures.test)
KW.RF = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_KW_9.rds$results$dataset$classif.randomForest.tuned$measures.test)
LDM.RF = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_LDM_24.rds$results$dataset$classif.randomForest.tuned$measures.test)

l = list(CLUST.GLMNET, FCBF.GLMNET, KW.RF, LDM.RF)
names = c("CLUST.GLMNET", "FCBF.GLMNET", "KW.RF", "LDM.RF")
names(l)= names
for (i in seq_along(l)){
  l[[i]] = comp.models(l[[i]])
}
for (i in seq_along(l)){
  Model = c(rep(names(l[i]),10))
  l[[i]] = cbind(l[[i]], Model)
}

df = rbind.data.frame(l[[1]],l[[2]], l[[3]], l[[4]])

my_comparaisons <- list(c("CLUST.GLMNET", "FCBF.GLMNET"),
                        c("CLUST.GLMNET", "KW.RF"),
                        c("CLUST.GLMNET", "LDM.RF"),
                        c("FCBF.GLMNET","KW.RF"),
                        c("FCBF.GLMNET","LDM.RF"),
                        c("KW.RF","LDM.RF"))

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

#--------------
#A) Feature importance
#Counts
phy.aglomerate = function(phyobject, rank){
  require(phyloseq)
  phy =  tax_glom(phyobject, rank)
  return(phy)
}
require(mlr)
models = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_LDM_29.rds)
model = models[models$learner.id == "classif.randomForest.tuned",]
index = which.max(model$acc)
# Get best model
modelos = getBMRModels(counts$Bmr_Sriniv_Amsel_Genus_C_train_LDM_29.rds)
alg_index = which(names(modelos$dataset)=="classif.randomForest.tuned")
best_mod_C = getLearnerModel(modelos$dataset[[alg_index]][[index]])
FI_C = getFeatureImportance(best_mod_C)
View(FI_C$res)
df_FIc = as.data.frame(FI_C$res)
kk = as.vector(df_FIc$variable)
Sriniv = readRDS("projects/Entropy/data/Sriniv_Amsel_phyloseq.rds")
Sriniv = phy.aglomerate(phyobject = Sriniv, rank = "Rank6")
taxa_names(Sriniv) <- make.names(taxa_names(Sriniv), unique=TRUE)

Sriniv.subC <- subset_taxa(Sriniv, rownames(tax_table(Sriniv)) %in% kk)
Sriniv.df = as.data.frame(tax_table(Sriniv.subC))
Sriniv.df = Sriniv.df[kk,]
Sriniv.df = as.data.frame(cbind(rownames(Sriniv.df),Sriniv.df$Rank6))
df_FIc$variable =Sriniv.df$V2
df_FIc$variable = substr(df_FIc$variable, 4,100)

plotFIC = ggplot(df_FIc, aes(x = reorder(variable, importance), y = importance))+
  geom_segment( aes(xend=variable, yend=0,), color = viridis(3)[2]) +
  geom_point( size=4, color=viridis(1)) +
  coord_flip() +
  theme_light()+
  theme( axis.text=element_text(size=10),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))
plotFIC





#AR
models = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_KW_9.rds)
model = models[models$learner.id == "classif.randomForest.tuned",]
index = which.max(model$auc)

# Get best model
modelos = getBMRModels(ar$Bmr_Sriniv_Amsel_Genus_AR_train_KW_9.rds)
alg_index = which(names(modelos$dataset)=="classif.randomForest.tuned")
best_mod_AR = getLearnerModel(modelos$dataset[[alg_index]][[index]])
FI_AR = getFeatureImportance(best_mod_AR)
View(FI_AR$res)
df_FIAR= as.data.frame(FI_AR$res)
kk = as.vector(df_FIAR$variable)
Sriniv = readRDS("projects/Entropy/data/Sriniv_Amsel_phyloseq.rds")
Sriniv = phy.aglomerate(phyobject = Sriniv, rank = "Rank6")
taxa_names(Sriniv) <- make.names(taxa_names(Sriniv), unique=TRUE)

Sriniv.subC <- subset_taxa(Sriniv, rownames(tax_table(Sriniv)) %in% kk)
Sriniv.df = as.data.frame(tax_table(Sriniv.subC))
Sriniv.df = Sriniv.df[kk,]
Sriniv.df = as.data.frame(cbind(rownames(Sriniv.df),Sriniv.df$Rank6))
df_FIAR$variable =Sriniv.df$V2
df_FIAR$variable = substr(df_FIAR$variable, 4,100)


plotFIAR = ggplot(df_FIAR, aes(x = reorder(variable, importance), y = importance))+
  geom_segment( aes(xend=variable, yend=0,), color = viridis(3)[2]) +
  geom_point( size=4, color=viridis(1)) +
  coord_flip() +
  theme_light()+
  theme( axis.text=element_text(size=10),axis.title.y = element_blank(),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=10), 
         legend.text=element_text(size=10))
plotFIAR


#PANEL
panel1 = ggarrange(comp_C, comp_AR,mod_C,mod_AR, plotFIC, plotFIAR,
                   ncol = 2, nrow = 3,widths = c(4,4,1),
                   labels = list("A", "", "B", "", "C"))
panel1

#--------
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
library(ggplotify)
confusionMatrix(data = prediccion_C$data$response, reference = as.factor(test_C$target),positive = "pos")
a = asROCRPrediction(prediccion_C)
p = ROCR::performance(a, "tpr", "fpr")
plotAUC_C = as.ggplot(~plot(p))
plotAUC_C
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
confusionMatrix(data = prediccionAR$data$response, reference = as.factor(test_AR$target),positive = "pos")
a2 = asROCRPrediction(prediccionAR)
p2 = ROCR::performance(a2, "tpr", "fpr")
plotAUC_AR = as.ggplot(~plot(p2))
plotAUC_AR










