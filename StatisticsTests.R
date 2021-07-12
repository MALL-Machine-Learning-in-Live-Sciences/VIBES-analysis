#### NUGENT @@@@@
# KRUSKAL WALLIS Y DUNN

library(tidyverse)
library(ggpubr)
library(rstatix)

counts = readRDS(file = "projects/Entropy/data/benchmarks/Nugent/Ravel_Genus_C_Second_Benchmarks.rds")
ar = readRDS(file = "projects/Entropy/data/benchmarks/Nugent/Ravel_Genus_AR_Second_Benchmarks.rds")

# COUNTS
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

#Per Algorithm
GBM = df[df$Algorithm == "GBM",]
GLMNET = df[df$Algorithm == "GLMNET",]
RF = df[df$Algorithm == "RF",]
SVM = df[df$Algorithm == "SVM",]
xGBOOST = df[df$Algorithm == "xGBOOST",]

Al.list = list(GBM,GLMNET,RF,SVM,xGBOOST)
names = c("GBM","GLMNET","RF", "SVM","xGBOOST")
names2 = c("Kruskal", "Dunn", "Wilcox") 
names(Al.list) = names
tests = list()
for (i in seq_along(Al.list)){
  tests[[i]] = list(kruskal_test(AUC ~ FS, data = Al.list[[i]]),
                    dunn_test(AUC ~ FS, p.adjust.method = "bonferroni",data = Al.list[[i]]),
                    wilcox_test(AUC ~ FS, p.adjust.method = "bonferroni", data = Al.list[[i]])
                    )
  names(tests[[i]]) = names2
}
names(tests) = names
tests$xGBOOST


#RA
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

#Per Algorithm
GBM = df[df$Algorithm == "GBM",]
GLMNET = df[df$Algorithm == "GLMNET",]
RF = df[df$Algorithm == "RF",]
SVM = df[df$Algorithm == "SVM",]
xGBOOST = df[df$Algorithm == "xGBOOST",]

Al.list = list(GBM,GLMNET,RF,SVM,xGBOOST)
names = c("GBM","GLMNET","RF", "SVM","xGBOOST")
names2 = c("Kruskal", "Dunn", "Wilcox") 
names(Al.list) = names
tests = list()
for (i in seq_along(Al.list)){
  tests[[i]] = list(kruskal_test(AUC ~ FS, data = Al.list[[i]]),
                    dunn_test(AUC ~ FS, p.adjust.method = "bonferroni",data = Al.list[[i]]),
                    wilcox_test(AUC ~ FS, p.adjust.method = "bonferroni", data = Al.list[[i]])
  )
  names(tests[[i]]) = names2
}
names(tests) = names


##### AMSEL ######
# KRUSKAL WALLIS Y DUNN

library(tidyverse)
library(ggpubr)
library(rstatix)

counts = readRDS(file = "projects/Entropy/data/benchmarks/Amsel/Sriniv_Amsel_Genus_C_Second_Benchmarks.rds")
ar = readRDS(file = "projects/Entropy/data/benchmarks/Amsel/Sriniv_Amsel_Genus_AR_THIRD_Benchmarks.rds")

# COUNTS
dfcbf = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_FCBF_4.rds)
dkw = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_KW_9.rds)
dldm = as.data.frame(counts$Bmr_Sriniv_Amsel_Genus_C_train_LDM_29.rds)

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

#Per Algorithm
GBM = df[df$Algorithm == "GBM",]
GLMNET = df[df$Algorithm == "GLMNET",]
RF = df[df$Algorithm == "RF",]
SVM = df[df$Algorithm == "SVM",]
xGBOOST = df[df$Algorithm == "xGBOOST",]

Al.list = list(GBM,GLMNET,RF,SVM,xGBOOST)
names = c("GBM","GLMNET","RF", "SVM","xGBOOST")
names2 = c("Kruskal", "Dunn", "Wilcox") 
names(Al.list) = names
tests = list()
for (i in seq_along(Al.list)){
  tests[[i]] = list(kruskal_test(AUC ~ FS, data = Al.list[[i]]),
                    dunn_test(AUC ~ FS, p.adjust.method = "bonferroni",data = Al.list[[i]]),
                    wilcox_test(AUC ~ FS, p.adjust.method = "bonferroni", data = Al.list[[i]])
  )
  names(tests[[i]]) = names2
}
names(tests) = names


#RA
dfcbf = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_FCBF_3.rds)
dkw = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_KW_9.rds)
dldm = as.data.frame(ar$Bmr_Sriniv_Amsel_Genus_AR_train_LDM_24.rds)

l = list(dfcbf,dkw, dldm )
fs = list("FCBF", "KW", "LDM")
lista= list()

for (i in seq_along(l)){
  lista[[i]]= comp.models(l[[i]], fs[[i]])
}
df = rbind(lista[[1]], lista[[2]], lista[[3]])

#Per Algorithm
GBM = df[df$Algorithm == "GBM",]
GLMNET = df[df$Algorithm == "GLMNET",]
RF = df[df$Algorithm == "RF",]
SVM = df[df$Algorithm == "SVM",]
xGBOOST = df[df$Algorithm == "xGBOOST",]

Al.list = list(GBM,GLMNET,RF,SVM,xGBOOST)
names = c("GBM","GLMNET","RF", "SVM","xGBOOST")
names2 = c("Kruskal", "Dunn", "Wilcox") 
names(Al.list) = names
tests = list()
for (i in seq_along(Al.list)){
  tests[[i]] = list(kruskal_test(AUC ~ FS, data = Al.list[[i]]),
                    dunn_test(AUC ~ FS, p.adjust.method = "bonferroni",data = Al.list[[i]]),
                    wilcox_test(AUC ~ FS, p.adjust.method = "bonferroni", data = Al.list[[i]])
  )
  names(tests[[i]]) = names2
}
names(tests) = names



