# Validation Models
# Enfrentar los datos del segundo articulo al modelo sacado con los del primero
# Load phyloseq
require(mlr)
require(phyloseq)
source("git/Entropy/functions/FunctionsDataFilter.R")
source("git/Entropy/functions/FunctionsValidation.R")
source("git/Entropy/functions/FunctionsGetSplitData.R")
source("git/Entropy/functions/FunctionsML.R")
# Load phyloseqs
Ravel = readRDS("projects/Entropy/data/Ravel_phyloseq.rds")
Ravel = phy.aglomerate(phyobject = Ravel, rank = "Rank6")
Sriniv = readRDS("projects/Entropy/data/Sriniv_Nugent_phyloseq.rds")
Sriniv = relat.abun(Sriniv)
Sriniv = phy.aglomerate(phyobject = Sriniv, rank = "Rank6")
#Load Bncmarks
#BMR.C = readRDS(file = "projects/Entropy/data/benchmarks/Nugent/Ravel_Genus_C_Second_Benchmarks.rds")
BMR.AR = readRDS(file = "projects/Entropy/data/benchmarks/Nugent/Ravel_Genus_AR_Second_Benchmarks.rds")


features = get.features(bench = BMR.AR, algoritmo = "classif.xgboost.tuned", indice = 4)

Ravel.df = as.data.frame(tax_table(Ravel))
Ravel.df = Ravel.df[features,]
features = as.data.frame(cbind(rownames(Ravel.df),Ravel.df$Rank6))
Sriniv_sub <- subset_taxa(Sriniv, Rank6 %in% features[,2]) 
Sriniv.df = as.data.frame(tax_table(Sriniv_sub))
features = check.features(Sriniv.df= Sriniv.df, features = features)

#We have to train again the model without this feature
data = readRDS("projects/Entropy/data/train/Ravel_Genus_AR_train_LDM_39.rds")
data <- data[ , !(names(data) %in% features)]
bmr = ML.exec(dataset = data)
best = get.best.model(bncmark = bmr)
#saveRDS(object = best, file = "projects/Entropy/data/validation/CM_Validation/Counts/Counts_Retrain_GBOOST_KW.rds")
#best = readRDS("projects/Entropy/data/validation/CM_Validation/Counts/Counts_Retrain_GBOOST_KW.rds")
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

#Validation xGBOOST
v = best$features
my_data = test_task$env$data[1:28]
my_data2 <- my_data[, v]
my_data2
prediccion <- predict(best, newdata =  my_data2, type = "prob")

#prediccion<- predict(best, task = test_task, type = "prob" )

#saveRDS(object = prediccion, file = "projects/Entropy/data/validation/PredictionAR_RF_FCBF.rds")
library(caret)
CM = confusionMatrix(data = prediccion$data$response, reference = as.factor(test.Sriniv$target))
saveRDS(CM, "projects/Entropy/data/validation/CM_Validation/RA/CM_xGBOOST_LDM.rds")


#BoxPlots
# Counts
l = list.files("projects/Entropy/data/validation/CM_Validation/Counts/")

ldata = list()
for (i in 1:length(l)) {
  ldata[[i]] = readRDS(paste("projects/Entropy/data/validation/CM_Validation/Counts/", l[[i]], sep=''))
}
names(ldata) = l

Algoritmo = rep(c(rep('GBM', 3), rep('GLMNET', 3), rep('RF', 3),rep('SVM', 3),rep('xGBOOST', 3)))
FS = rep(c("FCBF", "KW", "LDM"),5)
ACC = list()
low= list()
upp = list()

for (i in 1:length(l)){
  ACC[[i]] = ldata[[i]]$overall[[1]]
  low[[i]] = ldata[[i]]$overall[[3]]
  upp[[i]] =ldata[[i]]$overall[[4]]
}
ACC = unlist(ACC)
low = unlist(low)
upp = unlist(upp)
data  <- data.frame (Algoritmo  = Algoritmo,
                    FS = FS,
                    ACC = ACC,
                    low = low,
                    upp = upp)

require(ggplot2)
require(viridis)
require(ggpubr)

p = ggplot(data, aes(x=Algoritmo, y=ACC, fill=FS)) +
  geom_bar(position = 'dodge', stat='identity', width=.5) +
  xlab(NULL) + coord_cartesian(ylim=c(0.4,1)) +
  geom_errorbar(aes(ymin = low, ymax = upp),
                width = 0.2,
                position = position_dodge(width = 0.5),
                color="black", size=1)
p = change_palette(p, palette= viridis(3))
p_counts = p_counts+ ggtitle("Counts External Cohort Validation")+theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
p_counts 



#BoxPlots
l = list.files("projects/Entropy/data/validation/CM_Validation/RA/")

ldata = list()
for (i in 1:length(l)) {
  ldata[[i]] = readRDS(paste("projects/Entropy/data/validation/CM_Validation/RA//", l[[i]], sep=''))
}
names(ldata) = l

Algoritmo = rep(c(rep('GBM', 3), rep('GLMNET', 3), rep('RF', 3),rep('SVM', 3),rep('xGBOOST', 3)))
FS = rep(c("FCBF", "KW", "LDM"),5)
ACC = list()
low= list()
upp = list()

for (i in 1:length(l)){
  ACC[[i]] = ldata[[i]]$overall[[1]]
  low[[i]] = ldata[[i]]$overall[[3]]
  upp[[i]] =ldata[[i]]$overall[[4]]
}
ACC = unlist(ACC)
low = unlist(low)
upp = unlist(upp)
data  <- data.frame (Algoritmo  = Algoritmo,
                     FS = FS,
                     ACC = ACC,
                     low = low,
                     upp = upp)


p = ggplot(data, aes(x=Algoritmo, y=ACC, fill=FS)) +
  geom_bar(position = 'dodge', stat='identity', width=.5) +
  xlab(NULL) + coord_cartesian(ylim=c(0.4,1)) +
  geom_errorbar(aes(ymin = low, ymax = upp),
                width = 0.2,
                position = position_dodge(width = 0.5),
                color="black", size=1)
p = change_palette(p, palette= viridis(3))
p_AR = p + ggtitle("RA External Cohort Validation")+theme(axis.title.y =  element_blank(),plot.title = element_text(hjust = 0.5))
leg <- get_legend(p_AR)
p_AR = p_AR + theme(legend.position = "none")
# Convert to a ggplot and print
leg <-as_ggplot(leg)


panel_Val = ggarrange(p_counts, p_AR,leg,
                      ncol = 3, nrow = 1,widths =c(3,3,1))
