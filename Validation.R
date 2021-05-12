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
#BMR.C = readRDS("projects/Entropy/data/benchmarks/Ravel_Genus_C_Second_Benchmarks.rds")
BMR.AR = readRDS("projects/Entropy/data/benchmarks/Ravel_Genus_AR_Second_Benchmarks.rds")


features = get.features(bench = BMR.AR, algoritmo = "classif.randomForest.tuned", indice = 2)

Ravel.df = as.data.frame(tax_table(Ravel))
Ravel.df = Ravel.df[features,]
features = as.data.frame(cbind(rownames(Ravel.df),Ravel.df$Rank6))
Sriniv_sub <- subset_taxa(Sriniv, Rank6 %in% features[,2]) 
Sriniv.df = as.data.frame(tax_table(Sriniv_sub))
features = check.features(Sriniv.df= Sriniv.df, features = features)

#We have to train again the model without this feature
data = readRDS("projects/Entropy/data/train/Ravel_Genus_AR_train_FCBF_8.rds")
data <- data[ , !(names(data) %in% features)]
bmr = ML.exec(dataset = data)
best = get.best.model(bncmark = bmr)
saveRDS(object = best, file = "projects/Entropy/data/models/RetrainBest_AR_RF_FCBF_7Feat.rds")
best = readRDS("projects/Entropy/data/models/RetrainBest_C_RF_FCBF8Feat.rds")
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
saveRDS(object = prediccion, file = "projects/Entropy/data/validation/PredictionAR_RF_FCBF.rds")
library(caret)
table(test.Sriniv$target,prediccion$data$response)
confusionMatrix(data = prediccion$data$response, reference = as.factor(test.Sriniv$target))
a = asROCRPrediction(prediccion)
p = ROCR::performance(a, "tpr", "fpr")
plot(p)

