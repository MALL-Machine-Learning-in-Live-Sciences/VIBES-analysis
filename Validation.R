# Validation Models
# Enfrentar los datos del segundo articulo al modelo sacado con los del primero
# Load phyloseq
require(mlr)
require(phyloseq)
source("git/Entropy/functions/FunctionsDataFilter.R")
source("git/Entropy/functions/FunctionsGetSplitData.R")
Ravel = readRDS("projects/Entropy/data/Ravel_phyloseq.rds")
Ravel = phy.aglomerate(phyobject = Ravel, rank = "Rank6")
Sriniv = readRDS("projects/Entropy/data/Sriniv_Nugent_phyloseq.rds")
Sriniv = phy.aglomerate(phyobject = Sriniv, rank = "Rank6")

#------
#Load model
modelsC = readRDS("projects/Entropy/data/benchmarks/Ravel_Genus_C_Benchmarks.rds")
modelsAR = readRDS("projects/Entropy/data/benchmarks/Ravel_Genus_AR_Benchmarks.rds")

#Get algoruthm and index 
model = as.data.frame(modelsC$Bmr_Ravel_Genus_C_train_FCBF_8.rds)
model = model[model$learner.id == "classif.ksvm.tuned",]
index = which.max(model$auc)

# Get best model
modelos = getBMRModels(modelsC$Bmr_Ravel_Genus_C_train_FCBF_8.rds)
modelos = modelos$dataset$classif.ksvm.tuned
best = getLearnerModel(modelos$dataset$classif.ksvm.tuned[[index]])
#Get features
names = best$features
#-------

# Extract corespondence between OTU anf Genus
Ravel = as.data.frame(tax_table(Ravel))
features = Ravel[names,]
features = cbind(rownames(features),features$Rank6)
vg = features[,2]
vn = features[,1]
#Extract a subset from second phyloseq with only features from the model
Sriniv_sub <- subset_taxa(Sriniv, Rank6 %in% vg) # We have 7 from 8 features
tax_table(Sriniv_sub)

#We have to train again the model without this feature
FCBF = readRDS("projects/Entropy/data/train/Ravel_Genus_C_train_FCBF_8.rds")
FCBF <- subset(FCBF, select = -NR_028773.1 )
bmFCBF = ML.exec(dataset = FCBF)
modelo = getBMRModels(bmFCBF)
resFCBF = as.data.frame(bmFCBF)

# Check the model with best results(we asume here that we have only 1 algorithm)
# We can change this according to algorithms used
models_index = resFCBF[
  order( resFCBF[,5], resFCBF[,4],  resFCBF[,6] ),
]
index = as.numeric(tail(rownames(models_index), 1))
# Get the model
best = getLearnerModel(modelo$dataset$classif.ksvm.tuned[[index]])
# Get the featureas that we ll have to extract from test
features = best$features
#
Ravel = as.data.frame(tax_table(Ravel))
features = Ravel[features,]
features = cbind(rownames(features),features$Rank6)
vg = features[,2]
vn = features[,1]

# Extract a subset from second phyloseq with only features from the model
Sriniv_sub <- subset_taxa(Sriniv, Rank6 %in% vg) 
table = as.data.frame(tax_table(Sriniv_sub))
table = cbind(rownames(table), table$Rank6)
k = table[,2]
Sriniv_sub = get.dataset(phyobject = Sriniv_sub)
cols <- sapply(Sriniv_sub, is.numeric)
variables = Sriniv_sub[cols]
targets = Sriniv_sub$target
colnames(variables) = k
dataset = as.data.frame(variables)
dataset = as.data.frame(t(dataset))
features = as.data.frame(features)
require(dplyr)
df <- tibble::rownames_to_column(dataset, "V2")
aa = merge(x = df, y = features, by = "V2", all = TRUE)
dataset = as.data.frame(t(dataset))
dataset = as.data.frame(cbind(dataset, target =targets ))

#
library(tidyverse)
aa = aa %>% remove_rownames %>% column_to_rownames(var="V1")
test.Sriniv <- subset( aa, select = -V2)
test.Sriniv = t(test.Sriniv)
test.Sriniv = as.data.frame(cbind(test.Sriniv, target =targets ))

# Type numeric
test.Sriniv$NR_044929.2 = as.numeric(test.Sriniv$NR_044929.2)
test.Sriniv$NR_113356.1 = as.numeric(test.Sriniv$NR_113356.1)
test.Sriniv$NR_117757.1 = as.numeric(test.Sriniv$NR_117757.1)
test.Sriniv$NR_118377.1 = as.numeric(test.Sriniv$NR_118377.1)
test.Sriniv$NR_036982.1 = as.numeric(test.Sriniv$NR_036982.1)
test.Sriniv$NR_041796.1 = as.numeric(test.Sriniv$NR_041796.1)
test.Sriniv$NR_113093.1 = as.numeric(test.Sriniv$NR_113093.1)

test.Sriniv = norm.dataset(test.Sriniv)

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

