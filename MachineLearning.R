#ML CESGA
#setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data')
#path = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/train/"
source("git/Entropy/functions/FunctionsML.R")
path = "projects/Entropy/data/train/"
pattern = "Ravel_Genus_C"
l = list.files(path, pattern = pattern)

ldata = list()
for (i in 1:length(l)) {
  ldata[[i]] = readRDS(paste(path, l[[i]], sep=''))
}
names(ldata) = l
bmrs = list()
for (i in seq_along(ldata)){
  bmrs[[i]] = ML.exec(dataset = ldata[[i]])
}
#Rename Benchmarks by data input
namesbmr <- paste("Bmr_",l, sep="")
names(bmrs) = namesbmr
#Get models from dataset desirable
modelos_FCBF = getBMRModels(bmrs$Bmr_Ravel_Genus_C_train_FCBF_8.rds)
resFCBF = as.data.frame(bmrs$Bmr_Ravel_Genus_C_train_FCBF_8.rds)

# Check the model with best results(we asume here that we have only 1 algorithm)
# We can change this according to algorithms used
which.max(resFCBF$auc); which.max(resFCBF$acc); which.min(resFCBF$mmce)
models_index = resFCBF[
  order( resFCBF[,5], resFCBF[,4],  resFCBF[,6] ),
]
index = as.numeric(tail(rownames(models_index), 1))
# Get the model
best = getLearnerModel(modelos_FCBF$dataset$classif.randomForest.tuned[[index]])
# Get the featureas that we ll have to extract from test
features = best$features
#pathw = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/test//Ravel_Genus_C_test.rds"
path2 = "projects/Entropy/data/test/Ravel_Genus_C_test.rds"
test = readRDS(path2)
target = test$target
test = test[features]
test = cbind(test, target)
test_task = makeClassifTask(data = test, target = "target")
test_task = normalizeFeatures(
  test_task,
  method = "range",
  cols = NULL,
  range = c(0, 1),
  on.constant = "quiet")
prediccion = predict(best, task= test_task)
saveRDS(object = bmrs , file = "projects/Entropy/data/benchmarks/Ravel_Genus_C_Benchmarks.rds")
saveRDS(object = best, file = "projects/Entropy/data/models/RF_8.rds")
