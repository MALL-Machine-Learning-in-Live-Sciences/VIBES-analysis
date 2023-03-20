require(phyloseq)
setwd("~/git/BV_Microbiome/")
rank <-  "Species" # Genus or "Species"
nfeat <- "14" # 24 and 14 for Species or 37 and 21 for Genus

Ravel <- readRDS(paste0("extdata/", rank, "Intersect/Ravel_",
                        rank, "_pseq_", nfeat, ".rds"))
Sriniv <- readRDS(paste0("extdata/", rank, "Intersect/Sriniv_",
                        rank, "_pseq_", nfeat, ".rds"))
PR3020   <- readRDS(paste0("extdata/", rank, "Intersect/PRJNA3020_",
                        rank, "_pseq_", nfeat, ".rds"))
# Select Target
select.target <- function(pseq, tget) {
  require(dplyr)
  df <- sample_data(pseq)
  colnames(df)[colnames(df) == tget] <- "target"
  sample_data(pseq) <- df
  return(pseq)
}
norm.dataset <- function(pseq) {
  # Change columns by rows too, interested in maintain fts in columns
  require(phyloseq)
  otu <- data.frame(otu_table(pseq))
  # Normalize that variables
  otu <- apply(X = otu, FUN = function(x) log2(x + 1), MARGIN = 2)
  otu_table(pseq) <- otu_table(otu, taxa_are_rows = FALSE)
  return(pseq)
}

#3C
Ravel_3C = select.target(pseq = Ravel, tget ="Nugent_score_category")
Sriniv_3C = select.target(pseq = Sriniv, tget = "nugent")

#Norm Datasets
Ravel_3C = norm.dataset(Ravel_3C)
Sriniv_3C = norm.dataset(Sriniv_3C)
PR3020 = norm.dataset(PR3020)

#Clustering functions
make.cluster = function(pseq, n_clusters){
  require(factoextra)
  mat = as.matrix(data.frame(otu_table(pseq)))
  meta = sample_data(pseq)
  set.seed(1312)
  km.res <- kmeans(mat, n_clusters, iter.max = 100, nstart = 25)
  a = fviz_cluster(km.res,mat, geom = "point")
  identical(names(km.res$cluster), rownames(sample_data(pseq)))
  sample_data(pseq)$cluster <- km.res$cluster
  sample_data(pseq)$cluster <- sub("^", "C", sample_data(pseq)$cluster)
  table(sample_data(pseq)$target , sample_data(pseq)$cluster)
  lista = list(pseq, km.res, a)
  names = c("pseq", "km.res", "plot")
  names(lista) = names
  return(lista)
}
make.cl_predict= function(lista, pseq2){
  require(clue)
  pseq = lista$pseq
  km.res = lista$km.res
  mat = as.matrix(data.frame(otu_table(pseq)))
  meta = sample_data(pseq)
  ft = colnames(mat)
  mat2 = as.matrix(data.frame(otu_table(pseq2)))
  new_target = cl_predict(
    object = km.res,
    newdata = mat2[,ft])
  sample_data(pseq2)$cluster <- new_target
  sample_data(pseq2)$cluster <- sub("^", "C", sample_data(pseq2)$cluster)
  return(pseq2)
}

Clust_3C = make.cluster(pseq = Ravel_3C, n_clusters = 3)
Ravel = Clust_3C$pseq
table(sample_data(Ravel)$target , sample_data(Ravel)$cluster)
Sriniv = make.cl_predict(lista = Clust_3C, pseq2 = Sriniv_3C)
table(sample_data(Sriniv)$target , sample_data(Sriniv)$cluster)
PR3020 = make.cl_predict(lista = Clust_3C, pseq2 = PR3020)
table(sample_data(PR3020)$cluster)

# STAND BY
if (rank == "Genus") {
  # Rename clusters with biological meannig
  # Ravel
  sample_data(Ravel)$cluster[sample_data(Ravel)$cluster == "C1"] <- "H"
  sample_data(Ravel)$cluster[sample_data(Ravel)$cluster == "C2"] <- "D"
  sample_data(Ravel)$cluster[sample_data(Ravel)$cluster == "C3"] <- "ID"
  # Sriniv
  sample_data(Sriniv)$cluster[sample_data(Sriniv)$cluster == "C1"] <- "H"
  sample_data(Sriniv)$cluster[sample_data(Sriniv)$cluster == "C2"] <- "D"
  sample_data(Sriniv)$cluster[sample_data(Sriniv)$cluster == "C3"] <- "ID"
  
}else if (rank == "Species"){
  # Rename clusters with biological meanig
  # Ravel
  sample_data(Ravel)$cluster[sample_data(Ravel)$cluster == "C3"] <- "H"
  sample_data(Ravel)$cluster[sample_data(Ravel)$cluster == "C2"] <- "D"
  sample_data(Ravel)$cluster[sample_data(Ravel)$cluster == "C1"] <- "ID"
  # Sriniv
  sample_data(Sriniv)$cluster[sample_data(Sriniv)$cluster == "C3"] <- "H"
  sample_data(Sriniv)$cluster[sample_data(Sriniv)$cluster == "C2"] <- "D"
  sample_data(Sriniv)$cluster[sample_data(Sriniv)$cluster == "C1"] <- "ID"
}else{
  print("Introduce valid taxon level (Genus or Species)")
}

# Saving
saveRDS(object = Ravel, file = paste0("extdata/", rank, "Intersect/Cluster/Ravel_Cluster_", rank, "_pseq.rds"))
saveRDS(object = Sriniv, file = paste0("extdata/", rank, "Intersect/Cluster/Sriniv_Cluster_", rank, "_pseq.rds"))

#Train
ML.exec_C3 = function(dataset){
require(mlr)
require(methods)
require(parallel)
require(parallelMap)
drops <- c("target")
dataset = dataset[ , !(names(dataset) %in% drops)]
cores = detectCores()
task = makeClassifTask(data = dataset, target = 'cluster')
task = normalizeFeatures(
  task,
  method = "range",
  cols = NULL,
  range = c(0, 1),
  on.constant = "quiet")

# Hyperparameter tuning
ctrl<-makeTuneControlGrid()
inner<-makeResampleDesc("Holdout")

# Random Forest
psRF<-makeParamSet(
  makeDiscreteParam("mtry", values = c(round(sqrt(ncol(task$env$data)))-1,
                                       round(sqrt(ncol(task$env$data))),
                                       round(sqrt(ncol(task$env$data)))+1,
                                       round(sqrt(ncol(task$env$data)))+2,
                                       round(sqrt(ncol(task$env$data)))+3,
                                       round(sqrt(ncol(task$env$data)))+4)),
  makeDiscreteParam("ntree", values= 1000L),
  makeDiscreteParam("nodesize", values= c(1:10))
)
l1<-makeLearner("classif.randomForest", predict.type = "prob")
lrn_RF<-makeTuneWrapper(l1,  resampling = inner, par.set = psRF, measures = acc, control=ctrl,  show.info = T)

# GLMNET
psGL = makeParamSet(
  makeDiscreteParam("lambda", c(0.0001,0.001,0.01,0.1,1)),
  makeDiscreteParam("alpha",c(0,0.15,0.25,0.35,0.5,0.65,0.75,0.85,1))
)
l2<-makeLearner("classif.glmnet", predict.type = "prob")
lrn_glmnet<-makeTuneWrapper(l2, inner, psGL, measures = auc, ctrl, show.info=T)

# xGboost
psGB = makeParamSet(
  makeDiscreteLearnerParam("booster", values = c("gbtree", "gblinear", "dart")),
  makeNumericParam("eta", lower = 0, upper = 1),
  makeNumericLearnerParam("lambda", upper = 1, lower = 0),
  makeIntegerParam("max_depth", lower = 1, upper = 20),
  makeDiscreteParam("eval_metric", "logloss"))
l3 = makeLearner("classif.xgboost", predict.type = "prob", nrounds=10)
lrn_GB = makeTuneWrapper(learner = l3, resampling = inner, measures = auc, par.set = psGB, control = ctrl, show.info = T)

# SVM
psKSVM = makeParamSet(makeDiscreteParam('C', values = 2^c(-8, -6, -4, -2, 0, 2, 4, 6, 8)),
                      makeDiscreteParam('sigma', values = 2^c(-8, -6, -4, -2, 0,2 , 4, 6, 8)))
l4 = makeLearner("classif.ksvm", predict.type = "prob")
lrn_KSVM = makeTuneWrapper(learner = l4, resampling = inner, measures = auc, par.set = psKSVM, control = ctrl, show.info = T)

# GBM
psGBM = makeParamSet(makeDiscreteParam("distribution", values = c("bernoulli", "gaussian", "huberized")),
                     makeIntegerParam("n.trees", lower = 100, upper = 800), 
                     makeIntegerParam("interaction.depth", lower = 2, upper = 10),
                     makeNumericParam("bag.fraction", lower = 0.90, upper = 0.90))
l5 = makeLearner("classif.gbm", predict.type = "prob")
lrn_GBM =  makeTuneWrapper(learner = l5, resampling = inner, measures = auc, par.set = psGBM, control = ctrl, show.info = T)

#learners = list(lrn_RF, lrn_glmnet, lrn_GB, lrn_KSVM, lrn_GBM)
learners = (lrn_RF)
# Outer
outer = makeResampleDesc('RepCV' , reps = 10, folds = 5 , stratify = T)

# Benchmarking
#parallelStartSocket(cpus = detectCores()*0.5, level = 'mlr.tuneParams')
parallelStartMulticore(cores , level = 'mlr.tuneParams')
bmr = benchmark(learners, task, outer, measures =  list(acc,mmce), show.info = T, models = T)
parallelStop()
return(bmr)
}
set.seed(1312)
bmr_2C = ML.exec_C2(dataset = Clust_2C$dataframe)
bmr_3C = ML.exec_C3(dataset = Clust_3C$dataframe)

#Test
#Predict
require(caret)
get.best.model_C2 = function(bncmark){
  modelo = getBMRModels(bncmark)
  modelodf = as.data.frame(bncmark)
  
  # Check the model with best results(we asume here that we have only 1 algorithm)
  # We can change this according to algorithms used
  models_index = modelodf[
    order(modelodf[,4], modelodf[,5], modelodf[,6]),
  ]
  index = as.numeric(tail(rownames(models_index), 1))
  # Get the model
  best = getLearnerModel(modelo$dataset[[1]][[index]])
  return(best)
}
get.best.model_C3 = function(bncmark){
  modelo = getBMRModels(bncmark)
  modelodf = as.data.frame(bncmark)
  
  # Check the model with best results(we asume here that we have only 1 algorithm)
  # We can change this according to algorithms used
  models_index = modelodf[
    order( modelodf[,4], modelodf[,5]),
  ]
  index = as.numeric(tail(rownames(models_index), 1))
  # Get the model
  best = getLearnerModel(modelo$dataset[[1]][[index]])
  return(best)
}
bm2 =get.best.model_C2(bncmark = bmr_2C)
test_task2 = makeClassifTask(data = Sriniv.predict_2C, target = "cluster")
test_task2 = normalizeFeatures(
  test_task2,
  method = "range",
  cols = NULL,
  range = c(0, 1),
  on.constant = "quiet")
predict2 = predict(bm2, task = test_task2, type = "prob")
CM2 = confusionMatrix(data = predict2$data$response, reference = as.factor(Sriniv.predict_2C$cluster),positive = "C2")

bm3 = get.best.model_C3(bncmark = bmr_3C)
test_task3 = makeClassifTask(data = Sriniv.predict_3C, target = "cluster")
test_task3 = normalizeFeatures(
  test_task3,
  method = "range",
  cols = NULL,
  range = c(0, 1),
  on.constant = "quiet")
predict3 = predict(bm3, task = test_task3, type = "prob")
CM3 = confusionMatrix(data = predict3$data$response, reference = as.factor(Sriniv.predict_3C$cluster))

print(CM2)
print(CM3)