Ravel = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Ravel_Genus_Counts.rds")
Sriniv = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Sriniv_Genus_Counts.rds")

 
# Select Target
select.target = function(dataset, tget){
  require(dplyr)
  matr = dataset[1:(length(names(dataset)) - 7)] # 7 cause are clinical variables
  target = dataset  %>% select(starts_with(tget))
  dataframe = data.frame(cbind(matr,target))
  names(dataframe)[length(names(dataframe))] <- "target"
  return(dataframe)
}
remove.intermediate = function(dataset){
  df = dataset
  df2 = df[!(df$target == "intermediate"),]
  return(df2)
}
norm.dataset = function(data){
  # Retain only numerics variables
  cols <- sapply(data, is.numeric) 
  
  # Normalize that variables
  data[cols] <- apply(X = data[cols], FUN = function(x) log2(x+1), MARGIN = 2) 
  return(data) 
}

#3C
Ravel_3C = select.target(dataset = Ravel, tget ="Nugent_score_category")
Sriniv_3C = select.target(dataset = Sriniv, tget = "nugent")


#2C
#Ravel_2C = remove.intermediate(dataset = Ravel_3C)
#Sriniv_2C = remove.intermediate(dataset = Sriniv_3C)

#Norm Datasets
Ravel_3C = norm.dataset(Ravel_3C)
Sriniv_3C = norm.dataset(Sriniv_3C)

#Ravel_2C = norm.dataset(Ravel_2C)
#Sriniv_2C = norm.dataset(Sriniv_2C)
Ravel_2C = Ravel_3C
Sriniv_2C = Sriniv_3C

#Clustering
make.cluster = function(dataset, k){
  require(factoextra)
  dataframe = dataset
  cols <- sapply(dataset, is.numeric)
  var <- sapply(dataset, is.character)
  mat = dataset[cols]
  meta = dataset[var]
  set.seed(1312)
  opt = fviz_nbclust(mat, FUNcluster = kmeans) 
  print(opt)
  km.res <- kmeans(mat, k, iter.max = 100, nstart = 25 )
  a = fviz_cluster(km.res,mat, geom = "point")
  dataframe$cluster <- km.res$cluster
  dataframe$cluster <- sub("^", "C", dataframe$cluster)
  table(dataframe$target,dataframe$cluster)
  lista = list(dataframe, km.res, a)
  names = c("dataframe", "km.res", "plot")
  names(lista) = names
  return(lista)
}
make.cl_predict= function(lista, newdata){
  require(clue)
  dataset = lista$dataframe
  km.res = lista$km.res
  cols <- sapply(dataset, is.numeric)
  var <- sapply(dataset, is.character)
  mat = dataset[cols]
  meta = dataset[var]
  ft = colnames(mat)
  new_target = cl_predict(
    object = km.res,
    newdata = newdata[ft])
  cluster = as.factor(as.vector(new_target))
  target = as.factor(newdata$target)
  new_dataframe = cbind(newdata[ft], cluster, target)
  new_dataframe$cluster <- sub("^", "C", new_dataframe$cluster)
  return(new_dataframe)
}
Clust_2C = make.cluster(dataset = Ravel_2C, k = 2)
Clust_3C = make.cluster(dataset = Ravel_3C, k = 3)
table(Clust_2C$dataframe$target,Clust_2C$dataframe$cluster)
table(Clust_3C$dataframe$target, Clust_3C$dataframe$cluster)
Sriniv.predict_2C = make.cl_predict(lista = Clust_2C, newdata = Sriniv_2C)
Sriniv.predict_3C = make.cl_predict(lista = Clust_3C, newdata = Sriniv_3C)
table(Sriniv.predict_2C$target,Sriniv.predict_2C$cluster)
table(Sriniv.predict_3C$target, Sriniv.predict_3C$cluster)


#Train
ML.exec_C2 = function(dataset){
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
  lrn_RF<-makeTuneWrapper(l1,  resampling = inner, par.set = psRF, measures = auc, control=ctrl,  show.info = T)
  
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
  bmr = benchmark(learners, task, outer, measures =  list(acc,auc,mmce), show.info = T, models = T)
  parallelStop()
  return(bmr)
}
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