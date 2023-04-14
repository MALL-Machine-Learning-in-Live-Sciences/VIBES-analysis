
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