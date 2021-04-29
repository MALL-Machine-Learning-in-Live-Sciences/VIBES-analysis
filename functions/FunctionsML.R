# ML Functions
ML.exec = function(dataset){
  require(mlr)
  require(methods)
  require(parallel)
  require(parallelMap)
  cores = 4
  task = makeClassifTask(data = dataset, target = 'target')
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
  psGB = makeParamSet(makeNumericParam("eta", lower = 0, upper = 1),
                      makeIntegerParam("max_depth", lower = 1, upper = 20),
                      makeDiscreteParam("eval_metric", "logloss"))
  l3 = makeLearner("classif.xgboost", predict.type = "prob", nrounds=10)
  lrn_GB = makeTuneWrapper(learner = l3, resampling = inner, measures = auc, par.set = psGB, control = ctrl, show.info = T)
  
  # SVM
  psKSVM = makeParamSet(makeDiscreteParam('C', values = 2^c(-8, -4, -2, 0)),
                        makeDiscreteParam('sigma', values = 2^c(-8, -4, 0, 4)))
  l4 = makeLearner("classif.ksvm", predict.type = "prob")
  lrn_KSVM = makeTuneWrapper(learner = l4, resampling = inner, measures = auc, par.set = psKSVM, control = ctrl, show.info = T)
  
  # GBM
  psGBM = makeParamSet(makeDiscreteParam("distribution", values = "bernoulli"),
                       makeIntegerParam("n.trees", lower = 100, upper = 800), 
                       makeIntegerParam("interaction.depth", lower = 2, upper = 10),
                       makeNumericParam("bag.fraction", lower = 0.80, upper = 0.80))
  l5 = makeLearner("classif.gbm", predict.type = "prob")
  lrn_GBM =  makeTuneWrapper(learner = l5, resampling = inner, measures = auc, par.set = psGBM, control = ctrl, show.info = T)
  
  #learners = list(lrn_RF, lrn_glmnet, lrn_GB, lrn_KSVM, lrn_GBM)
  learners = (lrn_KSVM)
  # Outer
  outer = makeResampleDesc('RepCV' , reps = 5, folds = 3 , stratify = T)

  # Benchmarking
  #parallelStartSocket(cpus = detectCores()*0.5, level = 'mlr.tuneParams')
  parallelStartMulticore(cores , level = 'mlr.tuneParams')
  bmr = benchmark(learners, task, outer, measures =  list(acc,auc,mmce), show.info = T, models = T)
  parallelStop()
  return(bmr)
}










