# ML  functions
ML.exec = function(dataset){
  require(mlr)
  #require(parallel)
  #require(parallelMap)
  task = makeClassifTask(data = dataset, target = 'target')
  
  # Hyperparameter tuning
  ctrl<-makeTuneControlGrid()
  inner<-makeResampleDesc("Holdout")
 
 
  
  # Random Forest
  psrf<-makeParamSet(
    makeDiscreteParam("mtry", values = c(round(sqrt(ncol(task$env$data)))-2,
                                         round(sqrt(ncol(task$env$data))),
                                         round(sqrt(ncol(task$env$data)))+2)),
    makeDiscreteParam("ntree", values= 1000L),
    makeDiscreteParam("nodesize", values= c(1:3))
  )
  l<-makeLearner("classif.randomForest", predict.type = "prob")
  lrn_rf<-makeTuneWrapper(l,  resampling = inner, par.set = psrf, measures = acc, control=ctrl,  show.info = T)
  
  # GLMNET
  psglmnet = makeParamSet(
    makeDiscreteParam("lambda", c(0.001,0.01,0.1)),
    makeDiscreteParam("alpha",c(0.15,0.25,0.35,0.5))
  )
  l<-makeLearner("classif.glmnet", predict.type = "prob")
  lrn_glmnet<-makeTuneWrapper(l, inner, psglmnet, measures = acc, control = ctrl, show.info=T)
  
  # GBM
  psGBM = makeParamSet(makeDiscreteParam("distribution", values = "huberized"),
                       makeIntegerParam("n.trees", lower = 100, upper = 1000), #number of trees
                       makeIntegerParam("interaction.depth", lower = 2, upper = 10), #depth of tree
                       makeIntegerParam("n.minobsinnode", lower = 10, upper = 80),
                       makeNumericParam("shrinkage",lower = 0.01, upper = 1))
  lrn6 = makeLearner("classif.gbm", predict.type = "response") #ES PROB O RESPONSE??
  lrnGBM =  makeTuneWrapper(learner = lrn6, resampling = inner, measures = auc, par.set = psGBM, control = ctrl, show.info = T)
  learners = list( lrnGBM, lrn_rf,lrn_glmnet)
  
  # Outer
  outer = makeResampleDesc('RepCV' , reps = 5, folds = 10 , stratify = T)

  # Benchmarking
  #parallelStartSocket(cpus = detectCores()*0.5, level = 'mlr.tuneParams')
  bmr = benchmark(learners, task, outer, measures =  list(acc), show.info = T, models = T)
  #parallelStop()
  return(bmr)
}










