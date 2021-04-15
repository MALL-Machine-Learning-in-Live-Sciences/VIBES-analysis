# ML  functions
ML.exec = function(dataset){
  require(mlr)
  require(methods)
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
  learners = list(lrn_glmnet)
  
  # Outer
  outer = makeResampleDesc('RepCV' , reps = 3, folds = 10 , stratify = T)

  # Benchmarking
  #parallelStartSocket(cpus = detectCores()*0.5, level = 'mlr.tuneParams')
  bmr = benchmark(learners, task, outer, measures =  list(acc), show.info = T, models = T)
  #parallelStop()
  return(bmr)
}










