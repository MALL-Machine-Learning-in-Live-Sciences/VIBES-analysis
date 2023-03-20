set.seed(1312)
ML.exec_C3 = function(dataset, identi){
  require(mlr)
  require(methods)
  require(parallel)
  require(parallelMap)
  cores = detectCores()
  task = makeClassifTask(data = dataset, target = 'cluster', id = identi)
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
    makeDiscreteParam("mtry", values = c(round(sqrt(ncol(task$env$data))),
                                         round(sqrt(ncol(task$env$data)))+1,
                                         round(sqrt(ncol(task$env$data)))+2,
                                         round(sqrt(ncol(task$env$data)))+3,
                                         round(sqrt(ncol(task$env$data)))+4)),
    makeDiscreteParam("ntree", values= c(1,5,25,50,75,100,250,500,1000)),
    makeDiscreteParam("nodesize", values= c(1:10))
  )
  l1<-makeLearner("classif.randomForest", predict.type = "prob")
  lrn_RF<-makeTuneWrapper(l1,  resampling = inner, par.set = psRF, measures = acc, control=ctrl,  show.info = T)
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
