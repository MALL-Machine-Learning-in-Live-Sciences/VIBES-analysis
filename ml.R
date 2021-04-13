### ScriptMLParalllelRF&GLMNET ###
machnlearn = function(path, patron, outputfile){
  #CARGADATOS
  require(mlr)
  require(parallelMap)
  l = list.files(path, pattern = patron)
  ldata = list()
  for (i in 1:length(l)) {
    ldata[[i]] = readRDS(paste(path, l[[i]], sep=''))
  }
  names(ldata) = l
  
  #TASK
  lista.tareas = list()
  for (i in 1:length(ldata)){
    lista.tareas[[i]] = makeClassifTask(id = paste('nfeat_',ncol(ldata[[i]])-1,sep=""),data = ldata[[i]],target = "target")
  }
  
  #SIN CF
  lista.sinCF = list()
  for (i in 1:length(lista.tareas)){
    lista.sinCF[[i]] = removeConstantFeatures(lista.tareas[[i]])
  }
  lista.Norm = list()
  for (i in 1:length(lista.sinCF)){
    lista.Norm[[i]] = normalizeFeatures(
      lista.sinCF[[i]],
      method = "range",
      cols = NULL,
      range = c(0, 1),
      on.constant = "quiet")
  }
  ###################################################################################
  bnchmarks = list()
  #parallelStartMulticore(20L , level = 'mlr.tuneParams')
  #parallelStarSocket()   PARA  USARLO EN WINDOWS
  
  for (i in 1:length(lista.Norm)){
    ### Busqueda de hiperparámetros ###
    control = makeTuneControlGrid() # Establecer un control par los parametros del modelo
    inner = makeResampleDesc(method = "Holdout")#Optimizacion de los parametros del modelo usando un Holdout para despues utilizar un CV
    #RF#
    psRF = makeParamSet(makeDiscreteParam("mtry", values = c(round(sqrt(ncol(lista.Norm[[i]]$env$data)))-2,
                                                             round(sqrt(ncol(lista.Norm[[i]]$env$data))),
                                                             round(sqrt(ncol(lista.Norm[[i]]$env$data)))+2)),#Sumar dos y restar dos para tener un rango 
                        makeDiscreteParam("ntree", values= 1000L),
                        makeDiscreteParam("nodesize", values= c(1:3)))
    lrn2 = makeLearner("classif.randomForest",  predict.type = "prob")
    lrnRF = makeTuneWrapper(learner = lrn2, resampling = inner, measures = acc, par.set = psRF, control = control, show.info = T)
    
    #GLMNET#
    psGL = makeParamSet(makeDiscreteParam("lambda", values = c(0.0001,0.001,0.01,0.1,1)),
                        makeDiscreteParam("alpha",values = c(0,0.15,0.25,0.35,0.5,0.65,0.75,0.85,1)))
    lrn3 = makeLearner("classif.glmnet", predict.type = "prob")
    lrnGL = makeTuneWrapper(learner = lrn3, resampling = inner, measures = acc, par.set = psGL, control = control, show.info = T)
    
    learners = list(lrnRF)
    outer = makeResampleDesc(method = 'RepCV', predict = 'both', reps = 5, folds = 10, stratify = TRUE)
    #parallelStop()
    ## Benchmarking ##
    bnchmarks[[i]] = benchmark(learners, lista.Norm[[i]], outer, measures = list(acc), show.info = T , models = T)
  }
  #bmrk = mergeBenchmarkResults(bnchmarks)
  #saveRDS(bmrk, file= paste(outputfile,"Benchmark",patron,sep=""))
  return(bnchmarks)
}