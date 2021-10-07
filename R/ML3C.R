#ML3C
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
path = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy2/Data/GenusIntersect/Cluster/3C/FS/"
patron = "Ravel"
l = list.files(path, pattern = patron)
ldata = list()
for (i in 1:length(l)) {
  ldata[[i]] = readRDS(paste(path, l[[i]], sep=''))
}
df = data.frame("data"=l,"data2"=1:23, stringsAsFactors = FALSE)
df$data = substr(df$data,1,nchar(df$data)-4)
l = df$data
names(ldata) = l

#ML
bmrs = list()
for (i in seq_along(ldata)){
  bmrs[[i]] = ML.exec_C3(dataset = ldata[[i]], identi = names(ldata[i]))
}
names(bmrs) = l
saveRDS(bmrs, "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy2/Data/GenusIntersect/Cluster/3C/BMR1-23.rds")

#GetBR
get.best.RF = function(bncmark){
  require(mlr)
  best_models = list()
  c = names(bncmark)
  for (i in seq_along(bncmark)){
    modelo = getBMRModels(bncmark[[i]])
    modelodf = as.data.frame(bncmark[[i]])[1:50,]
    models_index = modelodf[
      order( modelodf[,4], modelodf[,5] ),
    ]
    index = as.numeric(tail(rownames(models_index), 1))
    # Get the model
    best_models[[i]] = getLearnerModel(modelo[[1]][[1]][[index]])
  }
  names(best_models) = c
  return(best_models)
}
best_RFs = get.best.RF(bncmark = bmrs)

#Val
make.val = function(best_models, data){
  require(caret)
  prds = list()
  CMs = list()
  c = names(best_models)
  for (i in seq_along(best_models)){
    features = best_models[[i]]$features
    features = c(features, "cluster")
    test = data[,features]
    test_task = makeClassifTask(data = test, target = "cluster")
    test_task = normalizeFeatures(
      test_task,
      method = "range",
      cols = NULL,
      range = c(0, 1),
      on.constant = "quiet")
    prds[[i]] = predict(best_models[[i]], task = test_task, type = "prob")
    CMs[[i]] = confusionMatrix(data = prds[[i]]$data$response, reference = as.factor(test$cluster))
  }
  names(prds) = c
  names(CMs) = c
  l = list(prds,CMs)
  names(l) = c("Predictions", "CM")
  return(l)
}
Sriniv = readRDS("/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy2/Data/GenusIntersect/Cluster/3C/Sriniv_Clust3.rds")
Sriniv$target <- NULL
val = make.val(best_models = best_RFs, data = Sriniv)

saveRDS(val,"/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy2/Data/GenusIntersect/Cluster/3C/VAL1-23.rds")

#  library(mlr)
#  b = readRDS("Entropy2/Data/SpeciesIntersect/Cluster/2C/BMR1-27.rds")
#  val = readRDS("Entropy2/Data/SpeciesIntersect/Cluster/2C/VAL1-27.rds")
#  #Train
#  n = names(b)
#  ACC = list()
#  for (i in seq_along(b)){
#    ACC[[i]] = b[[i]]$results[[1]]$classif.randomForest.tuned$aggr[[1]]
#  }
#  df = data.frame("N.Feat" = n, "ACC" = unlist(ACC))
#  df$N.Feat = as.double(substr(df$N.Feat,10,12))
#  df = df[order(df$N.Feat),]
# library(ggplot2)
# library(viridis)
#  p = ggplot(df, aes(x = N.Feat, y = ACC)) +
#    geom_line(color= viridis(3)[2]) + geom_point(color= viridis(1))+
#    theme_light()+
#    theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.ticks.x = element_blank())
#  p
# #
# # #VAL
#  ACC = list()
#  up = list()
#  low = list()
#  nomb = names(val$CM)
#  for (i in seq_along(val$CM)){
#    ACC[[i]] = val$CM[[i]]$overall[[1]]
#    low[[i]] = val$CM[[i]]$overall[[3]]
#    up[[i]] = val$CM[[i]]$overall[[4]]
#  }
#  df = data.frame("N.Feat" = nomb, "ACC" = unlist(ACC), Low = unlist(low), Upp = unlist(up))
#  df$N.Feat = as.double(substr(df$N.Feat,10,12))
#  df = df[order(df$N.Feat),]
# 
# 
#  p = ggplot(df, aes(x = N.Feat, y = ACC))+
#     geom_line(color= viridis(3)[2]) + geom_point(color= viridis(1))+
#    #  geom_segment(aes(x = 7, y = 0.98, xend = 7.9, yend = df$ACC[[8]]+0.005),
#    #              arrow = arrow(length = unit(0.3, "cm")),
#    #              colour = viridis(3)[3])+
#     #geom_point(aes(x=8, y=df$ACC[[8]]), shape=23,fill=viridis(3)[1],color=viridis(3)[3],size=3)+
#     theme_light()+
#     theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.ticks.x = element_blank())
#   p
# 
# val$CM$Ravel_2C_9
