#MMUPIH
norm.dataset = function(data){
  # Retain only numerics variables
  cols <- sapply(data, is.numeric) 
  
  # Normalize that variables
  data[cols] <- apply(X = data[cols], FUN = function(x) log2(x+1), MARGIN = 2) 
  return(data)
}
MMUPIH.combat = function(dataset){
  require(MMUPHin)
  target = dataset$target
  batch = dataset$Study
  meta = data.frame(study = as.factor(dataset$Study),
                    target = dataset$target)
  rownames(meta) = rownames(dataset)
  cols <- sapply(dataset, is.numeric)
  mat = as.matrix(dataset[cols])
  MMUPIH <- adjust_batch(feature_abd = t(mat),
                         batch = "study",
                         data = meta,
                         control = list(maxit = 5000, verbose = TRUE))
  MMUPIH.features <- MMUPIH$feature_abd_adj
  mat = as.data.frame(t(MMUPIH.features))
  dataframe = as.data.frame(cbind(mat, target, batch))
  #dataframe = norm.dataset(dataframe)
  return(dataframe)
}
make.cluster = function(dataset, k){
  require(factoextra)
  dataframe = dataset
  cols <- sapply(dataset, is.numeric)
  var <- sapply(dataset, is.character)
  mat = dataset[cols]
  meta = dataset[var]
  opt = fviz_nbclust(mat, FUNcluster = kmeans) 
  km.res <- kmeans(mat, k, iter.max = 1000, nstart = 25)
  a = fviz_cluster(km.res,mat, geom = "point" )
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
show.pca = function(dataset, batch = "batch", cond = "target"){
  require(PCAtools)
  cols <- sapply(dataset, is.numeric)
  var <- sapply(dataset, is.character)
  mat = dataset[cols]
  meta = dataset[var]
  p = PCAtools::pca(mat = t(mat), metadata = meta , removeVar = 0.1)
  biplot(p, drawConnectors = FALSE,labSize = 0,showLoadings = FALSE,
         colby = batch, shape = cond, legendPosition = "right")
}
ML.exec_C2 = function(dataset){
  require(mlr)
  require(methods)
  require(parallel)
  require(parallelMap)
  drops <- c("batch")
  dataset = dataset[ , !(names(dataset) %in% drops)]
  cores = detectCores()
  task = makeClassifTask(data = dataset, target = 'target')
  #task = normalizeFeatures(
  #  task,
  #  method = "range",
  #  cols = NULL,
  #  range = c(0, 1),
  #  on.constant = "quiet")
  
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
get.best.model_C2 = function(bncmark){
  modelo = getBMRModels(bncmark)
  modelodf = as.data.frame(bncmark)
  
  # Check the model with best results(we asume here that we have only 1 algorithm)
  # We can change this according to algorithms used
  models_index = modelodf[
    order(modelodf[,6], modelodf[,4], modelodf[,5]),
  ]
  index = as.numeric(tail(rownames(models_index), 1))
  # Get the model
  best = getLearnerModel(modelo$dataset[[1]][[index]])
  return(best)
}
draw_confusion_matrix <- function(cm) {
  require(viridis)
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  # create the matrix (Care with pos and  neg class, see in Confusion matrix the place ocupied)
  rect(150, 430, 240, 370, col=viridis(2)[1])
  text(195, 445, 'C1', cex=1.2)
  rect(250, 430, 340, 370, col=viridis(3)[2])
  text(295, 445, 'C2', cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=viridis(3)[2])
  rect(250, 305, 340, 365, col=viridis(2)[1])
  text(140, 400, 'C1', cex=1.2, srt=90)
  text(140, 335, 'C2', cex=1.2, srt=90)
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "Details", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1)
  text(30, 85, names(cm$byClass[2]), cex=1, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1)
  text(50, 85, names(cm$byClass[5]), cex=1, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1)
  text(70, 85, names(cm$byClass[6]), cex=1, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1)
  text(90, 85, names(cm$byClass[7]), cex=1, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1)
  text(70, 35, names(cm$overall[2]), cex=1, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1)
} 
KW.FS = function(data, fs.type = 'kruskal.test', pctge = NULL, thold = NULL){
  require(mlr)
  require(testthat)
  task = makeClassifTask(data = data, target = 'target')
  fv = generateFilterValuesData(task, method = fs.type)
  print(fv)
  if (is.null(pctge)){
    filtered.task = filterFeatures(task, fval = fv, threshold = thold)
    filtered.data = filtered.task$env$data
  } else if(is.null(thold)){
    filtered.task = filterFeatures(task, fval = fv, perc = pctge)
    filtered.data = filtered.task$env$data
  }else{
    print("Uncorrect format, choose pctge or thold")
  }
  return(filtered.data)
}
FCBF.FS = function(data, thold){
  require(FCBF)
  #Split target from other variables
  targets = as.factor(data$target)
  cols <- sapply(data, is.numeric)
  variables = data[cols]
  
  # Discretize RA from OTUS in High and LOW. The transpose is done since the function wants the OTUS in the rows.
  discretization = discretize_exprs(t(variables))
  #su_plot(discretization, targets)
  
  #Execute Algm
  fcbf = fcbf(discretization, targets, verbose = T, thresh = thold)
  
  #Select new features
  filtered.data = variables[,fcbf$index]
  filtered.data = as.data.frame(cbind(filtered.data, target = targets))
  # Return dataframe only w the features select from FCBF
  return(filtered.data)
}

# Load Data
All_data_21 = readRDS("/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/All_data_21_RA.rds")
#All_data_21 = readRDS("projects/Entropy/data/All_data_21_RA.rds")
Ravel_data = subset(All_data_21, Study == "Ravel")
Sriniv_data = subset(All_data_21, Study == "Srinivasan")

#Apply MMUPIH
mmu.data = MMUPIH.combat(All_data_21)
Ravel_data = subset(mmu.data, batch == "Ravel")
Sriniv_data = subset(mmu.data, batch == "Srinivasan")

#Train
set.seed(1312)
bmr_2C = ML.exec_C2(dataset = Ravel_data)
#Predict
require(caret)
bm2 =get.best.model_C2(bncmark = bmr_2C)
print(bm2)
test_task2 = makeClassifTask(data = Sriniv_data[1:22], target = "target")
predict2 = predict(bm2, task = test_task2, type = "prob")
CM2 = confusionMatrix(data = predict2$data$response,  reference = as.factor(Sriniv_data$target))
print(CM2)

#FS
R_FCBF = FCBF.FS(Ravel_data[1:22],thold = 0.0025)
R_KW = KW.FS(Ravel_data[1:22],thold = 90)
#train
bmrFC = ML.exec_C2(dataset = R_FCBF)
bmrKW = ML.exec_C2(dataset = R_KW)
#Predict
#FCBF
print("FCBF")
bm2 =get.best.model_C2(bncmark = bmrFC)
print(bm2)
c = colnames(R_FCBF)
test_task2 = makeClassifTask(data = Sriniv_data[c], target = "target")
predict2 = predict(bm2, task = test_task2, type = "prob")
CM2 = confusionMatrix(data = predict2$data$response, reference = as.factor(Sriniv_data$target))
print(CM2)

#KW
print("KW")
bm2 =get.best.model_C2(bncmark = bmrKW)
print(bm2)
c = colnames(R_KW)
test_task2 = makeClassifTask(data = Sriniv_data[c], target = "target")
predict2 = predict(bm2, task = test_task2, type = "prob")
CM2 = confusionMatrix(data = predict2$data$response, reference = as.factor(Sriniv_data$target))
print(CM2)
