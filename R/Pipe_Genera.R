# Data Load
path = "~/git/BVMetaGenomics/data/GenusIntersect/"
Ravel_Ge = readRDS(paste0(path,"Ravel_Genus_Counts.rds"))
Sriniv_Ge = readRDS(paste0(path,"Sriniv_Genus_Counts.rds"))

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

Ravel_Ge = select.target(dataset = Ravel_Ge, tget ="Nugent_score_category")
Ravel_Ge = remove.intermediate(dataset = Ravel_Ge)
Ravel_Ge = norm.dataset(data = Ravel_Ge)

Sriniv_Ge = select.target(dataset = Sriniv_Ge, tget = "nugent")
Sriniv_Ge = remove.intermediate(dataset = Sriniv_Ge)
Sriniv_Ge = norm.dataset(data = Sriniv_Ge)


# FS
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
LDM.FS = function(data, variables, targets, seed, method="bray", thold= 0.005){
  require(LDM)
  #ExecuteLDM
  fit.ldm = ldm(variables ~target,
                data = data,
                test.otu = TRUE,
                test.global = TRUE,
                dist.method = method,
                fdr.nominal = thold,
                n.perm.max = 10000,
                seed = seed)
  
  #Filter features
  w1 = which(fit.ldm$q.otu.omni[1,] < thold)
  (n.otu.omni.m1 = length(w1))
  features = (otu.omni.m1 = colnames(fit.ldm$q.otu.omni)[w1])
  
  #Select new features
  filtered.data = variables[,features]
  filtered.data = as.data.frame(cbind(filtered.data, target = targets))
  # Return dataframe only w the features select from LDM
  return(filtered.data)
}
#FCBF
train.FCBF = FCBF.FS(data = Ravel_Ge, thold = 0.025)
#KW
train.KW = KW.FS(data = Ravel_Ge, thold = 90)
#LDM
targets = as.factor(Ravel_Ge$target)
cols = sapply(Ravel_Ge, is.numeric)
variables = Ravel_Ge[cols]
train.LDM = LDM.FS(data = Ravel_Ge,variables = variables,targets = targets, seed = 123)
#Raw
train.raw = Ravel_Ge
nm = c("FCBF", "KW", "LDM", "Raw")
ldata = list(train.FCBF, train.KW, train.LDM, train.raw)
names(ldata) = nm
## Save dataframes
for (i in seq_along(ldata)) {
  saveRDS(ldata[[i]], file = paste0(path,"train/Ravel_Ge_",names(ldata[i]),"_",ncol(ldata[[i]])-1,".rds"))
}


# Benchmarks
set.seed(1312)
ML.exec = function(dataset, identificator){
  require(mlr)
  require(methods)
  require(parallel)
  require(parallelMap)
  cores = 2
  task = makeClassifTask(data = dataset, target = 'target',id = identificator)
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
    makeDiscreteParam("ntree", values= 1000),
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
  
  learners = list(lrn_RF, lrn_glmnet, lrn_GB, lrn_KSVM, lrn_GBM)
  #learners = (lrn_RF)
  # Outer
  outer = makeResampleDesc('RepCV' , reps = 10, folds = 5 , stratify = T)
  
  # Benchmarking
  #parallelStartSocket(cpus = detectCores()*0.5, level = 'mlr.tuneParams')
  parallelStartMulticore(cores , level = 'mlr.tuneParams')
  bmr = benchmark(learners, task, outer, measures =  list(acc,auc,mmce), show.info = T, models = T)
  parallelStop()
  return(bmr)
}
bmrs = list()
for (i in seq_along(ldata)){
  bmrs[[i]] = ML.exec(dataset = ldata[[i]], identi = names(ldata[i]))
}
names(bmrs) = nm
saveRDS(object = bmrs, paste0(path,"Benchmarks/BMR_Ge.rds"))

