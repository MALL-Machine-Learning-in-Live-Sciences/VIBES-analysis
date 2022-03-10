# HIL
Ravel_Ge = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Ravel_Sps_Counts.rds")
Sriniv_Ge = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Sriniv_Genus_Counts.rds")

select.target = function(dataset, tget){
  require(dplyr)
  matr = dataset[1:(length(names(dataset)) - 7)] # 7 cause are clinical variables
  target = dataset  %>% select(starts_with(tget))
  dataframe = data.frame(cbind(matr,target))
  names(dataframe)[length(names(dataframe))] <- "target"
  return(dataframe)
}
norm.dataset = function(data){
  # Retain only numerics variables
  cols <- sapply(data, is.numeric) 
  
  # Normalize that variables
  data[cols] <- apply(X = data[cols], FUN = function(x) log2(x+1), MARGIN = 2) 
  return(data) 
}

Ravel_Ge = select.target(dataset = Ravel_Ge, tget ="Nugent_score_category")
Ravel_Ge = norm.dataset(data = Ravel_Ge)

Sriniv_Ge = select.target(dataset = Sriniv_Ge, tget = "nugent")
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
for (i in seq_along(ldata)) {
  saveRDS(ldata[[i]], file = paste0("~/git/BVMetaGenomics/data/SpeciesIntersect/train/Ravel_Sp_HIL_",names(ldata[i]),"_",ncol(ldata[[i]])-1,".rds"))
}

# Benchmarks
set.seed(1312)
ML.exec.3C = function(dataset, identificator){
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
    makeDiscreteParam("mtry", values = c(round(sqrt(ncol(task$env$data))),
                                         round(sqrt(ncol(task$env$data)))+1,
                                         round(sqrt(ncol(task$env$data)))+2,
                                         round(sqrt(ncol(task$env$data)))+3,
                                         round(sqrt(ncol(task$env$data)))+4)),
    makeDiscreteParam("ntree", values= 1000L),
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
bmrs = list()
for (i in seq_along(ldata)){
  bmrs[[i]] = ML.exec.3C(dataset = ldata[[i]], identificator = names(ldata[i]))
}
names(bmrs) = nm
saveRDS(object = bmrs, "~/git/BVMetaGenomics/data/SpeciesIntersect/Benchmarks/BMR_Sp_HIL.rds")
