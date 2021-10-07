#Validation
a = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Benchmarks/BMR_Sp_HIL.rds")

#b = list(a)
# GetBM
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
best_RFs = get.best.RF(bncmark = a)

#Validation
#Test
make.val = function(best_models, data, pos.class = "high"){
  require(caret)
  prds = list()
  CMs = list()
  c = names(best_models)
  for (i in seq_along(best_models)){
    features = best_models[[i]]$features
    features = c(features, "target")
    test = data[,features]
    test_task = makeClassifTask(data = test, target = "target")
    test_task = normalizeFeatures(
      test_task,
      method = "range",
      cols = NULL,
      range = c(0, 1),
      on.constant = "quiet")
    prds[[i]] = predict(best_models[[i]], task = test_task, type = "prob")
    CMs[[i]] = confusionMatrix(data = prds[[i]]$data$response, reference = as.factor(test$target), positive = pos.class)
  }
  names(prds) = c
  names(CMs) = c
  l = list(prds,CMs)
  names(l) = c("Predictions", "CM")
  return(l)
}
test = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/test/Ravel_Counts_Test_HIL.rds")
validation = make.val(best_models = best_RFs, data = test )


#External Validation
Sriniv = readRDS(file = "~/git/BVMetaGenomics/data/SpeciesIntersect/test/Sriniv_Counts_HIL.rds")
norm.dataset = function(data){
  # Retain only numerics variables
  cols <- sapply(data, is.numeric) 
  
  # Normalize that variables
  data[cols] <- apply(X = data[cols], FUN = function(x) log2(x+1), MARGIN = 2) 
  return(data) 
}
Sriniv = norm.dataset(Sriniv)
validation = make.val(best_models = best_RFs, data = Sriniv )
validation$CM

for (i in seq_along(validation$CM)){
  saveRDS(validation$CM[i], paste0("~/git/BVMetaGenomics/data/SpeciesIntersect/Validation/HIL/CM_RF_",names(validation$CM)[i],".rds"))
}

