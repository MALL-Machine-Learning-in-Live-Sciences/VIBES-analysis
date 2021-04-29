# Funciones Validacion
bench = readRDS("projects/Entropy/data/benchmarks/Ravel_Genus_C_Benchmarks.rds")
modelsAR = readRDS("projects/Entropy/data/benchmarks/Ravel_Genus_AR_Benchmarks.rds")
indice = 2
algoritmo = "classif.ksvm.tuned"
get.features = function(bench, algoritmo, indice){
  require(mlr)
  #Get algorithm and index 
  models = as.data.frame(bench[[indice]])
  model = model[model$learner.id == algoritmo,]
  index = which.max(model$auc)
  
  # Get best model
  modelos = getBMRModels(bench[[indice]])
  alg_index = which(names(modelos$dataset)==algoritmo)
  best = getLearnerModel(modelos$dataset[[alg_index]][[index]])
  
  #Get features
  names = best$features
  return(names)
}


Ravel = readRDS("projects/Entropy/data/Ravel_phyloseq.rds")
Ravel = phy.aglomerate(phyobject = Ravel, rank = "Rank6")
Sriniv = readRDS("projects/Entropy/data/Sriniv_Nugent_phyloseq.rds")
Sriniv = phy.aglomerate(phyobject = Sriniv, rank = "Rank6")


check.features = function(Ravel, Sriniv, FT){
  require(phyloseq)
  Ravel.df = as.data.frame(tax_table(Ravel))
  Ravel.df = Ravel.df[FT,]
  features = as.data.frame(cbind(rownames(Ravel.df),Ravel.df$Rank6))
  Sriniv_sub <- subset_taxa(Sriniv, Rank6 %in% features[,2]) 
  Sriniv.df = as.data.frame(tax_table(Sriniv_sub))
  if (length(features[,2]) == length(Sriniv.df$Rank6)){
    print("All OK")
    return(FT)
  } else if (length(features[,2]) != length(Sriniv.df$Rank6)){
  length(features[,2]) == length(Sriniv.df$Rank6)
  dif = setdiff(features[,2], Sriniv.df$Rank6)
  list.features = list()
  
  for (i in seq_along(dif)){
    list.features[[i]] = features$V1[(features$V2 == dif[i])]
  }
  features = unlist(list.features)
  print("Have to retrain the model without features returned")
  return(features)
  }
}
bncmark = bmFCBF

get.best.model = function(bncmark){
  modelo = getBMRModels(bncmark)
  modelodf = as.data.frame(bncmark)
  
  # Check the model with best results(we asume here that we have only 1 algorithm)
  # We can change this according to algorithms used
  models_index = modelodf[
    order( modelodf[,5], modelodf[,4],  modelodf[,6] ),
  ]
  index = as.numeric(tail(rownames(models_index), 1))
  # Get the model
  best = getLearnerModel(modelo$dataset[[1]][[index]])
  return(best)
}













