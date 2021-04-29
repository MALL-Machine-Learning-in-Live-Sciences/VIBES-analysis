# Funciones Validacion

get.features = function(bench, algoritmo, indice){
  require(mlr)
  #Get algorithm and index 
  models = as.data.frame(bench[[indice]])
  model = models[models$learner.id == algoritmo,]
  index = which.max(model$auc)
  
  # Get best model
  modelos = getBMRModels(bench[[indice]])
  alg_index = which(names(modelos$dataset)==algoritmo)
  best = getLearnerModel(modelos$dataset[[alg_index]][[index]])
  
  #Get features
  names = best$features
  return(names)
}

require(phyloseq)


check.features = function(Sriniv.df, features){

  if (length(features[,2]) == length(Sriniv.df$Rank6)){
    print("All OK")
    return(features)
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

get.sriniv.test = function(Sriniv_sub, Sriniv.df, features){
  require(phyloseq)
  require(dplyr)
  require(tidyverse)
  source("git/Entropy/functions/FunctionsGetSplitData.R")
  #Creo una tabla de correspondencia entre EL OTU y los generos traidos de Ravel
  features2 = as.data.frame(cbind(rownames(Sriniv.df), Sriniv.df$Rank6))
  #Recupero el dataset Matcheo las dos tablas de correspondencia y me quedo con los number acces
  Sriniv_sub = get.dataset(phyobject = Sriniv_sub)
  # Me quedo con los targets
  targets = Sriniv_sub$target
  # Cojo solo las columnas que sean numericas
  cols <- sapply(Sriniv_sub, is.numeric)
  tt = Sriniv_sub[cols]
  tt = as.data.frame(t(tt))
  # Matcheo los features de ambos datasets
  equi = merge(features2,features, by = "V2",sort = FALSE)
  tt <- tibble::rownames_to_column(tt, "V1.x")
  j = merge(tt, equi,by = "V1.x",sort = FALSE)
  aa = j %>% remove_rownames %>% column_to_rownames(var="V1.y")
  #Borro las columnas innecesarias
  test.Sriniv <- subset( aa, select = -c(V2,V1.x))
  # Añado los targets
  test.Sriniv = as.data.frame(cbind(as.data.frame(t(test.Sriniv)),
                                    target =targets ))
  test.Sriniv = norm.dataset(test.Sriniv)
  return(test.Sriniv)
}













