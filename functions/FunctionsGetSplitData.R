#Functions get/split Dataset

get.dataset = function(phyobject){
  require(phyloseq)
  clinics = as.data.frame(sample_data(phyobject))
  otus = as.data.frame(t(get_taxa(phyobject)))
  data = cbind(otus, clinics)
  cols <- sapply(data, is.character)
  names(data)[cols] <- "target"
  return(data)
} 

# Preprocess Generated Dataset
norm.dataset = function(data){
  # Retain only numerics variables
  cols <- sapply(data, is.numeric) 
  
  # Normalize that variables
  data[cols] <- apply(X = data[cols], FUN = function(x) log2(x+1), MARGIN = 2) 
  return(data) 
}


split.data = function(data,seed = 123, pctg = 0.90, path, project, id){ #set a seed to reproducibility and  % of train/test
  require(caret)
  set.seed(seed)
  tr_size <- floor(pctg * nrow(data))
  train_ind = sample(nrow(data), size =tr_size)
  train <- data[train_ind, ]
  test <- data[-train_ind, ]
  data_pool <- list(train,test)
  c = c("train", "test")
  names(data_pool) = c
  # Return a list w train/test
  saveRDS(train, file =paste0(path,project,"_",id,"_","train.rds"))
  saveRDS(test, file =paste0(path,project,"_",id,"_","test.rds"))
  return(data_pool)
  
}

# save.splits = function(path, list, project, id){
#   train = list[[1]]
#   test = list[[2]]
#   saveRDS(train, file =paste0(path,project,"_",id,"_","train.rds"))
#   saveRDS(test, file =paste0(path,project,"_",id,"_","test.rds"))
# }
