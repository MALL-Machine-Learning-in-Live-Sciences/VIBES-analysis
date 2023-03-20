require(phyloseq)
setwd("~/git/BV_Microbiome/")
rank = "Genus" # Genus or "Species"
Ravel = readRDS(paste0("extdata/",rank,"Intersect/Cluster/Ravel_Cluster_",rank,"_pseq.rds"))

# FS
fs.abs = function(pseq, target = 'target', fs.type= 'kruskal.test', ref='KT', nfeat, nombres, destino){
  require(mlr)
  get.dataset = function(phyobject){
    require(phyloseq)
    clinics = as.data.frame(sample_data(phyobject))
    maintain <- c("cluster")
    clinics <- clinics[, maintain]
    otus = as.data.frame(get_taxa(phyobject))
    data = cbind(otus, clinics)
    return(data)
  }
  
  data <- get.dataset(pseq)
  task = makeClassifTask(data = data, target = "cluster")
  tasks = lapply(nfeat, function(x) filterFeatures(task, method = fs.type, abs = x))
  
  for (i in 1:length(nfeat)) {
    tasks[[i]]$task.desc$id =  paste(fs.type, ncol(tasks[[i]]$env$data) - 1 , sep = "_")
  }
  
  t = list()
  for (i in 1:length(tasks)) {
    t[[i]] = tasks[[i]]$env$data
    names(t)[[i]] = paste0(nombres,"_",nfeat[i])
  }
  for (i in 1:length(t)) {
    saveRDS(t[[i]], file = paste0(destino, names(t)[[i]],".rds"))
  }
}

