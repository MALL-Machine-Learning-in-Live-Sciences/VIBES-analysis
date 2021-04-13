# Phyloseq preprocessing functions

# Preprocess Phyloseq Object
aglomerate = function(phyobject, rank){
  require(phyloseq)
  phy =  tax_glom(phyobject, rank)
  return(phy)
}

relat.abun = function(phyobject){
  require(phyloseq)
  phy = transform_sample_counts(phyobject, function(x) x / sum(x))
  return(phy)
}

prune.OTUs = function(phyobject, pctg = 0.05, count = 0, vari = 0){ 
  require(phyloseq)
  print("Care: u have to change the rank in the fuction if u arent working with genus")
  #Delete unidentified Ranks
  ps <- subset_taxa(phyobject, !Rank6 %in% c("g__"))
  
  #Make and apply the filter
  filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > count), 
                                        A = pctg*nsamples(ps))
  phy <- prune_taxa(filter, ps)
  
  # If we want filter by variance 
  phy = filter_taxa(phy, function(x) var(x) > vari, TRUE)
  return(phy)
}

get.dataset = function(phyobject){
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


split.data = function(data, seed = 123, pctg = 0.90){ #set a seed to reproducibility and  % of train/test
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
  return(data_pool)
}

save.splits = function(path, list, project, id){
  train = list[[1]]
  test = list[[2]]
  saveRDS(train, file =paste0(path,project,"_",id,"_","train.rds"))
  saveRDS(test, file =paste0(path,project,"_",id,"_","test.rds"))
}


# Feature Selection

kruskal.FS = function(data, fs.type = 'kruskal.test', nfeat){
  
  stopifnot('target' %in% names(data))
  
  require(mlr)
  task = makeClassifTask(data = data, target = 'target')
  
  tasks = lapply(nfeat, function(x) filterFeatures(task, method = fs.type, abs = x))
  
  for (i in 1:length(nfeat)) {
    tasks[[i]]$task.desc$id =  paste(fs.type, ncol(tasks[[i]]$env$data) - 1 , sep = "_")
  }
  
  t = list()
  for (i in 1:length(tasks)) {
    t[[i]] = tasks[[i]]$env$data
    names(t)[[i]] = paste(fs.type, nfeat[i], sep = '_')
  }
  
  return(t[[1]])
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

LDM.FS = function(data, seed, method="bray", thold= 0.1){
  targets = as.factor(data$target)
  cols = sapply(data, is.numeric)
  variables = data[cols]
  require(LDM)
  #ExecuteLDM 
  fit.ldm = ldm(variables ~ (target),
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

#cluster.FS = function(){}



# Cosas varias
#Dividimos el phyloseq total en train y test
#set.seed(45)
#size <- floor(0.90 * nrow(data))
#train_ind  = sample(phylo100@sam_data@row.names, size = size)
#a = as.character(phylo100@sam_data@row.names)
#test_ind = setdiff(a, train_ind)
#train_phylo = prune_samples(x = phylo100, samples = train_ind )
#test_phylo = prune_samples(x = phylo100, samples = test_ind )
#saveRDS(train_phylo, file = paste0(path,"train_phylo.rds" ))
#saveRDS(test_phylo, file = paste0(path,"test_phylo.rds" ))


#"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
#names(taxonomyTable)[1] <- "Kingdom"
#names(taxonomyTable)[2] <- "Phylum"
#names(taxonomyTable)[3] <- "Class"
#names(taxonomyTable)[4] <- "Order"
#names(taxonomyTable)[5] <- "Family"
#names(taxonomyTable)[6] <- "Genus"
#names(taxonomyTable)[7] <- "Species"

#Delete marks from names of taxonomy table
#for (j in 1:ncol(taxonomyTable)){
#  nn = list()
#  taxonomyTable[,j] = as.vector(as.character(taxonomyTable[,j]))
#  for (i in 1:nrow(taxonomyTable)) {
#    nn = substr(taxonomyTable[,j], 4,100)
#  }
#  taxonomyTable[,j] = nn
#}

# Indicators of BV: 0-3 (Low), 4-6 (Intermediate), 7-10 (High)
# For numeric values, named arguments can also be used






