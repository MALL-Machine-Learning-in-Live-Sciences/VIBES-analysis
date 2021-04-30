# Functions Feature Selection 

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

LDM.FS = function(data, variables, targets, seed, method="bray", thold= 0.05){
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


CLUST.FS = function(data){
  require(dynamicTreeCut)
  # Remove target column
  datos = select(data, -(target))
  
  #Data Scale
  sdata = scale(x = as.matrix(datos))
  
  # Dissimilarity matrix
  d <- dist(t(sdata), method = "euclidean")
  
  # Hierarchical clustering using Complete Linkage
  hc1 <- hclust(d, method = "complete" )
  
  # Plot the obtained dendrogram
  plot(hc1, cex = 0.6, hang = -1)
  
  #Dinamic tree
  dtree = cutreeDynamic(dendro = hc1, cutHeight = 20, minClusterSize = 3,method = "tree")
  # Extract names for plots
  kk = as.logical(dtree)
  fs_names = colnames(datos)
  fs_names = fs_names[kk]
  # Extract names for plots
  colnames(datos) = as.character(dtree)
  prefixes = unique(sub("\\..*", "", colnames(datos)))
  data_clustered = as.data.frame(sapply(prefixes, function(x) rowSums(datos[,startsWith(colnames(datos), x)])))
  data_clustered = select(data_clustered, -("0"))
  colnames(data_clustered) <- paste("Clu", colnames(data_clustered), sep = "_")
  data_clustered = as.data.frame(cbind(data_clustered, target = data$target))
  #Save names for plots
  #saveRDS(object = fs_names, file = "projects/Entropy/data/CLUST_C_names.rds")
  return(data_clustered)
}
