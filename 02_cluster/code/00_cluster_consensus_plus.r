packgs <- c("phyloseq", "ConsensusClusterPlus")
lapply(packgs, require, character.only = TRUE)
setwd("~/git/BV_Microbiome/")
rank <-  "Species" # Genus or "Species"
nfeat <- "22" # 24 and 14 for Species or 37 and 21 for Genus
trans <- "clr" # "log2", "clr" or "alr"
cl <- "km" # "km", "pam" or "hc"

# Select Target
select_target <- function(pseq, tget) {
  require(dplyr)
  require(phyloseq)
  df <- sample_data(pseq)
  colnames(df)[colnames(df) == tget] <- "target"
  sample_data(pseq) <- df
  return(pseq)
}
# Normalization
if (trans == "log2") {
  norm_dataset <- function(pseq) {
    # Change columns by rows too, interested in maintain fts in columns
    require(phyloseq)
    otu <- data.frame(otu_table(pseq))
    # Normalize that variables
    otu <- apply(X = otu, FUN = function(x) log2(x + 1), MARGIN = 2)
    otu_table(pseq) <- otu_table(otu, taxa_are_rows = FALSE)
    print(paste("log2", "normalization selected"))
    return(pseq)
  }
  
}else if (trans == "clr") {
  norm_dataset <- function(pseq){
    require(microbiome)
    # Note that small pseudocount is added if data contains zeroes
    pseq_clr <-  microbiome::transform(pseq, transform =  'clr', shift=1)
    print(paste("clr", "normalization selected"))
    return(pseq_clr)
  }
}else if (trans == "alr") {
  norm_dataset = function(pseq){
    require(microbiome)
    pseq_alr =  microbiome::transform(pseq, transform =  'alr', shift=1,
                                      reference=1)
    print(paste("alr", "normalization selected"))
    return(pseq_alr)
  }
}else{
  print("Introduce valid normalization (log2, clr or alr)")
}

#Load data
Ravel <- readRDS(paste0("00_preprocess_cohorts/data/", rank, "Intersect/Ravel_",
                        rank, "_pseq_", nfeat, ".rds"))
Ravel <- select_target(pseq = Ravel, tget = "Nugent_score_category")
Ravel <- norm_dataset(pseq = Ravel)
ravel_mat <- t(otu_table(Ravel)@.Data) # samples on columns for CCP

# CCP
require(ConsensusClusterPlus)
title <- paste0("CCP_Ravel_", rank, "_", nfeat, "_", trans, "_", cl)
setwd("~/git/BV_Microbiome/02_cluster/res/")
ccp = ConsensusClusterPlus(d = ravel_mat, maxK = 6, reps = 1500, pItem = 0.8,
                           pFeature = 1, title = title, clusterAlg = cl,
                           distance = "euclidean", seed = 1580 ,
                           plot = "pdf")
icl <- calcICL(res = ccp, title = title, plot = "pdf")

# Select k
select_k2 <- function(data, res, maxK){
  require(cluster)
  require(fpc)
  d <- data
  clusters <- res
  
  stcl <- lapply(2:maxK, function(i) cluster.stats(dist(t(d)),clusters[[i]]$consensusClass))
  ldunn <- sapply(1:(maxK-1), function(i) stcl[[i]]$dunn )
  lwbr  <- sapply(1:(maxK-1), function(i) stcl[[i]]$wb.ratio )
  lch   <- sapply(1:(maxK-1), function(i) stcl[[i]]$ch)
  lsil <- vector("list",(maxK-1))
  
  for(i in 2:maxK){
    sil <- silhouette(clusters[[i]]$consensusClass,dist(t(d),method = "euclidean"))
    sizes <- table(clusters[[i]]$consensusClass)
    lsil[[i-1]] <- sil
  }
  
  msil <- sapply(1:(maxK-1), function(i) mean( lsil[[i]][,3] ) )
  cdl <- lapply(2:maxK, function(i) as.dist(1-clusters[[i]]$consensusMatrix ) )
  md <- dist(t(d),method = "euclidean")
  corl <- sapply(cdl, cor,md)
  
  co <- rep(1,(maxK-1))
  nclust.co <- which.max(corl) + 1  # cophenetic distance
  
  co <- rep(1,(maxK-1))
  nclust.sil <- which.max(msil) + 1  # silhouette distance
  
  co <- rep(1,(maxK-1))
  nclust.dunn <- which.max(ldunn) + 1  # dunn index
  
  indexes <- c(nclust.co, nclust.sil, nclust.dunn)
  names(indexes) <- c('cophenetic', 'silhouette', 'dunn')
  
  nclust <- as.numeric(table(indexes)[as.numeric(which.max(table(indexes)))])
  if (nclust == 1) {
    nclust <- nclust.co
  }
  
  return(nclust)
  
}
a <- select_k2(data = ravel_mat, res = ccp, maxK = 6)
a
