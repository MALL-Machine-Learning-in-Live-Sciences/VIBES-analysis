packgs <- c("phyloseq", "ConsensusClusterPlus")
lapply(packgs, require, character.only = TRUE)
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
setwd("~/git/BV_Microbiome/figures/plots/")
ccp = ConsensusClusterPlus(d = ravel_mat, maxK = 6, reps = 1500, pItem = 0.8,
                           pFeature = 1, title = title, clusterAlg = cl,
                           distance = "euclidean", seed = 1580 ,
                           plot = "pdf")
icl <- calcICL(res = ccp, title = title, plot = "pdf")
