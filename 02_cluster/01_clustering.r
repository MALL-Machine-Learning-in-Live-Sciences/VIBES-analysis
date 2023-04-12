packgs <- c("phyloseq")
lapply(packgs, require, character.only = TRUE)
setwd("~/git/BV_Microbiome/")
rank <-  "Species" # Genus or "Species"
nfeat <- "22" # 24 and 14 for Species or 37 and 21 for Genus
trans <- "clr" # "log2", "clr" or "alr"
k <- 3

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
    pseq_clr <-  microbiome::transform(pseq, transform =  'clr')
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

Ravel <- readRDS(paste0("extdata/", rank, "Intersect/Ravel_",
                        rank, "_pseq_", nfeat, ".rds"))
Sriniv <-readRDS(paste0("extdata/", rank, "Intersect/Sriniv_",
                        rank, "_pseq_", nfeat, ".rds"))
PR3020 <-readRDS(paste0("extdata/", rank, "Intersect/PRJNA3020_",
                        rank, "_pseq_", nfeat, ".rds"))
PR7977 <-readRDS(paste0("extdata/", rank, "Intersect/PRJNA7977D0_",
                        rank, "_pseq_", nfeat, ".rds"))
PR2085 <-readRDS(paste0("extdata/", rank, "Intersect/PRJNA2085D0_",
                        rank, "_pseq_", nfeat, ".rds"))

#3C
Ravel = select_target(pseq = Ravel, tget ="Nugent_score_category")
Sriniv = select_target(pseq = Sriniv, tget = "nugent")

#Norm Datasets
Ravel <-  norm_dataset(Ravel)
Sriniv <-  norm_dataset(Sriniv)
PR3020 <-  norm_dataset(PR3020)
PR7977 <- norm_dataset(PR7977)
PR2085 <- norm_dataset(PR2085)

#Clustering functions
make_cluster <- function(pseq, n_clusters){
  require(factoextra)
  mat = as.matrix(data.frame(otu_table(pseq)))
  meta = sample_data(pseq)
  set.seed(1580)
  km.res <- kmeans(mat, n_clusters, iter.max = 1000, nstart = 25)
  a = fviz_cluster(km.res,mat, geom = "point")
  identical(names(km.res$cluster), rownames(sample_data(pseq)))
  sample_data(pseq)$cluster <- km.res$cluster
  sample_data(pseq)$cluster <- sub("^", "C", sample_data(pseq)$cluster)
  table(sample_data(pseq)$target , sample_data(pseq)$cluster)
  lista = list(pseq, km.res, a)
  names = c("pseq", "km.res", "plot")
  names(lista) = names
  return(lista)
}
make_cl_predict <- function(lista, pseq2){
  require(clue)
  pseq = lista$pseq
  km.res = lista$km.res
  mat = as.matrix(data.frame(otu_table(pseq)))
  meta = sample_data(pseq)
  ft = colnames(mat)
  mat2 = as.matrix(data.frame(otu_table(pseq2)))
  new_target = cl_predict(
    object = km.res,
    newdata = mat2[,ft])
  sample_data(pseq2)$cluster <- new_target
  sample_data(pseq2)$cluster <- sub("^", "C", sample_data(pseq2)$cluster)
  return(pseq2)
}
# Make clusters
Clust_3C <- make_cluster(pseq = Ravel, n_clusters = k)
Ravel <- Clust_3C$pseq
#Predict clusters
table(sample_data(Ravel)$target , sample_data(Ravel)$cluster)
Sriniv <- make_cl_predict(lista = Clust_3C, pseq2 = Sriniv)
table(sample_data(Sriniv)$target , sample_data(Sriniv)$cluster)
PR3020 <- make_cl_predict(lista = Clust_3C, pseq2 = PR3020)
table(sample_data(PR3020)$cluster)
PR7977 <- make_cl_predict(lista = Clust_3C, pseq2 = PR7977)
table(sample_data(PR7977)$cluster)
PR2085 <- make_cl_predict(lista = Clust_3C, pseq2 = PR2085)
table(sample_data(PR2085)$cluster, sample_data(PR2085)$ABV)


#Rename clusters
# Ravel
sample_data(Ravel)$cluster[sample_data(Ravel)$cluster == "C1"] <- "D"
sample_data(Ravel)$cluster[sample_data(Ravel)$cluster == "C2"] <- "N"
sample_data(Ravel)$cluster[sample_data(Ravel)$cluster == "C3"] <- "ID"
#sample_data(Ravel)$cluster[sample_data(Ravel)$cluster == "C4"] <- "N"
# Sriniv
sample_data(Sriniv)$cluster[sample_data(Sriniv)$cluster == "C1"] <- "D"
sample_data(Sriniv)$cluster[sample_data(Sriniv)$cluster == "C2"] <- "N"
sample_data(Sriniv)$cluster[sample_data(Sriniv)$cluster == "C3"] <- "ID"
#sample_data(Sriniv)$cluster[sample_data(Sriniv)$cluster == "C4"] <- "N"
# PR3020
sample_data(PR3020)$cluster[sample_data(PR3020)$cluster == "C1"] <- "D"
sample_data(PR3020)$cluster[sample_data(PR3020)$cluster == "C2"] <- "N"
sample_data(PR3020)$cluster[sample_data(PR3020)$cluster == "C3"] <- "ID"
#sample_data(PR3020)$cluster[sample_data(PR3020)$cluster == "C4"] <- "N"
# PR7977
sample_data(PR7977)$cluster[sample_data(PR7977)$cluster == "C1"] <- "D"
sample_data(PR7977)$cluster[sample_data(PR7977)$cluster == "C2"] <- "N"
sample_data(PR7977)$cluster[sample_data(PR7977)$cluster == "C3"] <- "ID"
#sample_data(PR7977)$cluster[sample_data(PR7977)$cluster == "C4"] <- "N"
# PR2085
sample_data(PR2085)$cluster[sample_data(PR2085)$cluster == "C1"] <- "D"
sample_data(PR2085)$cluster[sample_data(PR2085)$cluster == "C2"] <- "N"
sample_data(PR2085)$cluster[sample_data(PR2085)$cluster == "C3"] <- "ID"
#sample_data(PR2085)$cluster[sample_data(PR2085)$cluster == "C4"] <- "N"

# Saving
saveRDS(object = Ravel, file = paste0("extdata/", rank, "Intersect/Cluster/C",
                                      k,"/Ravel_Cluster_", rank, "_", nfeat,
                                      "_pseq.rds"))
saveRDS(object = Sriniv, file = paste0("extdata/", rank, "Intersect/Cluster/C",
                                       k, "/Sriniv_Cluster_", rank, "_", nfeat,
                                       "_pseq.rds"))
saveRDS(object = PR3020, file = paste0("extdata/", rank, "Intersect/Cluster/C",
                                       k, "/PRJNA3020_Cluster_", rank, "_",
                                       nfeat, "_pseq.rds"))
saveRDS(object = PR7977, file = paste0("extdata/", rank, "Intersect/Cluster/C",
                                       k, "/PRJNA7977D0_Cluster_", rank, "_",
                                       nfeat, "_pseq.rds"))
saveRDS(object = PR2085, file = paste0("extdata/", rank, "Intersect/Cluster/C",
                                       k, "/PR2085D0_Cluster_", rank, "_",
                                       nfeat, "_pseq.rds"))
