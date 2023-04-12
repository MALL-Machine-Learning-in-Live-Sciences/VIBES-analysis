packgs <- c("phyloseq")
lapply(packgs, require, character.only = TRUE)
setwd("~/git/BV_Microbiome/")
#setwd("/mnt/netapp2/Store_uni/home/ulc/co/dfe/git/BV_Microbiome")
t1 = Sys.time()

# 0.Declare params and functions
rank <-  "Species" # Genus or "Species"
nfeat <- "24" # 24 and 14 for Species or 37 and 21 for Genus
nclust <- 3 # 3 or 4 
project <- "Ravel" # Ravel, Sriniv, PRJNA7977, PRJNA3020
co_abundances <- function(pseq) {
  require(igraph)
  require(SpiecEasi)
  require(phyloseq)
  pargs <- list(seed = 666, ncores = 4, rep.num = 50)
  spiec <- spiec.easi(pseq,
                      method = "glasso",
                      sel.criterion = "stars",
                      lambda.min.ratio = 1e-2,
                      nlambda = 20,
                      pulsar.params = pargs)
  print("listo1")
  cor <- cov2cor(apply(getOptCov(spiec), 2, as.numeric))
  weighted_adj_mat <- abs(cor) * getRefit(spiec)
  grph <- adj2igraph(weighted_adj_mat,
                     vertex.attr=list(name=taxa_names(pseq)))
  
  #Remove edges with very low weight 
  weight_threshold <- 0.01
  grph <- delete.edges(grph, which(abs(E(grph)$weight) < weight_threshold))
  
  #Remove unconnected vertices
  grph <- delete.vertices(grph, which(degree(grph) < 1))
  
  grph_deg <- degree(grph, v = V(grph), mode = "all")
  fine <- 500 # this will adjust the resolving power.
  print("listo2")
  #Clustering
  grph_louvain <- cluster_louvain(grph,
                                  weights = E(grph)$weight,
                                  resolution = 4)
  
  V(grph)$cluster <- grph_louvain$membership
  vertex_attr(grph, index = V(grph))
  print("listo3")
  nodes <- V(grph)$name
  taxa <- tax_table(pseq)
  taxa <- as.data.frame(taxa[nodes,])
  taxa$clusters <- V(grph)$cluster
  print("listo4")
  n <- plot_network(g = grph, physeq = pseq, type ='taxa',
               color = "Species", point_size = 6,
               title = paste(project,"Network with",ntaxa(pseq),"taxa for cluster", cls[i]))
  plot(n)
  l <- list(pseq, taxa, grph, n)
  names(l) <- c("filt_pseq", "taxa", "graph", "network")
  return(l)
}

# 1.Read phyloseqs
pseq_clusters <- readRDS(paste0("extdata/", rank, "Intersect/Cluster/C", nclust,
                                "/",project,"_Cluster_", rank, "_", nfeat,
                                "_pseq.rds"))

pseq <- readRDS(paste0("extdata/Phyloseqs/processed/", project,"_", rank, "_pseq.rds"))
cls <- unique(sample_data(pseq_clusters)$cluster)

# 2.Bucle for cluster class and run coab
list_cl <- list()
for (i in seq_along(cls)) {
  pseq_clusters_samples <- sample_names(subset_samples(pseq_clusters, cluster == cls[i]))
  sub_pseq <- subset_samples(pseq,(sample_names(pseq) %in% pseq_clusters_samples))
  sub_pseq <- tax_glom(physeq = sub_pseq, taxrank = rank)
  pseq_fil <- filter_taxa(sub_pseq, function(x) sum(x > 2) > (0.05*length(x)), TRUE)
  list_cl[[i]] <- co_abundances(pseq = pseq_fil)
}
names(list_cl) <- cls

# 3. Save list with coab for each cluster class
saveRDS(object = list_cl, file = paste0("extdata/", rank, "Intersect/Cluster/C",
                                      nclust,"/", project,"_Cluster_", rank,
                                      "_", nfeat, "_CoAb.rds"))
t2 = Sys.time()
runtime = t2-t1
print(runtime)