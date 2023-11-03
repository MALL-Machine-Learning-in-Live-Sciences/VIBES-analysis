##### Check indexes #####
# 0.Load packages
packgs <- c("phyloseq", "factoextra", "NbClust")
lapply(packgs, require, character.only = TRUE)

# 1.Declare several variables to perform analysis on several configurations
rank <-  "Species" 
nfeat <- "22"
trans <- "clr"

# 2.2.Normalization
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

# 2.Load data
Ravel <- readRDS(paste0("00_preprocess_cohorts/data/", rank, "Intersect/Ravel_",
                        rank, "_pseq_", nfeat, ".rds"))
Ravel <- norm_dataset(pseq = Ravel)
ravel_mat <- as.data.frame(otu_table(Ravel)@.Data)

# Elbow and shilouette
silhouette_plot <- fviz_nbclust(ravel_mat, kmeans, method = "silhouette") + 
  labs(subtitle = "Silhouette method") # 2,4
silhouette_plot$layers[3] <- NULL
silhouette_plot


elbow_plot <- fviz_nbclust(ravel_mat, kmeans, method = "wss") + 
  labs(subtitle = "Elbow method") # 4
elbow_plot


# Bucle with each index

vector = c("kl", "ch", "hartigan", "cindex", "db", "duda", "pseudot2",
           "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain",
           "gamma", "gplus", "tau", "dunn", "sdindex", "sdbw")
set.seed(1234)
res = list()
for (i in seq_along(vector)) {
  print(vector[i])
  res[[i]] = NbClust(data = ravel_mat, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 10, 
                   method = "kmeans", index = vector[i])
  print("DONE!")
}
names(res) <- vector
