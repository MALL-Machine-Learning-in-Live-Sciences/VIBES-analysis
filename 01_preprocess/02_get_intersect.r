setwd("~/git/BV_Microbiome/")
require(phyloseq)
rank <- "Species" # "Genus" or "Species"
Ravel <- readRDS("extdata/Phyloseqs/Ravel_phyloseq.rds")
Sriniv <- readRDS("extdata/Phyloseqs/Sriniv_phyloseq.rds")
PRJNA3020 <- readRDS("extdata/Phyloseqs/PRJNA3020_phyloseq.rds")
PRJNA7977 <- readRDS("extdata/Phyloseqs/PRJNA7977_phyloseq.rds")
PRJNA2085 <- readRDS("extdata/Phyloseqs/PRJNA2085_phyloseq.rds")

# 1.Rename tax names for Ravel and Sriniv
taxnames <- colnames(tax_table(PRJNA3020))
colnames(tax_table(Ravel)) <- taxnames
colnames(tax_table(Sriniv)) <- taxnames
# 2.1.Preprocess PRJNA3020 to retain patients without treatment
#meta <- data.frame(sample_data(PRJNA3020))
#d0 <- meta[grep("D0", meta$sample_alias), ]
#PRJNA3020 <- subset_samples(PRJNA3020,
#                            (sample_names(PRJNA3020) %in% rownames(d0)))
# 2.2.Preprocess PRJNA7977 to retain first sample of each patient
#meta2 <- data.frame(sample_data(PRJNA7977))
#meta2 <- meta2[order(meta2$sample_alias, decreasing = TRUE),]
#meta2 <- meta2[-grep("W10", meta2$sample_alias), ]
#meta2$sample_alias <- sapply(strsplit(basename(meta2$sample_alias), "_"), `[`, 1)
#meta2 <- meta2[!rev(duplicated(rev(meta2$sample_alias))),]
#PRJNA7977 <- subset_samples(PRJNA7977,
#                            (sample_names(PRJNA7977) %in% rownames(meta2)))
# 2.2.Preprocess PRJNA7977 to retain first sample of each patient
#meta3 <- data.frame(sample_data(PRJNA2085))
#meta3 <- meta3[order(meta3$sample_title, decreasing = TRUE),]
#meta3$sample_title <- sapply(strsplit(basename(meta3$sample_title), "_"), `[`, 1)
#meta3 <- meta3[!rev(duplicated(rev(meta3$sample_title))),]
#PRJNA2085 <- subset_samples(PRJNA2085,
#                            sample_names(PRJNA2085) %in% rownames(meta3))

# 3.Standardisation names of Ravel and Sriniv
if (rank == "Genus") {
  Ravel.df = as.data.frame(tax_table(Ravel))
  Sriniv.df = as.data.frame(tax_table(Sriniv))
  # Rename Fannyhessea vaginae for Atopobium vaginae, same speceis 
  Ravel.df["NR_117757.1",] <-c("k__Bacteria", "p__Actinobacteria", "c__Coriobacteriia",
                               "o__Coriobacteriales", "f__Atopobiaceae", "g__Atopobium",
                               "s__Atopobium vaginae")
  Sriniv.df["Atopobium vaginae",] <-c("k__Bacteria", "p__Actinobacteria", "c__Coriobacteriia",
                                      "o__Coriobacteriales", "f__Atopobiaceae", "g__Atopobium",
                                      "s__Atopobium vaginae")
  Ravel.df$Genus  <- substr(Ravel.df$Genus , start = 4, stop = 100)
  Sriniv.df$Genus <- substr(Sriniv.df$Genus, start = 4, stop = 100)
  require(stringr)
  Ravel.df$Genus <- str_replace_all(string = Ravel.df$Genus, "[ ]", "_" )
  Sriniv.df$Genus <- str_replace_all(string = Sriniv.df$Genus , "[ ]", "_" )
  tax_table(Ravel) <- as.matrix(Ravel.df)
  tax_table(Sriniv) <- as.matrix(Sriniv.df)
}else if(rank == "Species"){
  Ravel.df = as.data.frame(tax_table(Ravel))
  Sriniv.df = as.data.frame(tax_table(Sriniv))
  #Change also genus name
  # Rename Fannyhessea vaginae for Atopobium vaginae, same specie 
  Ravel.df["NR_117757.1",] <-c("k__Bacteria", "p__Actinobacteria", "c__Coriobacteriia",
                               "o__Coriobacteriales", "f__Atopobiaceae", "g__Atopobium",
                               "s__Atopobium vaginae")
  Sriniv.df["Atopobium vaginae",] <-c("k__Bacteria", "p__Actinobacteria", "c__Coriobacteriia",
                                      "o__Coriobacteriales", "f__Atopobiaceae", "g__Atopobium",
                                      "s__Atopobium vaginae")
  Ravel.df$Genus  <- substr(Ravel.df$Genus , start = 4, stop = 100)
  Sriniv.df$Genus <- substr(Sriniv.df$Genus, start = 4, stop = 100)
  require(stringr)
  Ravel.df$Genus <- str_replace_all(string = Ravel.df$Genus, "[ ]", "_" )
  Sriniv.df$Genus <- str_replace_all(string = Sriniv.df$Genus , "[ ]", "_" )
  # Rename Fannyhessea for Atopobium 
  Ravel.df$Species  <- substr(Ravel.df$Species , start = 4, stop = 100)
  Sriniv.df$Species <- substr(Sriniv.df$Species, start = 4, stop = 100)
  Ravel.df$Species <- str_replace_all(string = Ravel.df$Species, "[ ]", "_" )
  Sriniv.df$Species <- str_replace_all(string = Sriniv.df$Species , "[ ]", "_" )
  tax_table(Ravel) <- as.matrix(Ravel.df)
  tax_table(Sriniv) <- as.matrix(Sriniv.df)
}
#Save psqs with tax and subsa,ples preprocessed for futurte analysis
saveRDS(object = Ravel, file = paste0("extdata/Phyloseqs/processed/Ravel_", rank, "_pseq.rds"))
saveRDS(object = Sriniv, file = paste0("extdata/Phyloseqs/processed/Sriniv_", rank, "_pseq.rds"))
saveRDS(object = PRJNA7977, file = paste0("extdata/Phyloseqs/processed/PRJNA7977_", rank, "_pseq.rds"))
saveRDS(object = PRJNA3020, file = paste0("extdata/Phyloseqs/processed/PRJNA3020_", rank, "_pseq.rds"))
saveRDS(object = PRJNA2085, file = paste0("extdata/Phyloseqs/processed/PRJNA2085_", rank, "_pseq.rds"))

#Funciones
phy.aglomerate <- function(phyobject, rank){
  require(phyloseq)
  phy =  tax_glom(phyobject, rank)
  return(phy)
}
get.dataset <- function(phyobject){
  require(phyloseq)
  clinics = as.data.frame(sample_data(phyobject))
  otus = as.data.frame(t(get_taxa(phyobject)))
  identical(rownames(clinics), rownames(otus))
  data = cbind(otus, clinics)
  return(data)
}
prune.OTUs <- function(phyobject, Rank, pctg = 0.05, count = 0, vari = 0){
  require(phyloseq)
  if (Rank == "Genus"){
    #Delete unidentified Ranks
    ps <- subset_taxa(phyobject, !Genus %in% c("g__"))
    
    #Make and apply the filter
    filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > count), 
                                          A = pctg*nsamples(ps))
    phy <- prune_taxa(filter, ps)
    
    # If we want filter by variance 
    phy = filter_taxa(phy, function(x) var(x) > vari, TRUE)
    
  }else if(Rank == "Species"){
    #Delete unidentified Ranks
    ps <- subset_taxa(phyobject, !Species %in% c("s__"))
    
    #Make and apply the filter
    filter <- phyloseq::genefilter_sample(ps, filterfun_sample(function(x) x > count), 
                                          A = pctg*nsamples(ps))
    phy <- prune_taxa(filter, ps)
    
    # If we want filter by variance 
    phy = filter_taxa(phy, function(x) var(x) > vari, TRUE)
  }else{
    print("Introduce valid rank (Genus or Species)")
  }
  return(phy)
}

# 4.Aglomerate by Genus/Species
Ravel <- phy.aglomerate(phyobject = Ravel, rank = rank)
Sriniv <- phy.aglomerate(phyobject = Sriniv, rank = rank)
PRJNA3020 <- phy.aglomerate(phyobject = PRJNA3020, rank = rank)
PRJNA7977 <- phy.aglomerate(phyobject = PRJNA7977, rank = rank)
PRJNA2085 <- phy.aglomerate(phyobject = PRJNA2085, rank = rank)

# 5.Prune OTUs
#Ravel <- prune.OTUs(phyobject = Ravel, Rank = rank)
#Sriniv <- prune.OTUs(phyobject = Sriniv, Rank = rank)
#PRJNA3020 <- prune.OTUs(phyobject = PRJNA3020, Rank = rank)
#PRJNA7977 <- prune.OTUs(phyobject = PRJNA7977, Rank = rank)

# 6.Extract common taxs
if (rank == "Genus") {
  common_tax <- Reduce(intersect, list(
    data.frame(tax_table(Ravel))$Genus,
    data.frame(tax_table(Sriniv))$Genus,
    data.frame(tax_table(PRJNA3020))$Genus,
    data.frame(tax_table(PRJNA7977))$Genus,
    data.frame(tax_table(PRJNA2085))$Genus
  ))
}else if (rank == "Species") {
  common_tax <- Reduce(intersect, list(
    data.frame(tax_table(Ravel))$Species,
    data.frame(tax_table(Sriniv))$Species,
    data.frame(tax_table(PRJNA3020))$Species,
    data.frame(tax_table(PRJNA7977))$Species,
    data.frame(tax_table(PRJNA2085))$Species
  ))
}

require(microViz)
Ravel <- tax_select(ps = Ravel, tax_list = common_tax, ranks_searched = rank,
                    strict_matches = TRUE, n_typos = 1, deselect = FALSE)
Sriniv <- tax_select(ps = Sriniv, tax_list = common_tax, ranks_searched = rank,
                    strict_matches = TRUE, n_typos = 1, deselect = FALSE)
PRJNA3020 <- tax_select(ps = PRJNA3020, tax_list = common_tax, ranks_searched = rank,
                    strict_matches = TRUE, n_typos = 1, deselect = FALSE)
PRJNA7977 <- tax_select(ps = PRJNA7977, tax_list = common_tax, ranks_searched = rank,
                    strict_matches = TRUE, n_typos = 1, deselect = FALSE)
PRJNA2085 <- tax_select(ps = PRJNA2085, tax_list = common_tax, ranks_searched = rank,
                        strict_matches = TRUE, n_typos = 1, deselect = FALSE)

# 7.Ordering phyloseqs and reanme OTUs
# Return otu_table with OTUs on columns
Ravel <- tax_sort(data = Ravel, by = "name", at = rank, tree_warn = TRUE,
                  verbose = TRUE, trans = "identity", use_counts = TRUE)
Sriniv <- tax_sort(data = Sriniv, by = "name", at = rank, tree_warn = TRUE,
                  verbose = TRUE, trans = "identity", use_counts = TRUE)
PRJNA3020 <- tax_sort(data = PRJNA3020, by = "name", at = rank, tree_warn = TRUE,
                  verbose = TRUE, trans = "identity", use_counts = TRUE)
PRJNA7977 <- tax_sort(data = PRJNA7977, by = "name", at = rank, tree_warn = TRUE,
                  verbose = TRUE, trans = "identity", use_counts = TRUE)
PRJNA2085 <- tax_sort(data = PRJNA2085, by = "name", at = rank, tree_warn = TRUE,
                      verbose = TRUE, trans = "identity", use_counts = TRUE)
# Rename for same ID
if (rank == "Genus") {
  t_names <- data.frame(tax_table(Ravel))$Genus
}else if(rank == "Species"){
  t_names <- data.frame(tax_table(Ravel))$Species
}
taxa_names(Ravel) <- t_names
taxa_names(Sriniv) <- t_names
taxa_names(PRJNA3020) <- t_names
taxa_names(PRJNA7977) <- t_names
taxa_names(PRJNA2085) <- t_names

# 8.Saving phyloseqs
saveRDS(object = Ravel, file = paste0("extdata/",rank,"Intersect/Ravel_", rank,
                                      "_pseq_", length(common_tax), ".rds"))
saveRDS(object = Sriniv, file = paste0("extdata/",rank,"Intersect/Sriniv_",rank,
                                       "_pseq_", length(common_tax), ".rds"))
saveRDS(object = PRJNA3020, file = paste0("extdata/",rank,"Intersect/PRJNA3020_",
                                          rank, "_pseq_", length(common_tax), ".rds"))
saveRDS(object = PRJNA7977, file = paste0("extdata/",rank,"Intersect/PRJNA7977_",
                                          rank, "_pseq_", length(common_tax), ".rds"))
saveRDS(object = PRJNA2085, file = paste0("extdata/",rank,"Intersect/PRJNA2085_",
                                          rank, "_pseq_", length(common_tax), ".rds"))
