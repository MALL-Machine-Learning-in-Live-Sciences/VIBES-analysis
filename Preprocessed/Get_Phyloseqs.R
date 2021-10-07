#Ravel
library(dplyr)
library(tibble)
library(readxl)
TaxonomyTableRavel <- read_excel("~/git/BVMetaGenomics/data/TaxonomyTableRavel.xlsx")
TaxonomyTableRavel = data.frame(TaxonomyTableRavel)
library(tidyverse)
TaxonomyTableRavel = TaxonomyTableRavel %>% remove_rownames %>% column_to_rownames(var="...1")
View(TaxonomyTableRavel)
otu = read.delim2(paste0("~/git/BVMetaGenomics/data/Ravel_otutable.txt"), header = T, sep = '\t')
clin = read.delim2(paste0("~/git/BVMetaGenomics/data/Ravel_Metadata.txt"), header = T, sep = '\t')
get.phylo.Ravel = function(otu, clin, taxonomyTable){
  require(phyloseq)
  require(dplyr)
  ## Retain only the access numbers
  otu$X.OTU.ID = as.vector(otu$X.OTU.ID)
  splitted = strsplit(otu$X.OTU.ID, '_')
  ncbi = list()
  for (i in seq_along(splitted)) {
    ncbi[[i]] = paste(splitted[[i]][1], splitted[[i]][2], sep = '_')
  }
  ncbi = unlist(ncbi)
  otu$X.OTU.ID = ncbi
  clinics = data_frame(clin)
  otu = data_frame(otu)
  ## Phyloseq format
  clinics <- clinics %>%
    tibble::column_to_rownames("X.SampleID") 
  
  otu <- otu %>%
    tibble::column_to_rownames("X.OTU.ID") 
  
  tax_mat <- as.matrix(taxonomyTable)
  otu_mat <- as.matrix(otu)
  
  # Make the main phyloseq object
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(clinics)
  BV_phyloseq <- phyloseq(OTU, TAX, samples)
  return(BV_phyloseq)
}
Ravel_phyloseq = get.phylo.Ravel(otu = otu, clin = clin, taxonomyTable = TaxonomyTableRavel)
tax_table(Ravel_phyloseq)
otu_table(Ravel_phyloseq)
sample_data(Ravel_phyloseq)

#Srinivasan
library(readxl)
TaxonomyTableSriniv <- read_excel("Entropy2/Data/Origin/TaxonomyTableSriniv.xlsx")
TaxonomyTableSriniv = data.frame(TaxonomyTableSriniv)
TaxonomyTableSriniv = TaxonomyTableSriniv %>% remove_rownames %>% column_to_rownames(var="...1")
View(TaxonomyTableSriniv)
otu2 = read.delim2(paste0("Entropy2/Data/Origin/Srinivasan_otutable.txt"), header = T, sep = '\t')
clin2 = read.delim2(paste0("Entropy2/Data/Origin/Sriniv_Metadata.txt"), header = T, sep = '\t')
get.phylo.Sriniv = function(otu, clin, taxonomyTable,target,path){
  require(phyloseq)
  require(dplyr)
  ## Categorize in: low, intermediate, high 
  clin$samplename = paste0("X",clin$samplename)
  clin$nugent = replace(clin$nugent, clin$nugent > 6, "high")
  clin$nugent = replace(clin$nugent, clin$nugent > 3 & clin$nugent < 7, "intermediate")
  clin$nugent = replace(clin$nugent, clin$nugent < 4, "low")
  table(clin$nugent)
  clinics = data_frame(clin)
  otu = data_frame(otu)
  
  ## Phyloseq format
  clinics <- clinics %>%
    tibble::column_to_rownames("samplename")
  
  otu <- otu %>%
    tibble::column_to_rownames("X.OTU.ID") 
  
  tax_mat <- as.matrix(taxonomyTable)
  otu_mat <- as.matrix(otu)
  
  # Make the main phyloseq object
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = tax_table(tax_mat)
  samples = sample_data(clinics)
  BV_phyloseq <- phyloseq(OTU, TAX, samples)
  
  return(BV_phyloseq)
}
Sriniv_phyloseq = get.phylo.Sriniv(otu = otu2, clin = clin2, taxonomyTable = TaxonomyTableSriniv)
tax_table(Sriniv_phyloseq)
otu_table(Sriniv_phyloseq)
sample_data(Sriniv_phyloseq)


saveRDS(Ravel_phyloseq, file = "~/git/BVMetaGenomics/data/Phyloseqs/Phyloseqs/Ravel_phyloseq.rds")
saveRDS(Sriniv_phyloseq, file = "~/git/BVMetaGenomics/data/Phyloseqs/Sriniv_phyloseq.rds")


# Dataframe Intersect All

phy.aglomerate = function(phyobject, rank){
  require(phyloseq)
  phy =  tax_glom(phyobject, rank)
  return(phy)
}
relat.abun = function(phyobject){
  require(phyloseq)
  phy = transform_sample_counts(phyobject, function(x) x / sum(x))
  return(phy)
}
## Counts
# We do the "aglomerate" in Ravel because different accession numbers may refer to the same species (but different strain).
RP = phy.aglomerate(phyobject = Ravel_phyloseq, rank = "Rank7") 
SP = phy.aglomerate(phyobject = Sriniv_phyloseq, rank = "Rank7")

#Save originals datasets
Ravel_OG = get.dataset(phyobject = RP)
Sriniv_OG = get.dataset(phyobject = SP)
saveRDS(object = Ravel_OG, file = "~/git/BVMetaGenomics/data/Original_Datasets/Ravel_C.rds")
saveRDS(object = Sriniv_OG, file = "~/git/BVMetaGenomics/data/Original_Datasets/Sriniv_C.rds")
Ravel_RA = relat.abun(phyobject = RP)
Sriniv_RA = relat.abun(phyobject = SP)
Ravel_RA = get.dataset(phyobject = Ravel_RA)
Sriniv_RA = get.dataset(phyobject = Sriniv_RA)
saveRDS(object = Ravel_RA, file = "~/git/BVMetaGenomics/data/Original_Datasets/Ravel_RA.rds")
saveRDS(object = Sriniv_RA, file = "~/git/BVMetaGenomics/data/Original_Datasets/Sriniv_RA.rds")
