## Project: Entropy ##
## PC Work git directory: ~/git/BVMetaGenomics ##
## PC Work project directory: ~/git/BVMetaGenomics ##
## CESGA Work git direcotry: /home/ulc/co/dfe/git/Entropy/
## CESGA Work project directory: /mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data)

# Different directories according to the place of work
setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data') ## CESGA directory
# setwd("projects/Entropy/data") ## PC Directory

# Required libraries
library(phyloseq)
library(data.table)
library(rentrez)
library(taxonomizr)

path = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/"
prepareDatabase('accessionTaxa.sql')

# In/out paths
# CESGA
# in.path = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data"
# out.path = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data"

# Load data (OTU + clinical)
otu = read.delim2(paste0(path,"otutableRefSeq.txt" ), header = T, sep = '\t')
clin = read.delim2(paste0(path,"task-nugent-score.txt"), header = T, sep = '\t')

# Retain NCBI ID?s for each OTU.Keeping us with a vector containing only the access numbers.
otu$X.OTU.ID = as.vector(otu$X.OTU.ID)
splitted = strsplit(otu$X.OTU.ID, '_')
ncbi = list()
for (i in seq_along(splitted)) {
  ncbi[[i]] = paste(splitted[[i]][1], splitted[[i]][2], sep = '_')
}
ncbi = unlist(ncbi)

# Convert NCBI ID?s to taxid. Each accession number corresponds to a taxonomy id.
taxid = list()
for (j in seq_along(ncbi)) {
  res = entrez_search(db = "nucleotide", term = ncbi[j])
  esums = entrez_summary(db = "nucleotide", id = res$ids)
  print(ncbi[j])
  taxid[[j]] = extract_from_esummary(esums, "taxid")
  print(taxid[[j]])
}
taxid = unlist(taxid)

# Convert taxid to taxonomy table
taxa = list()
for (t in seq_along(taxid)) {
  tt = getTaxonomy(
    ids = taxid[t],
    sqlFile = "accessionTaxa.sql",
    desiredTaxa = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))
  
  taxa[[t]] = data.frame(
    Rank1 = paste0('k__',tt[1]),
    Rank2 = paste0('p__',tt[2]),
    Rank3 = paste0('c__',tt[3]),
    Rank4 = paste0('o__',tt[4]),
    Rank5 = paste0('f__',tt[5]),
    Rank6 = paste0('g__',tt[6]),
    Rank7 = paste0('s__',tt[7])
  )
}
taxTable = as.data.frame(rbindlist(taxa))
rownames(taxTable) = ncbi

# Convert NA?s to 'g__'
taxTable$Rank1 = replace(taxTable$Rank1, taxTable$Rank1 == 'k__NA', 'k__')
taxTable$Rank2 = replace(taxTable$Rank2, taxTable$Rank2 == 'p__NA', 'p__')
taxTable$Rank3 = replace(taxTable$Rank3, taxTable$Rank3 == 'c__NA', 'c__')
taxTable$Rank4 = replace(taxTable$Rank4, taxTable$Rank4 == 'o__NA', 'o__')
taxTable$Rank5 = replace(taxTable$Rank5, taxTable$Rank5 == 'f__NA', 'f__')
taxTable$Rank6 = replace(taxTable$Rank6, taxTable$Rank6 == 'g__NA', 'g__')
taxTable$Rank7 = replace(taxTable$Rank7, taxTable$Rank6 == 's__NA', 's__')

saveRDS(taxTable, file = paste0(path,"TaxonomyTableRavel.rds" ))

