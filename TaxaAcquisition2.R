# Taxa Acquisition 
setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data') 
path = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/"
#path = "projects/Entropy/data/"
otu = read.delim2(paste0(path,"otutable-Amsel.txt" ), header = T, sep = '\t')
clin = read.delim2(paste0(path,"task-amsel.txt"), header = T, sep = '\t')
# Required libraries
library(phyloseq)
library(data.table)
library(rentrez)
library(taxonomizr)
library(taxize)
prepareDatabase('accessionTaxa.sql')

# Retain names for each OTU.Keeping us with a vector containing only the name.
ncbi = as.vector(otu$X.OTU.ID)
taxid = list()
for (j in seq_along(ncbi)) {
  res = get_ids(ncbi[j],db = "ncbi")
  #print(ncbi[j])
  taxid[[j]] = res$ncbi[[1]]
  #print(taxid[[j]])
}
taxid = unlist(taxid)
# We have search manually taxid for corrupt/uncomplete OTUs 
kk = cbind(otu, taxid)
# Extract names from OTUs that have NA tax id
a = subset(kk,is.na(taxid))
# Extract index from data
index = rownames(a)
index = as.numeric(index)
# Make by hand a vector  with taxids(we coudnt find 2 of them)
vector = as.character(c("699240", "84111", NA, "838", "39948",
                        "884684", "165779", "838", "1301", NA,
                        "836", "838","1350","1301", "838", "29465",
                        "2129", "201174"))
# Replace NAs with taxIDs
taxid = replace(x = taxid, list = index, values = vector)
kk$taxid = taxid


# Convert taxid to taxonomy tables
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

saveRDS(taxTable, file = paste0(path,"TaxonomyTableSriniv.rds" ))





