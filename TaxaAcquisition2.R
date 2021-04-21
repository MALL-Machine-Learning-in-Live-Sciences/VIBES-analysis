# Taxa Acquisition 
setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data') 
path = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/"
#path = "projects/Entropy/data/"
otu = read.delim2(paste0(path,"otutable-Amsel.txt" ), header = T, sep = '\t')
clin = read.delim2(paste0(path,"task-amsel.txt"), header = T, sep = '\t')
x = c("phyloseq", "data.table", "rentrez","taxnomizr", "dplyr", "tibble")
lapply(x, require, character.only = TRUE)
prepareDatabase('accessionTaxa.sql')


otu<-otu[!(otu$X.OTU.ID=="Prevotella genogroup 3" | otu$X.OTU.ID=="Prevotella genogroup 4" | 
              otu$X.OTU.ID=="Prevotella genogroup 7" | otu$X.OTU.ID=="Prevotella genogroup 6"|
              otu$X.OTU.ID=="Veillonella atypica/parvula"),]


# Retain names for each OTU.Keeping us with a vector containing only the name.
otu$X.OTU.ID = as.vector(otu$X.OTU.ID)
ncbi = as.vector(otu$X.OTU.ID)
# Convert Names to taxid. We hae to modify this fnct cause in Amstel have more than 1 id cause we havnt humber acces.
taxid = list()
for (j in seq_along(ncbi)) {
  res = entrez_search(db = "nucleotide", term = ncbi[j])
  esums = entrez_summary(db = "nucleotide", id = res$ids)
  print(ncbi[j])
  taxid[[j]] = extract_from_esummary(esums, "taxid")
  taxid[[j]] = taxid[[j]][[1]]
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

saveRDS(taxTable, file = paste0(path,"TaxonomyTableSriniv.rds" ))







