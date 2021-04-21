# Taxa Acquisition 
"projects/Entropy/data/"
path = "projects/Entropy/data/task-amsel.txt"
otu = read.delim2(paste0(path,"otutable-Amsel.txt" ), header = T, sep = '\t')
clin = read.delim2(paste0(path,"task-amsel.txt"), header = T, sep = '\t')

x = c("phyloseq", "data.table", "rentrez", "dplyr", "tibble")
lapply(x, require, character.only = TRUE)

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
return(taxid)







