## Project: Entropy ##
## PC Work git directory: ~/git/Entropy ##
## PC Work project directory: ~/projects/Entropy ##
## CESGA Work git direcotry: /home/ulc/co/dfe/git/Entropy/
## CESGA Work project directory: /mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data)

# Different directories according to the place of work
setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data') ## CESGA directory
# setwd("projects/Entropy/data") ## PC Directory

# Required libraries
x = c("phyloseq", "data.table", "rentrez", "taxnomizr", "dplyr", "tibble")
lapply(x, require, character.only = TRUE)

path = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/"
project ="BV"
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

saveRDS(taxTable, file = paste0(path,"TaxonomyTable.rds" ))
taxonomyTable = taxTable
# Retain the access numbers, rename target, categorize variables and adecuate all for a phy object.
#path2 = "projects/Entropy/data/" #En el cesga usaria el paht del ppio
#taxonomyTable = readRDS(file = "projects/Entropy/data/TaxonomyTable.rds")
#otu = read.delim2(paste0(path2,"otutableRefSeq.txt" ), header = T, sep = '\t') #En el cesga usaria el paht del ppio
#clin = read.delim2(paste0(path2,"task-nugent-score.txt"), header = T, sep = '\t') #En el cesga usaria el paht del ppio

readRDS("projects/Entropy/data/TaxonomyTable.rds")

## Retain only the access numbers
otu$X.OTU.ID = as.vector(otu$X.OTU.ID)
splitted = strsplit(otu$X.OTU.ID, '_')
ncbi = list()
for (i in seq_along(splitted)) {
  ncbi[[i]] = paste(splitted[[i]][1], splitted[[i]][2], sep = '_')
}
ncbi = unlist(ncbi)
otu$X.OTU.ID = ncbi

## Categorize in: low, intermediate, high 
high = clin[which(clin$Var > 6),]
intermediate = clin[which(clin$Var >= 4 & clin$Var <= 6),]
low = clin[which(clin$Var < 4),]
high$Var = "High"
intermediate$Var = "Intermediate"
low$Var = "Low"
clinics = rbind(high,intermediate, low)
clinics = arrange(clinics, X.SampleID)
clin = data_frame(clinics)

## Phyloseq format
clin <- clin %>%
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
saveRDS(BV_phyloseq, file = paste0(path,project,"_phyloseq.rds"))

source("git/Entropy/Preprocess_Functions.R")

# Prueba Preprocessing phyloseq 0 -> FS
BV_phyloseq <- readRDS("projects/Entropy/data/BV_phyloseq.rds")

BV_Genus =  aglomerate(BV_phyloseq, rank = "Rank6")

BV_Genus_AR = relat.abun(BV_Genus)

Pruned_Genus_AR = prune.OTUs(phyobject = BV_Genus_AR, pctg = 0.05, count = 0, vari = 0)

Data_BV_Gn_AR = get.dataset(Pruned_Genus_AR)

Data_BV_Gn_AR_Norm = norm.dataset(Data_BV_Gn_AR)

train_test = split.data(data = Data_BV_Gn_AR_Norm, seed = 17, pctg = 0.90)

save.splits(path = "git/Entropy/", list = train_test, project = "BV_Genus", id = "AR")

train = readRDS(file = "projects/Entropy/data/train/BV_Genus_AR_train.rds")

#FS
# Kruskalwalis
Genus_train_KW = kruskal.FS(data = train, fs.type ='kruskal.test', nfeat = 15 )

#FCBF
Genus_train_FCBF = FCBF.FS(data = train, thold = 0.005)

#LDM(I cant implement a function due to incompatibilities with function ldm)
targets = as.factor(train$target)
cols = sapply(train, is.numeric)
variables = train[cols]
require(LDM)
#ExecuteLDM 
fit.ldm = ldm(variables ~target,
              data = train,
              test.otu = TRUE, 
              test.global = TRUE,
              dist.method = "bray",
              fdr.nominal = 0.1,
              n.perm.max = 10000,
              seed = 19)

#Filter features 
w1 = which(fit.ldm$q.otu.omni[1,] < 0.1)
(n.otu.omni.m1 = length(w1))
features = (otu.omni.m1 = colnames(fit.ldm$q.otu.omni)[w1])

#Select new features
filtered.data = variables[,features]
Genus_train_LDM = as.data.frame(cbind(filtered.data, target = targets))

fs = list(Genus_train_FCBF, Genus_train_KW, Genus_train_LDM)
names(fs)= c("FCBF", "KW", "LDM")

for (i in seq_along(fs)) {
  saveRDS(fs[[i]], file = paste0("projects/Entropy/data/train/Genus","_",names(fs[i]),"_",ncol(fs[[i]])-1,".rds"))
}

