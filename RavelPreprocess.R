# 1. Paths and data loading
path = "projects/Entropy/data/"
otu = read.delim2(paste0(path,"otutableRefSeq.txt" ), header = T, sep = '\t')
clin = read.delim2(paste0(path,"task-nugent-score.txt"), header = T, sep = '\t')
taxa = readRDS(paste0(path,"TaxonomyTableRavel.rds"))

# 2. Functions Load
source("git/Entropy/functions/FunctionsDataFilter.R")
source("git/Entropy/functions/FunctionsGetSplitData.R")
source("git/Entropy/functions/FunctionsFeatSel.R")
source("git/Entropy/functions/FunctionsML.R")

# 3. Filter Phyloseq
## Get phy object and aglomerate by desired Rank 
Ravel_phyloseq = get.phylo.Ravel(otu = otu, clin = clin, path = path, taxonomyTable = taxa)
Ravel_Genus_C = phy.aglomerate(phyobject = Ravel_phyloseq, rank = "Rank6")

## Here we can choose whether to continue with the counts or their relative abundance. AR as comment
# Ravel_Genus_AR = relat.abun(Ravel_Genus_C) 
  
Ravel_Genus_C = prune.OTUs(phyobject = Ravel_Genus_C, Rank = "Rank6", pctg = 0.05, count = 0, vari = 0 )
# Ravel_Genus_AR = prune.OTUs(phyobject = Ravel_Genus_AR, Rank = "Rank6", pctg = 0.05, count = 0, vari = 0 )

Ravel_Genus_C = get.dataset(Ravel_Genus_C)
#Ravel_Genus_AR = get.dataset(Ravel_Genus_AR)


# 4. Data Processing 

## Split in train/test
Ravel_Genus_C_train_test = split.data(data = Ravel_Genus_C, seed = 123, pctg = 0.90,
                        path = "projects/Entropy/data/",project = "Ravel", id = "Genus_C")# project = Ravel or Sriniv. id = Rank and type of data

Ravel_Genus_C_train = Ravel_Genus_C_train_test[[1]]
Ravel_Genus_C_test = Ravel_Genus_C_train_test[[2]]
table(Ravel_Genus_C_test$target)


## Data normalization (Remember do it when use test)
Ravel_Genus_C_train = norm.dataset(Ravel_Genus_C_train)

# 5. FS
## Cluster(With counts only)
Ravel_Genus_C_train.CLUST = CLUST.FS(Ravel_Genus_C_train)
## Kruskal Wallis
Ravel_Genus_C_train.KW = KW.FS(Ravel_Genus_C_train, fs.type = "kruskal.test", pctge = 0.25)
## FCBF
Ravel_Genus_C_train.FCBF = FCBF.FS(Ravel_Genus_C_train, thold =0.005 )
## LDM
targets = as.factor(Ravel_Genus_C_train$target)
cols = sapply(Ravel_Genus_C_train, is.numeric)
variables = Ravel_Genus_C_train[cols]
Ravel_Genus_C_train.LDM = LDM.FS(data = Ravel_Genus_C_train,variables = variables,targets = targets, seed = 123)

## Save dataframes
fs = list(Ravel_Genus_C_train.CLUST, Ravel_Genus_C_train.KW, Ravel_Genus_C_train.FCBF,Ravel_Genus_C_train.LDM)
names(fs)= c("CLUST", "KW", "FCBF", "LDM")

for (i in seq_along(fs)) {
  saveRDS(fs[[i]], file = paste0("projects/Entropy/data/train/Ravel_Genus_C_train","_",names(fs[i]),"_",ncol(fs[[i]])-1,".rds"))
}


