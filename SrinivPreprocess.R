# 1. Paths and data loading
path = "projects/Entropy/data/"
otu = read.delim2(paste0(path,"otutable-Amsel.txt" ), header = T, sep = '\t')
clin = read.delim2(paste0(path,"task-amsel.txt"), header = T, sep = '\t')
taxa = readRDS(paste0(path,"TaxonomyTableSriniv.rds"))

# 2. Functions Load
source("git/Entropy/functions/FunctionsDataFilter.R")
source("git/Entropy/functions/FunctionsGetSplitData.R")
source("git/Entropy/functions/FunctionsFeatSel.R")

#Name genus as otu cause we do not want to lose information
taxa$Rank6[7] = "g__BVAB1"
taxa$Rank6[14] = "g__BVAB2"
taxa$Rank6[32] = "g__BVAB3"
taxa$Rank6[40] = "g__TM7"

Sriniv_phyloseq = get.phylo.Sriniv(otu = otu, clin = clin,taxonomyTable = taxa,
                               target = "Amsel",path = path)

Sriniv_Genus_C = phy.aglomerate(phyobject = Sriniv_phyloseq, rank = "Rank6")

Sriniv_Genus_C = relat.abun(Sriniv_Genus_C) 

Sriniv_Genus_C = prune.OTUs(phyobject = Sriniv_Genus_C, Rank = "Rank6", pctg = 0.05, count = 0, vari = 0 )

Sriniv_Genus_C = get.dataset(Sriniv_Genus_C)

colnames(Sriniv_Genus_C) <- make.names(colnames(Sriniv_Genus_C), unique=TRUE)

Sriniv_Genus_C_train_test = split.data(data = Sriniv_Genus_C, seed = 123, pctg = 0.90,
                                      path = "projects/Entropy/data/",project = "Sriniv_Amsel", id = "Genus_AR")# project = Ravel or Sriniv. id = Rank and type of data

Sriniv_Genus_C_train = Sriniv_Genus_C_train_test[[1]]
Sriniv_Genus_C_test = Sriniv_Genus_C_train_test[[2]]
table(Sriniv_Genus_C_test$target)

Sriniv_Genus_C_train = norm.dataset(Sriniv_Genus_C_train)

# 5. FS
Sriniv_Genus_C_train.CLUST = CLUST.FS(Sriniv_Genus_C_train)
## Kruskal Wallis
Sriniv_Genus_C_train.KW = KW.FS(Sriniv_Genus_C_train, fs.type = "kruskal.test", pctge = 0.30)
## FCBF
Sriniv_Genus_C_train.FCBF = FCBF.FS(Sriniv_Genus_C_train, thold =0.0025 )
## LDM
targets = as.factor(Sriniv_Genus_C_train$target)
cols = sapply(Sriniv_Genus_C_train, is.numeric)
variables = Sriniv_Genus_C_train[cols]
Sriniv_Genus_C_train.LDM = LDM.FS(data = Sriniv_Genus_C_train,variables = variables,targets = targets, seed = 123)


## Save dataframes
fs = list(Sriniv_Genus_C_train.CLUST, Sriniv_Genus_C_train.KW, Sriniv_Genus_C_train.FCBF,Sriniv_Genus_C_train.LDM)
names(fs)= c("CLUST", "KW", "FCBF", "LDM")

for (i in seq_along(fs)) {
  saveRDS(fs[[i]], file = paste0("projects/Entropy/data/train/Sriniv_Amsel_Genus_AR_train","_",names(fs[i]),"_",ncol(fs[[i]])-1,".rds"))
}

