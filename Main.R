path = "projects/Entropy/data/"
otu = read.delim2(paste0(path,"otutableRefSeq.txt" ), header = T, sep = '\t')
clin = read.delim2(paste0(path,"task-nugent-score.txt"), header = T, sep = '\t')
taxa = readRDS("projects/Entropy/data/TaxonomyTable.rds")

source("git/Entropy/functions/FunctionsTaxaAcquisition.R")
source("git/Entropy/functions/FunctionsDataFilter.R")
source("git/Entropy/functions/FunctionsGetSplitData.R")
source("git/Entropy/functions/FunctionsFeatSel.R")


phyloseq = get.phylo(otu = otu, clin = clin, path = path, taxonomyTable = taxa, id = "BV")
genus_phy = phy.aglomerate(phyobject = phyloseq, rank = "Rank6")
pruned_genus_phy = prune.OTUs(phyobject = genus_phy, Rank = "Rank6", pctg = 0.05, count = 0, vari = 0 )
BV_data_Gen = get.dataset(pruned_genus_phy)
BV_Gen_Norm = norm.dataset(BV_data_Gen)
train_test = split.data(data = BV_Gen_Norm, seed = 123, pctg = 0.90,path = "projects/Entropy/data/",project = "BV", id = "Genus")

train = train_test[[1]]
test = train_test[[2]]
table(test$target)

# Kruskal Wallis
train.KW = kruskal.FS(train, fs.type = "kruskal.test", nfeat = 20)
# FCBF
train.FCBF = FCBF.FS(train, thold =0.005 )
# LDM
targets = as.factor(train$target)
cols = sapply(train, is.numeric)
variables = train[cols]
train.LDM = LDM.FS(data = train,variables = variables,targets = targets, seed = 123, )




