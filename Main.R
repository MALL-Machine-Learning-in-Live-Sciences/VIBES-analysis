path = "projects/Entropy/data/"
otu = read.delim2(paste0(path,"otutableRefSeq.txt" ), header = T, sep = '\t')
clin = read.delim2(paste0(path,"task-nugent-score.txt"), header = T, sep = '\t')
taxa = readRDS("projects/Entropy/data/TaxonomyTable.rds")

source("git/Entropy/functions/FunctionsTaxaAcquisition.R")
source("git/Entropy/functions/FunctionsDataFilter.R")
source("git/Entropy/functions/FunctionsGetSplitData.R")
source("git/Entropy/functions/FunctionsFeatSel.R")
source("git/Entropy/functions/FunctionsML.R")


phyloseq = get.phylo(otu = otu, clin = clin, path = path, taxonomyTable = taxa, id = "BV")
genus_phy_C = phy.aglomerate(phyobject = phyloseq, rank = "Rank6")

#SACAR ABUNDANCIA REALTIVA ANTES DE ELIMINAR OTUS
genus_phy_AR = relat.abun(genus_phy_C) 
  
pruned_genus_phy_C = prune.OTUs(phyobject = genus_phy_C, Rank = "Rank6", pctg = 0.05, count = 0, vari = 0 )
pruned_genus_phy_AR = prune.OTUs(phyobject = genus_phy_AR, Rank = "Rank6", pctg = 0.05, count = 0, vari = 0 )

BV_Gen_C = get.dataset(pruned_genus_phy_C)
BV_Gen_AR = get.dataset(pruned_genus_phy_AR)



#LA NORMALIZACION  HACERLA DESPUES ANTES DEL ML Y ESCALADO EN EL ML
train_test = split.data(data = BV_Gen_C, seed = 123, pctg = 0.90,path = "projects/Entropy/data/",project = "BV", id = "Genus")
train = train_test[[1]]
test = train_test[[2]]
table(test$target)

# Kruskal Wallis
train.KW = KW.FS(train, fs.type = "kruskal.test", nfeat = 20)
# FCBF
train.FCBF = FCBF.FS(train, thold =0.005 )
# LDM
targets = as.factor(train$target)
cols = sapply(train, is.numeric)
variables = train[cols]
train.LDM = LDM.FS(data = train,variables = variables,targets = targets, seed = 123)



bmr = ML.exec(dataset = train.FCBF)



