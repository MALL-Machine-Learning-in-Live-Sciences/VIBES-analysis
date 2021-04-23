# 1. Paths and data loading
path = "projects/Entropy/data/"
otu = read.delim2(paste0(path,"otutable-Amsel.txt" ), header = T, sep = '\t')
clin = read.delim2(paste0(path,"task-amsel.txt"), header = T, sep = '\t')
taxa = readRDS(paste0(path,"TaxonomyTableSriniv.rds"))

# 2. Functions Load
source("git/Entropy/functions/FunctionsDataFilter.R")
source("git/Entropy/functions/FunctionsGetSplitData.R")
source("git/Entropy/functions/FunctionsFeatSel.R")
source("git/Entropy/functions/FunctionsML.R")

#Name genus as otu cause we do not want to lose information
taxa$Rank6[7] = "g__BVAB1"
taxa$Rank6[14] = "g__BVAB2"
taxa$Rank6[32] = "g__BVAB3"
taxa$Rank6[40] = "g__TM7"

Sriniv_phyloseq = get.phylo.Sriniv(otu = otu, clin = clin,taxonomyTable = taxa,
                               target = "Nugent",path = path)

Sriniv_Genus_C = phy.aglomerate(phyobject = Sriniv_phyloseq, rank = "Rank6")
Sriniv_Genus_C = prune.OTUs(phyobject = Sriniv_Genus_C, Rank = "Rank6", pctg = 0.05, count = 0, vari = 0 )
tax_table(Sriniv_Genus_C)
