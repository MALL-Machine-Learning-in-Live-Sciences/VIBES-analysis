# Validation
# Enfrentar los datos del segundo articulo al modelo sacado con los del primero
Ravel = readRDS("projects/Entropy/data/Ravel_phyloseq.rds")
Ravel = phy.aglomerate(phyobject = Ravel, rank = "Rank6")
Sriniv = readRDS("projects/Entropy/data/Sriniv_Nugent_phyloseq.rds")
Sriniv = phy.aglomerate(phyobject = Sriniv, rank = "Rank6")


model= readRDS("projects/Entropy/data/models/RF_8.rds")
names = model$features
Ravel = as.data.frame(tax_table(Ravel))
features = Ravel[names,]
features = cbind(rownames(features),features$Rank6)
vg = features[,2]
vn = features[,1]
Sriniv_sub <- subset_taxa(Sriniv, Rank6 %in% vg)
tax_table(Sriniv_sub)
Sriniv_sub = get.dataset(phyobject = Sriniv_sub)



