# Validation
# Enfrentar los datos del segundo articulo al modelo sacado con los del primero
# Load phyloseq
Ravel = readRDS("projects/Entropy/data/Ravel_phyloseq.rds")
Ravel = phy.aglomerate(phyobject = Ravel, rank = "Rank6")
Sriniv = readRDS("projects/Entropy/data/Sriniv_Nugent_phyloseq.rds")
Sriniv = phy.aglomerate(phyobject = Sriniv, rank = "Rank6")

#Load model
model= readRDS("projects/Entropy/data/models/RF_8.rds")
# Extract features
names = model$features
# Extract corespondence between OTU anf Genus
Ravel = as.data.frame(tax_table(Ravel))
features = Ravel[names,]
features = cbind(rownames(features),features$Rank6)
vg = features[,2]
vn = features[,1]

#Extract a subset from second phyloseq with only features from the model
Sriniv_sub <- subset_taxa(Sriniv, Rank6 %in% vg) # We have 7 from 8 features
table = as.data.frame(tax_table(Sriniv_sub))
table = cbind(rownames(table), table$Rank6)
k = table[,2]
Sriniv_sub = get.dataset(phyobject = Sriniv_sub)
cols <- sapply(Sriniv_sub, is.numeric)
variables = Sriniv_sub[cols]
targets = Sriniv_sub$target
colnames(variables) = k
dataset = as.data.frame(variables)
dataset = as.data.frame(t(dataset))
features = as.data.frame(features)
library(dplyr)
df <- tibble::rownames_to_column(dataset, "V2")

merge(x = df, y = features, by = "V2", all = TRUE)
dataset = as.data.frame(cbind(dataset, target =targets ))

