#ML 
source("git/Entropy/functions/FunctionsML.R")
path = "projects/Entropy/data/train/"
pattern = "Ravel_Genus_AR"
l = list.files(path, pattern = pattern)
l = l[-5]
l = l[-1]

ldata = list()
for (i in 1:length(l)) {
  ldata[[i]] = readRDS(paste(path, l[[i]], sep=''))
}
names(ldata) = l
bmr = list()
#identificador = list("Counts_FCBF", "Counts_KW", "Counts_LDM")
#identificador = list("RA_Counts", "RA_KW", "RA_LDM")
for (i in seq_along(ldata)){
  bmr[[i]] = ML.exec(dataset = ldata[[i]])# MEter aqui el parametro ID
}

#Rename Benchmarks by data input
namesbmr <- paste("Bmr_",l, sep="")
names(bmrs) = namesbmr
saveRDS(object = bmrs , file = paste0("projects/Entropy/data/benchmarks/",pattern,"_Benchmarks.rds"))
kk = readRDS("projects/Entropy/data/benchmarks/Sriniv_Amsel_Genus_C_Benchmarks.rds")
bmrs$Bmr_Sriniv_Amsel_Genus_AR_train_CLUST_5.rds$results$dataset$classif.ksvm.tuned$measures.test
