#ML CESGA
setwd('/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data')
path = "/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/train/"
source("/home/ulc/co/dfe/git/Entropy/functions/FunctionsML.R")
#path = "projects/Entropy/data/train/"
pattern = "Ravel_Genus_C"
l = list.files(path, pattern = pattern)

ldata = list()
for (i in 1:length(l)) {
  ldata[[i]] = readRDS(paste(path, l[[i]], sep=''))
}
names(ldata) = l
bmrs = list()
for (i in seq_along(ldata)){
  bmrs[[i]] = ML.exec(dataset = ldata[[i]])
}
#Rename Benchmarks by data input
namesbmr <- paste("Bmr_",l, sep="")
names(bmrs) = namesbmr

saveRDS(object = bmrs , file = paste0("/mnt/netapp2/Store_uni/home/ulc/co/dfe/projects/Entropy/data/benchmarks/",pattern,"Benchmarks.rds"))
