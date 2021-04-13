source("git/Entropy/ML_Functions.R")
source("git/Entropy/Preprocess_Functions.R")

Genus_train_KW = kruskal.FS(data = train, fs.type ='kruskal.test', nfeat = 15 )
Genus_train_FCBF = FCBF.FS(data = train, thold = 0.005)
Genus_train_LDM = LDM.FS(data = train, seed = 17, method = "bray", thold = 0.1)
Genus

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

BMR_FCBF = ML.exec(Genus_train_FCBF)
BMR_KW = ML.exec(Genus_train_KW)
BMR_LDM = ML.exec(Genus_train_LDM)
BMR_Train = ML.exec(train)
