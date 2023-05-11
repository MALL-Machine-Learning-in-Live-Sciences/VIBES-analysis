##### Experiment D0 clusters vs CSTs #####
PRJNA3020 = readRDS("02_cluster/data/C4/PRJNA3020_Cluster_Species_22_pseq.rds")

meta <- data.frame(sample_data(PRJNA3020))
d0 <- meta[grep("D0", meta$sample_alias), ]
PRJNA3020 <- subset_samples(PRJNA3020,
                            (sample_names(PRJNA3020) %in% rownames(d0)))
# Solo clusters
d = as.data.frame(get_taxa(PRJNA3020))
d <- subset(d,
            select = -c(Mageeibacillus_indolicus, Peptoniphilus_lacrimalis))
model <- readRDS("03_machine_learning/model/pruned_glmnet.rds")
pruned_predict <- function(pruned_glmnet, newdata){
  require(glmnet)
  # 1.Pre-processing of data for betas multiplication
  dd = dim(newdata)
  if (inherits(newdata, "sparseMatrix"))
    newdata = as(newdata, "dMatrix")
  npred = dd[[1]]
  dn = list(names(pruned_glmnet$nbeta), dimnames(pruned_glmnet$nbeta[[1]])[[2]], dimnames(newdata)[[1]])
  dp = array(0, c(pruned_glmnet$nclass, pruned_glmnet$nlambda, npred), dimnames = dn)
  # 2.Multiplication by betas
  for (i in seq(pruned_glmnet$nclass)) {
    rn <- rownames(pruned_glmnet$nbeta[[i]])
    rn <- rn[-1]
    fitk = cbind2(1, newdata[,rn]) %*% (pruned_glmnet$nbeta[[i]])
    dp[i, , ] = dp[i, , ] + t(as.matrix(fitk))
  }
  # 3.Results are extracted
  pp = exp(dp)
  psum = apply(pp, c(2, 3), sum)
  # Response
  response = aperm(pp/rep(psum, rep(pruned_glmnet$nclass, pruned_glmnet$nlambda * npred)), c(3,1, 2))
  response <- as.data.frame(response)
  colnames(response) <- substr(colnames(response), 1, nchar(colnames(response))-2)
  cp = aperm(dp, c(3, 1, 2))
  # Class
  class <- apply(cp, 3, glmnet:::glmnet_softmax)
  #Bind both
  response$p_cluster <- class[,]
  return(response)
}
r <- pruned_predict(pruned_glmnet = model, newdata = as.matrix(d))

identical(rownames(r), rownames(PRJNA3020@sam_data))
doc <- cbind(r, PRJNA3020@sam_data$status)
names(doc)[6] <- "status"
saveRDS(object = doc, file = "04_treatment/data/D0_only_clusters.rds")

# Solo valencias
valencias <- data.frame(read_csv("01_get_valencias/res/PRJNA3020.csv"))
valencias <- data.frame(valencias[,-1], row.names = valencias[,1])
valencias <- valencias[rownames(d0),]
valencias <- valencias[,99:ncol(valencias)]
valencias <- valencias[,!names(valencias) %in% c("subCST")]
identical(rownames(valencias), rownames(PRJNA3020@sam_data))
dov <- cbind(valencias, PRJNA3020@sam_data$status)
names(dov)[16] <- "status" 
saveRDS(object = dov, file = "04_treatment/data/D0_only_valencias.rds")

# Ambos
identical(rownames(r), rownames(valencias))
dall <- cbind(r, valencias, PRJNA3020@sam_data$status)
names(dall)[21] <- "status"
saveRDS(object = dall, file = "04_treatment/data/D0_all.rds")
