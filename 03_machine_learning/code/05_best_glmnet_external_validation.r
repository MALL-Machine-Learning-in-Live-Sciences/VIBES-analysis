#### External Validation ####
# 1.Load model and declare function
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

# 2.Delcare paths to data 
path_test <- "03_machine_learning/data/"
lf <- list.files(path = path_test, pattern = "clr")
lf <- lf[-grep("Ravel", lf)]
pl <- list()
for (i in seq_along(lf)) {
  # 3.Load data and remove species
  d <- readRDS(paste0("03_machine_learning/data/", lf[i]))
  d <- subset(d,
              select = -c(Mageeibacillus_indolicus, Peptoniphilus_lacrimalis))
  # 4.Predict new classes
  p <- pruned_predict(pruned_glmnet = model,
                      newdata = as.matrix(d[,1:20]))
  identical(rownames(p), rownames(d))
  dp <- cbind(d, p)
  pl[[i]] <- dp
}
names(pl) <- sapply(strsplit(lf, "_"), "[", 1)
# 4.Extract performance 
l_msr <- list()
for (i in seq_along(pl)) {
  require(caret)
  cm <- confusionMatrix(data = as.factor(pl[[i]]$p_cluster),
                        reference = as.factor(pl[[i]]$cluster),
                        positive = "N")
  acc <- cm$overall[[1]]
  kap <- cm$overall[[2]]
  acc_low <- cm$overall[[3]]
  acc_upp <- cm$overall[[4]]
  
  bacc <- mlr3measures::bacc(truth = as.factor(pl[[i]]$cluster),
                     response = as.factor(pl[[i]]$p_cluster))
  brier <- mlr3measures::mbrier(truth = as.factor(pl[[i]]$cluster),
                       prob = as.matrix(pl[[i]][,22:25]))
  df <- data.frame (
    cohort = names(pl)[i],
    acc  = acc,
    acc_low = acc_low,
    acc_upp = acc_upp,
    kappa = kap,
    bacc = bacc,
    brier = brier)
  l_msr[[i]] <- df
}
ext_val <- data.table::rbindlist(l_msr)
# Save summary
saveRDS(ext_val, file = "03_machine_learning/res/performances_extval.rds")
