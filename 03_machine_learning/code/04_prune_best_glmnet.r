#### Prune Glmnet Model ####
require(glmnet)
# 1.Load glmnet model
#FALTA METER UN SCRIPT DONDE SE VEA COMO JOSE REENTRNA Y EXTRAE EL MODELO (ASÍ COMO LA PERFORMANCE)
object = readRDS("03_machine_learning/res/glmnet.rds")

# 2.Prepare data before processing for each class
s = object$lambda[98]
a0 = object$a0
rownames(a0) = rep("(Intercept)", nrow(a0))
nbeta = object$beta
klam = dim(a0)
nclass = klam[[1]]
nlambda = length(s)
lambda = object$lambda
lamlist = glmnet:::lambda.interp(lambda, s)

# 3.Intercept and betas of each species are retained for each class
for (i in seq(nclass)) {
  # Species whose betas are not 0 in s98 are selected from each class
  b0 <- nbeta[[i]][,98]
  b0 <- b0[b0 != 0]
  nbeta[[i]] <- nbeta[[i]][names(b0),]
  # Intercepts are binded to species
  kbeta = methods::rbind2(a0[i, , drop = FALSE], nbeta[[i]])
  vnames = dimnames(kbeta)[[1]]
  dimnames(kbeta) = list(NULL, NULL)
  kbeta = kbeta[, lamlist$left, drop = FALSE] %*% Diagonal(x=lamlist$frac) +
    kbeta[, lamlist$right, drop = FALSE] %*% Diagonal(x=1 - lamlist$frac)
  dimnames(kbeta) = list(vnames, paste(seq(along = s)))
  # Species with the selected lambda intercept is reassigned to the nbeta object
  nbeta[[i]] = kbeta
}

# 4.The nbetas, nclass and nlambdas are stored in a new object.
glmnet_pruned <- list(nbeta, nclass, nlambda)
names(glmnet_pruned) <- c("nbeta", "nclass", "nlambda")
saveRDS(object = glmnet_pruned, "03_machine_learning/model/pruned_glmnet.rds")