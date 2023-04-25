# Predict with glmnet pruned
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
  response$class <- class[,]
  return(response)
}

a = readRDS("00_preprocess_cohorts/data/SpeciesIntersect/PRJNA2085_Species_pseq_22.rds")
# CLR con 20 species en phyloseq
mntin <- taxa_names(a)
rmv <- c("Mageeibacillus_indolicus", "Peptoniphilus_lacrimalis")
mntin <- mntin [! mntin %in% rmv]

a2 = prune_taxa(taxa = mntin,x = a)

a2 = microbiome::transform(a2, transform =  'clr')
get_dataset <- function(phyobject){
  require(phyloseq)
  clinics = as.data.frame(sample_data(phyobject))
  otus = as.data.frame(get_taxa(phyobject))
  identical(rownames(clinics), rownames(otus))
  data = cbind(otus, clinics)
  return(data)
}

a2 <- get_dataset(a2)
a2 <- a2 [,1:20]
a2 = as.matrix(a2)
pp2 <- pruned_predict(pruned_glmnet =pr_gl,
                      newdata = a2 )



# CLR con 20 species en marix

# CLR con 22
ab = readRDS("03_machine_learning/data/PRJNA2085_clr.rds")



object = readRDS("glmnet.rds")
test = readRDS("~/git/BV_Microbiome/03_machine_learning/data/PRJNA2085_clr.rds")
test = test[,1:22]
newx = as.matrix(test)
type = "class"
pp = glmnet:::predict.multnet(object = object,
                              newx = newx,
                              s = object$lambda[98],
                              type = type)
pp
table(pp)

test2 <- subset(test,
                select = -c(Mageeibacillus_indolicus,
                            Peptoniphilus_lacrimalis))

pr_gl <- readRDS(file = "03_machine_learning/model/pruned_glmnet.rds")
newx2 = as.matrix(test2)
pp2 <- pruned_predict(pruned_glmnet =pr_gl,
                      newdata = newx2 )
pp2
table(pp2$class)
