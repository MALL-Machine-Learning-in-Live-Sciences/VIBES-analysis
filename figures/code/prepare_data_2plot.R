# Figure 1
packgs <- c("phyloseq")
lapply(packgs, require, character.only = TRUE)
setwd("~/git/BV_Microbiome/")
rank = "Species"
nfeat = "22"
cl <- list.files(path = "02_cluster/data/C4/", pattern = paste(rank,nfeat,sep = "_"))
vl <- list.files(path = "01_get_valencias/res/")
pseqs_list <- list()

for (i in seq_along(cl)) {
  # load pseq and valencias
  pseq <- readRDS(file = paste0("02_cluster/data/C4/", cl[i]))
  val <- read.csv2(file = paste0("01_get_valencias/res/", vl[i]),
                   header = TRUE, sep =",")
  # check same order
  val <- data.frame(val[,-1], row.names = val[,1])
  identical(rownames(pseq@sam_data), rownames(val))
  # add valencia 
  pseq@sam_data$CST <- val$CST
  # order pseq by target
  df <- pseq@sam_data
  df$cluster[df$cluster == "D"] <- "VCS-IV"
  df$cluster[df$cluster== "IDD"] <- "VCS-III"
  df$cluster[df$cluster == "IDN"] <- "VCS-II"
  df$cluster[df$cluster == "N"] <- "VCS-I"
  sample_data(pseq) <- df
  ordered_df <- df[order(df$cluster,decreasing = FALSE), ]
  pseq_sorted <- ps_reorder(ps = pseq, sample_order = rownames(ordered_df))
  pseqs_list[[i]] <- pseq_sorted
}
names(pseqs_list) <- sapply(strsplit(cl, "_"), "[", 1)
# get NS categories for prjna2085
pseqs_list$PRJNA2085@sam_data$Nugent_Score = replace(pseqs_list$PRJNA2085@sam_data$Nugent_Score,
                                                     pseqs_list$PRJNA2085@sam_data$Nugent_Score > 6, "high")
pseqs_list$PRJNA2085@sam_data$Nugent_Score = replace(pseqs_list$PRJNA2085@sam_data$Nugent_Score,
                                                     pseqs_list$PRJNA2085@sam_data$Nugent_Score > 3 & pseqs_list$PRJNA2085@sam_data$Nugent_Score < 7, "intermediate")
pseqs_list$PRJNA2085@sam_data$Nugent_Score = replace(pseqs_list$PRJNA2085@sam_data$Nugent_Score, pseqs_list$PRJNA2085@sam_data$Nugent_Score < 4, "low")

# Rename with new package names clusters
saveRDS(pseqs_list,file = "figures/data/f1_data.rds")

# Figure 2
# Forestplot
require(dplyr)
ext_val = readRDS("03_machine_learning/res/performances_extval.rds")
ext_val$cohort[ext_val$cohort == 'Ravel'] <- 'Train'
ext_val$cohort[ext_val$cohort == 'Sriniv'] <- 'Validation 1'
ext_val$cohort[ext_val$cohort == 'PRJNA2085'] <- 'Validation 2'
ext_val$cohort[ext_val$cohort == 'PRJNA7977'] <- 'Validation 3'
ext_val$cohort[ext_val$cohort == 'PRJNA3020'] <- 'Validation 4'
ext_val <- ext_val[c(4,5,1,3,2),]
ext_val$n <-c(394, 220, 1657, 220, 201)
names(ext_val)[2] <- "mean"
names(ext_val)[3] <- "lower"
names(ext_val)[4] <- "upper"
ext_val <- ext_val %>% 
  mutate_if(is.numeric, round, digits = 3)
base_data <- tibble::tibble(ext_val)
base_data <- base_data[-1,]
saveRDS(base_data,"figures/data/f2_forest_data.rds")

# Figure 2 Supp

toplot1 <- readRDS("03_machine_learning/res/summary_performances.rds")
toplot2 <- readRDS("03_machine_learning/res/summary_performances_extval.rds")
# Rename algorithms 
toplot2[toplot2 == "classif.glmnet.tuned"] <- "glmnet"
toplot2[toplot2 == "classif.ranger.tuned"] <- "rf"
toplot2[toplot2 == "classif.svm.tuned"] <- "svm"
toplot2[toplot2 == "classif.xgboost.tuned"] <- "xgboost"
# Rename variable
toplot2$variable <- as.character(toplot2$variable)
toplot2[toplot2 == "Brier Error"] <- "Brier"
toplot2[toplot2 == "Balanced Accuracy"] <- "BAccuracy"
toplot2$variable <- as.factor(toplot2$variable)
toplot <- rbind(toplot1, toplot2)
# Rename Studies
toplot[toplot == "Ravel"] <- "Discovery"
toplot[toplot == "Sriniv"] <- "Validation 1"
toplot[toplot == "PRJNA2085"] <- "Validation 2"
toplot[toplot == "PRJNA7977"] <- "Validation 3"
toplot[toplot == "PRJNA3020"] <- "Validation 4"
toplot$data <- factor(toplot$data,
                      levels = c("Validation 4","Validation 3","Validation 2",
                                 "Validation 1","Discovery"))
toplot$value <- round(toplot$value, digits = 2)

saveRDS(toplot,"figures/data/f2_supp_data.rds")

# Figure 3
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

d <- readRDS("03_machine_learning/data/PRJNA3020_clr.rds")
d <- subset(d,
            select = -c(Mageeibacillus_indolicus, Peptoniphilus_lacrimalis))
p <- pruned_predict(pruned_glmnet = model,
                    newdata = as.matrix(d[,1:20]))
identical(rownames(p), rownames(d))
dp <- cbind(d, p)
packgs <- c("phyloseq", "ComplexHeatmap", "microViz")
lapply(packgs, require, character.only = TRUE)
cl <- readRDS("02_cluster/data/C4/PRJNA3020_Cluster_Species_22_pseq.rds")
identical(rownames(dp), rownames(cl@sam_data))
sample_data(cl) <- cbind(cl@sam_data, dp[,22:26])

vl <- read.csv2(file  = "01_get_valencias/res/PRJNA3020.csv", sep = ",", header = TRUE)
# check same order
vl <- data.frame(vl[,-1], row.names = vl[,1])
identical(rownames(cl@sam_data), rownames(vl))
# add valencia 
cl@sam_data$CST <- vl$CST
# order pseq by probbilities of belong to own class
df <- cl@sam_data

dfn <- df[df$p_cluster == "N",]
dfn <- dfn[order(dfn$N, decreasing = TRUE),]
dfidn <- df[df$p_cluster == "IDN",]
dfidn <- dfidn[order(dfidn$IDN, decreasing = TRUE),]
dfidd <- df[df$p_cluster == "IDD",]
dfidd <- dfidd[order(dfidd$IDD, decreasing = TRUE),]
dfd <- df[df$p_cluster == "D",]
dfd <- dfd[order(dfd$D, decreasing = TRUE),]
pseq_sorted <- ps_reorder(ps = cl,
                          sample_order = rownames(rbind(dfn,
                                                        dfidn,
                                                        dfidd,
                                                        dfd
                                                        )))
df <- pseq_sorted@sam_data
df$cluster[df$cluster == "D"] <- "VCS-IV"
df$cluster[df$cluster== "IDD"] <- "VCS-III"
df$cluster[df$cluster == "IDN"] <- "VCS-II"
df$cluster[df$cluster == "N"] <- "VCS-I"

df$p_cluster[df$p_cluster == "D"] <- "VCS-IV"
df$p_cluster[df$p_cluster== "IDD"] <- "VCS-III"
df$p_cluster[df$p_cluster == "IDN"] <- "VCS-II"
df$p_cluster[df$p_cluster == "N"] <- "VCS-I"
sample_data(pseq_sorted) <- df
saveRDS(object = pseq_sorted, file = "figures/data/f3_data.rds")
