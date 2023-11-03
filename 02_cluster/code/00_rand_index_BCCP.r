# ARI
library(mclust)
# Load BootCCP results
r_0.5 <- readRDS(file = "02_cluster/res/Bootstrap_0.5_CCP_Ravel_Species_22_clr_km/Bootstrap_0.5_CCP_Ravel_Species_22_clr_km.rds")
r_0.6 <- readRDS(file = "02_cluster/res/Bootstrap_0.625_CCP_Ravel_Species_22_clr_km/Bootstrap_0.625_CCP_Ravel_Species_22_clr_km.rds")
r_0.7 <- readRDS(file = "02_cluster/res/Bootstrap_0.75_CCP_Ravel_Species_22_clr_km/Bootstrap_0.75_CCP_Ravel_Species_22_clr_km.rds")
r_0.9 <- readRDS(file = "02_cluster/res/Bootstrap_0.9_CCP_Ravel_Species_22_clr_km/Bootstrap_0.9_CCP_Ravel_Species_22_clr_km.rds")
r_1 <- readRDS(file = "02_cluster/res/Bootstrap_1_CCP_Ravel_Species_22_clr_km/Bootstrap_1_CCP_Ravel_Species_22_clr_km.rds")


# Mean cluster-consensus values

resultados_0.5 <- aggregate(clusterConsensus ~ k, data = r_0.5$icl$clusterConsensus,
                        FUN = function(x) c(Mean = mean(x), Median = median(x)))
resultados_0.6 <- aggregate(clusterConsensus ~ k, data = r_0.6$icl$clusterConsensus,
                         FUN = function(x) c(Mean = mean(x), Median = median(x)))
resultados_0.7 <- aggregate(clusterConsensus ~ k, data = r_0.7$icl$clusterConsensus,
                            FUN = function(x) c(Mean = mean(x), Median = median(x)))
resultados_0.9 <- aggregate(clusterConsensus ~ k, data = r_0.9$icl$clusterConsensus,
                            FUN = function(x) c(Mean = mean(x), Median = median(x)))


# Extracting cluster labels vectors
c5 <- r_0.5$ccp[[4]]$consensusClass
c6 <- r_0.6$ccp[[4]]$consensusClass
c7 <- r_0.7$ccp[[4]]$consensusClass
c9 <- r_0.9$ccp[[4]]$consensusClass
c1 <- r_1$ccp[[4]]$consensusClass

rbind()
# c5 vs c1 
shared <- intersect(names(c5), names(c1))
c1_ordered <- c1[names(c1) %in% shared]
identical(sort(names(c1_ordered)), sort(names(c5)))
c1_ordered <- c1_ordered[order(names(c1_ordered))]
c5 <- c5[order(names(c5))]
ri_1 <- adjustedRandIndex(x = c5, y = c1_ordered)
ri_1

# c6 vs c1 
shared <- intersect(names(c6), names(c1))
c1_ordered <- c1[names(c1) %in% shared]
identical(sort(names(c1_ordered)), sort(names(c6)))
c1_ordered <- c1_ordered[order(names(c1_ordered))]
c6 <- c6[order(names(c6))]
ri_2 <- adjustedRandIndex(x = c6, y = c1_ordered)
ri_2


# c7 vs c1 
shared <- intersect(names(c7), names(c1))
c1_ordered <- c1[names(c1) %in% shared]
identical(sort(names(c1_ordered)), sort(names(c7)))
c1_ordered <- c1_ordered[order(names(c1_ordered))]
c7 <- c7[order(names(c7))]
ri_3 <- adjustedRandIndex(x = c7, y = c1_ordered)
ri_3


# c9 vs c1 
shared <- intersect(names(c9), names(c1))
c1_ordered <- c1[names(c1) %in% shared]
identical(sort(names(c1_ordered)), sort(names(c9)))
c1_ordered <- c1_ordered[order(names(c1_ordered))]
c9 <- c9[order(names(c9))]
ri_4 <- adjustedRandIndex(x = c9, y = c1_ordered)
ri_4


ri <- c(ri_1, ri_2, ri_3, ri_4)
print(paste( "RI mean across all subsampling:", round(mean(ri),digits = 3)))
print(paste( "RI median across all subsampling:", round(median(ri),digits = 3)))
