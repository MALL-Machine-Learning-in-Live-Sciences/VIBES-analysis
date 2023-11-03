# Jaccard
setwd("02_cluster/res/")
# Load BootCCP results
r_0.5 <- readRDS(file = "Bootstrap_0.5_CCP_Ravel_Species_22_clr_km/Bootstrap_0.5_CCP_Ravel_Species_22_clr_km.rds")
r_0.6 <- readRDS(file = "Bootstrap_0.625_CCP_Ravel_Species_22_clr_km/Bootstrap_0.625_CCP_Ravel_Species_22_clr_km.rds")
r_0.7 <- readRDS(file = "Bootstrap_0.75_CCP_Ravel_Species_22_clr_km/Bootstrap_0.75_CCP_Ravel_Species_22_clr_km.rds")
r_0.9 <- readRDS(file = "Bootstrap_0.9_CCP_Ravel_Species_22_clr_km/Bootstrap_0.9_CCP_Ravel_Species_22_clr_km.rds")
r_1 <- readRDS(file = "Bootstrap_1_CCP_Ravel_Species_22_clr_km/Bootstrap_1_CCP_Ravel_Species_22_clr_km.rds")

# Extracting cluster labels vectors
c5 <- r_0.5$ccp[[4]]$consensusClass
c6 <- r_0.6$ccp[[4]]$consensusClass
c7 <- r_0.7$ccp[[4]]$consensusClass
c9 <- r_0.9$ccp[[4]]$consensusClass
c1 <- r_1$ccp[[4]]$consensusClass


jaccard_similarity <- function(A, B) { 
  intersection = length(intersect(A, B)) 
  union = length(A) + length(B) - intersection 
  return (intersection/union) 
}

# c5 vs c1 
shared <- intersect(names(c5), names(c1))
c1_ordered <- c1[names(c1) %in% shared]
identical(sort(names(c1_ordered)), sort(names(c5)))
c1_ordered <- c1_ordered[order(names(c1_ordered))]
c5 <- c5[order(names(c5))]

# Contar coincidencias
coincidencias <- sum(c5 == c1_ordered)
# Imprimir el número de coincidencias
cat("Número de coincidencias:", coincidencias, "\n")
table(c5, c1_ordered) # Están cambiados el 1 por el 3
ri_1 <- adjustedRandIndex(x = c5, y = c1_ordered)
ri_1
jaccard_similarity(c5,c1_ordered)




# c6 vs c1 
shared <- intersect(names(c6), names(c1))
c1_ordered <- c1[names(c1) %in% shared]
identical(sort(names(c1_ordered)), sort(names(c6)))
c1_ordered <- c1_ordered[order(names(c1_ordered))]
c6 <- c6[order(names(c6))]
# Contar coincidencias
coincidencias <- sum(c6 == c1_ordered)
# Imprimir el número de coincidencias
cat("Número de coincidencias:", coincidencias, "\n")
table(c6, c1_ordered) # Están cambiados el 1 por el 3
ri_2 <- adjustedRandIndex(x = c6, y = c1_ordered)
ri_2
jaccard_similarity(c6,c1_ordered)


# c7 vs c1 
shared <- intersect(names(c7), names(c1))
c1_ordered <- c1[names(c1) %in% shared]
identical(sort(names(c1_ordered)), sort(names(c7)))
c1_ordered <- c1_ordered[order(names(c1_ordered))]
c7 <- c7[order(names(c7))]
# Contar coincidencias
coincidencias <- sum(c7 == c1_ordered)
# Imprimir el número de coincidencias
cat("Número de coincidencias:", coincidencias, "\n")
table(c7, c1_ordered) # Están cambiados el 1 por el 3, el 3 por el 4 y el 4 por el 1
ri_3 <- adjustedRandIndex(x = c7, y = c1_ordered)
ri_3
jaccard_similarity(c7,c1_ordered)

# c9 vs c1 
shared <- intersect(names(c9), names(c1))
c1_ordered <- c1[names(c1) %in% shared]
identical(sort(names(c1_ordered)), sort(names(c9)))
c1_ordered <- c1_ordered[order(names(c1_ordered))]
c9 <- c9[order(names(c9))]
# Contar coincidencias
coincidencias <- sum(c9 == c1_ordered)
# Imprimir el número de coincidencias
cat("Número de coincidencias:", coincidencias, "\n")
table(c9, c1_ordered) # Están cambiados el 1 por el 3, el 3 por el 4
ri_4 <- adjustedRandIndex(x = c9, y = c1_ordered)
ri_4
jaccard_similarity(c9,c1_ordered) 

ri <- c(ri_1, ri_2, ri_3, ri_4)
print(paste( "RI mean across all subsampling:", mean(ri)))
print(paste( "RI median across all subsampling:", median(ri)))
