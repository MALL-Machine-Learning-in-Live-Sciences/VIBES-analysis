# Figure 4
library(pheatmap)
library(dplyr)

# HeatMap Ravel 8 feat Genera
a = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Cluster/3C/FS/Ravel_3C_8.rds")
ordenG = c("Prevotella", "Fannyhessea","Dialister","Sneathia","Anaerococcus","Lactobacillus","Peptoniphilus", "Gardnerella")
feat = c(ordenG,"cluster","target")
b = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Cluster/3C/Ravel_Clust3.rds")
data = b[feat]
names(data)[2] <- 'Atopobium'

data$cluster = as.factor(data$cluster)
levels(data$cluster) <- c("C3", "C2", "C1")
data$cluster = as.character(data$cluster)
data = data[order(data[,9]),]


mat = data[1:8]
mat.scale <- scale(mat, center = T, scale = T) 
mat.scale <- scale(t(mat.scale), center = T, scale = T)
names(data)[9] <- "Clusters" 
names(data)[10] <- "Nugent"
anno_col <- data.frame(Clusters = data[9], Nugent = data[10])
anno_metabo_colors <- list(Nugent = c("high" = "#440154FF", "intermediate" = "#21908CFF", "low" = "#FDE725FF"),
                            Clusters =c("C1" ="#21908CFF" , "C2"= "#440154FF", "C3" = "#FDE725FF"))
pheatmap(mat.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = F, 
         fontsize_row = 11, fontsize_col = 8,
         fontsize = 10,
         show_colnames = FALSE,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = anno_col,
         annotation_colors = anno_metabo_colors,
         border_color = 'NA')


# HeatMap Srinivasan 8 feat Genera
a = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Cluster/3C/FS/Ravel_3C_8.rds")
ordenG = c("Prevotella", "Fannyhessea","Dialister","Sneathia","Anaerococcus","Lactobacillus","Peptoniphilus", "Gardnerella")
feat = c(ordenG,"cluster", "target", "amsel")
b = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Cluster/3C/Sriniv_Clust3.rds")
c = readRDS("~/git/BVMetaGenomics/data/GenusIntersect/Sriniv_Genus_Counts.rds")
identical(rownames(b), rownames(c))
amsel = c$amsel
data = cbind(b, amsel)
data = data[feat]
data = data[order(data[,9]),]
remove = c("X1130")
data = data[ !(rownames(data) %in% remove), ]
names(data)[2] <- 'Atopobium'
data$cluster = as.factor(data$cluster)
levels(data$cluster) <- c("C3", "C2", "C1")
data$cluster = as.character(data$cluster)
data = data[order(data[,9]),]

mat = data[1:8]
# scale on OTUs
mat.scale <- scale(mat, center = T, scale = T) 
# scale on samples
mat.scale <- scale(t(mat.scale), center = T, scale = T) 
names(data)[9] <- "Clusters" 
names(data)[10] <- "Nugent"
names(data)[11] <- "Amsel"
anno_col <- data.frame(Clusters = data[9], Nugent = data[10], Amsel = data[11])
anno_metabo_colors <- list( Amsel = c("pos" = "#440154FF", "neg" = "#FDE725FF"),
                            Nugent = c("high" = "#440154FF", "intermediate" = "#21908CFF", "low" = "#FDE725FF"),
                            Clusters =c("C1" ="#21908CFF" , "C2"= "#440154FF", "C3" = "#FDE725FF"))


pheatmap(mat.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = F, 
         fontsize_row = 11, fontsize_col = 8,
         fontsize = 10,
         show_colnames = FALSE,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = anno_col,
         annotation_colors = anno_metabo_colors,
         border_color = 'NA')



# HeatMap Ravel 8 feat Species
a = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/FS/Ravel_3C_8.rds")
ordenS = c("Lactobacillus.crispatus", "Fannyhessea.vaginae", "Prevotella.buccalis",
           "Sneathia.sanguinegens", "Dialister.micraerophilus", "Prevotella.timonensis",
           "Prevotella.amnii", "Gardnerella.vaginalis")
feat = c(ordenS,"cluster", "target")
b = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/Ravel_Clust3.rds")
data = b[feat]
data = data[order(data[,9]),]
names(data)[2] <- 'Atopobium.vaginae'

mat = data[1:8]
mat.scale <- scale(mat, center = T, scale = T) 
mat.scale <- scale(t(mat.scale), center = T, scale = T)
names(data)[9] <- "Clusters" 
names(data)[10] <- "Nugent"
anno_col <- data.frame(Clusters = data[9], Nugent = data[10])
anno_metabo_colors <- list(Nugent = c("high" = "#440154FF", "intermediate" = "#21908CFF", "low" = "#FDE725FF"),
                           Clusters =c("C1" ="#21908CFF" , "C2"= "#440154FF", "C3" ="#FDE725FF" ))
pheatmap(mat.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = F, 
         fontsize_row = 11, fontsize_col = 8,
         fontsize = 10,
         show_colnames = FALSE,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = anno_col,
         annotation_colors = anno_metabo_colors,
         border_color = 'NA')



#HeatMap Srinivasan 8 feat Species
a = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/FS/Ravel_3C_8.rds")
ordenS = c("Lactobacillus.crispatus", "Fannyhessea.vaginae", "Prevotella.buccalis",
           "Sneathia.sanguinegens", "Dialister.micraerophilus", "Prevotella.timonensis",
           "Prevotella.amnii", "Gardnerella.vaginalis")
feat = c(ordenS,"cluster", "target","amsel")
b = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/Sriniv_Clust3.rds")
c = readRDS("~/git/BVMetaGenomics/data/SpeciesIntersect/Sriniv_Sps_Counts.rds")
identical(rownames(b), rownames(c))
amsel = c$amsel
data = cbind(b, amsel)
data = data[feat]
data = data[order(data[,9]),]
remove = c("X1130")
data = data[ !(rownames(data) %in% remove), ]
names(data)[2] <- 'Atopobium.vaginae'
mat = data[1:8]
# scale on OTUs
mat.scale <- scale(mat, center = T, scale = T) 
# scale on samples
mat.scale <- scale(t(mat.scale), center = T, scale = T) 
names(data)[9] <- "Clusters" 
names(data)[10] <- "Nugent"
names(data)[11] <- "Amsel"
anno_col <- data.frame(Clusters = data[9], Nugent = data[10], Amsel = data[11])
anno_metabo_colors <- list( Amsel = c("pos" = "#440154FF", "neg" = "#FDE725FF"),
                            Nugent = c("high" = "#440154FF", "intermediate" = "#21908CFF", "low" = "#FDE725FF"),
                            Clusters =c("C1" ="#21908CFF", "C2"= "#440154FF", "C3" = "#FDE725FF"))


pheatmap(mat.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = F, 
         fontsize_row = 11, fontsize_col = 8,
         fontsize = 10,
         show_colnames = FALSE,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = anno_col,
         annotation_colors = anno_metabo_colors,
         border_color = 'NA')
