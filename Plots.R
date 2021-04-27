require(ggplot2)
require(viridis)
library(ggvenn)
require(plotly)
## Fig1
# A) PLOTS COUNTS VS RELA ABUN
comp_C = ggplot(Counts, aes(x = FS, y = AUC, color = Algorithm, shape = Algorithm)) +
  geom_point() +
  geom_line(aes(group = Algorithm, color = Algorithm))+
  scale_color_manual(values = viridis(3))+
  theme_light()+
  theme( axis.text=element_text(size=20),axis.title.y = element_text(size = 30),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=25), 
         legend.text=element_text(size=19))
comp_C

comp_AR = ggplot(RA, aes(x = FS, y = AUC, color = Algorithm, shape = Algorithm)) +
  geom_point() +
  geom_line(aes(group = Algorithm, color = Algorithm))+
  scale_color_manual(values = viridis(3))+
  theme_light()+
  theme( axis.text=element_text(size=20),axis.title.y = element_text(size = 30),axis.title.x = element_blank(), axis.ticks.x = element_blank(),legend.title=element_text(size=25), 
         legend.text=element_text(size=19))
comp_AR

# B) Plot comparacon de features
CLUST_C = readRDS("projects/Entropy/data/CLUST_C_names.rds")
KW_C = readRDS("projects/Entropy/data/train/Ravel_Genus_C_train_KW_12.rds")
FCBF_C = readRDS("projects/Entropy/data/train/Ravel_Genus_C_train_FCBF_8.rds")
LDM_C = readRDS("projects/Entropy/data/train/Ravel_Genus_C_train_LDM_40.rds")
list_C = list(CLUST_C, names(KW_C), names(FCBF_C), names(LDM_C))
names = c("CLUST","KW", "FCBF","LDM")
remove = "target"
names(list_C) = names
for (i in seq_along(list_C)){
  list_C[[i]] = list_C[[i]][!list_C[[i]] %in% remove]
}
library("ggVennDiagram")
venn_C  = ggVennDiagram(list_C,color = viridis(1),lty = 1, label_alpha = 0, label= "count" ) 
venn_C  = venn_C  + scale_fill_gradient(high="#21908CFF",low = "#FDE725FF")
venn_C  = venn_C  + labs(fill = "Features")
venn_C

CLUST_AR = readRDS("projects/Entropy/data/CLUST_AR_names.rds")
KW_AR = readRDS("projects/Entropy/data/train/Ravel_Genus_AR_train_KW_12.rds")
FCBF_AR = readRDS("projects/Entropy/data/train/Ravel_Genus_AR_train_FCBF_9.rds")
LDM_AR = readRDS("projects/Entropy/data/train/Ravel_Genus_AR_train_LDM_39.rds")
list_AR = list(CLUST_AR, names(KW_AR), names(FCBF_AR), names(LDM_AR))
names = c("CLUST","KW", "FCBF","LDM")
remove = "target"
names(list_AR) = names
for (i in seq_along(list_AR)){
  list_AR[[i]] = list_AR[[i]][!list_AR[[i]] %in% remove]
}
venn_AR  = ggVennDiagram(list_AR,color = viridis(1),lty = 1, label_alpha = 0, label= "count" ) 
venn_AR  = venn_AR  + scale_fill_gradient(high="#21908CFF",low = "#FDE725FF")
venn_AR  = venn_AR  + labs(fill = "Features")
venn_AR


# C) Comparacion de los mejores modelos de cada tipo de datos
# D) Variable importance del mejor modelo en general

## Fig2
# AUROC Curve Predict

## Fig3
# Heatmap Phyloseq













 



# Curva ROC
a = asROCRPrediction(prediccion)
p = ROCR::performance(a, "tpr", "fpr")
plot(p)




