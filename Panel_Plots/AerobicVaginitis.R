# Boxplots bichejos
library(phyloseq)
library(data.table)
library(ggpubr)
library(viridis)
library(mlr)
norm.dataset = function(data){
  # Retain only numerics variables
  cols <- sapply(data, is.numeric) 
  
  # Normalize that variables
  data[cols] <- apply(X = data[cols], FUN = function(x) log2(x+1), MARGIN = 2) 
  return(data) 
}

# Ravel
Ravel.phy = readRDS("git/BVMetaGenomics/data/Phyloseqs/Ravel_phyloseq.rds")
Ravel.phy = tax_glom(Ravel.phy, taxrank = "Rank7")

# Preparar los nombres de los bichejos
otus = as.data.frame(Ravel.phy@otu_table)
corres = as.data.frame(Ravel.phy@tax_table@.Data[,7])
rownames(otus) = corres$`Ravel.phy@tax_table@.Data[, 7]`
otus = as.data.frame(t(otus))

# Sacar las etiquetas del clusterizado
etiq = readRDS("git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/Ravel_Clust3.rds")
identical(rownames(otus), rownames(etiq))
cluster = etiq$cluster
otus = cbind(otus,cluster)
colnames(otus)[1:557] = substr(colnames(otus)[1:557], start = 4, stop = 100) 
sort(colnames(otus))
feat = c("Enterococcus faecalis","Staphylococcus aureus","Streptococcus agalactiae", "cluster")
otus = otus[feat]
otus = norm.dataset(otus)
otus = normalizeFeatures(
  otus,
  method = "range",
  cols = NULL,
  range = c(0, 1),
  on.constant = "quiet")

xx = otus
n = names(xx)
l = list()
for (i in 1:(length(n) - 1)) {
  l[[i]] = cbind.data.frame(species = rep(n[i], nrow(xx)), values = xx[,i], Cluster = xx$cluster)
}
# As?, para cada item de la list tengo un dataframe de dim = (nrow(xx), 3), y puedo por lo tanto hacer un rbind de los elementos de la lista
data = rbindlist(l)
data$species = as.vector(as.character(data$species))
# Ploteo directamente
p = ggplot(data = data, aes(x = species, y = values)) +
  geom_boxplot(aes(fill=Cluster))+ theme_light()+ theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p = p + facet_wrap( ~ species, scales = 'free') + stat_compare_means(aes(group = Cluster), label = 'p.format',vjust = 0.7)
p = change_palette(p, palette= c("#21908CFF", "#440154FF","#FDE725FF"))
# Extract the legend. Returns a gtable
leg <- get_legend(p)
# Convert to a ggplot and print
as_ggplot(leg)
p = p + theme(legend.position = "none") 
p

# Srinivasan
Srin.phy = readRDS("git/BVMetaGenomics/data/Phyloseqs/Sriniv_phyloseq.rds")
Srin.phy = tax_glom(Srin.phy, taxrank = "Rank7")
Srin.phy =  subset_taxa(Srin.phy, !Rank7 %in% c("s__NA"))
ts = as.data.frame(Srin.phy@otu_table)
co = as.data.frame(Srin.phy@tax_table@.Data[,7])
rownames(ts)= co$`Srin.phy@tax_table@.Data[, 7]`
ts = as.data.frame(t(ts))

# Sacar las etiquetas del clusterizado
labels = readRDS("git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/Sriniv_Clust3.rds")
identical(rownames(ts), rownames(labels))
cluster = labels$cluster
ts = cbind(ts,cluster)
colnames(ts)[1:132] = substr(colnames(ts)[1:132], start = 4, stop = 100) 
sort(colnames(ts))
cols = c("Escherichia coli","Streptococcus agalactiae", "cluster")
ts = ts[cols]
ts = norm.dataset(ts)
ts = normalizeFeatures(
  ts,
  method = "range",
  cols = NULL,
  range = c(0, 1),
  on.constant = "quiet")
ns = names(ts)
l2 = list()
for (i in 1:(length(ns) - 1)) {
  l2[[i]] = cbind.data.frame(species = rep(ns[i], nrow(ts)), values = ts[,i], Cluster = ts$cluster)
}
# As?, para cada item de la list tengo un dataframe de dim = (nrow(xx), 3), y puedo por lo tanto hacer un rbind de los elementos de la lista
data2 = rbindlist(l2)
data2$species = as.vector(as.character(data2$species))
# Ploteo directamente
p = ggplot(data = data2, aes(x = species, y = values)) +
  geom_boxplot(aes(fill=Cluster))+ theme_light()+ theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p = p + facet_wrap( ~ species, scales = 'free') + stat_compare_means(aes(group = Cluster), label = 'p.format',vjust = 0.7)
p = change_palette(p, palette= c("#21908CFF", "#440154FF","#FDE725FF"))
p = p + theme(legend.position = "none") 
p