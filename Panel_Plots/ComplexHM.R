# Packages Load 
library(phyloseq)
library(data.table)
library(BVML) #Predict package
library(metagMisc) #phyloseq_filters
# Complex Heatmap packages
require(ComplexHeatmap)
library(circlize)
require(viridis)

# Relative Abundances Function
relat.abun = function(phyobject){
  require(phyloseq)
  phy = transform_sample_counts(phyobject, function(x) x / sum(x))
  return(phy)
}

# Samples Ordering by prob.
#####
# Ravel
data = readRDS("git/BVMetaGenomics/data/SpeciesIntersect/Ravel_Sps_Counts.rds")
pre = BVML::BVClassify(dataset = data,type = "species",plot_HM = TRUE)  
R_order = pre$results

R_C1 = subset(R_order,response == "C1") #C1 Omly
R_C1 = R_C1[order(R_C1[,2], decreasing = FALSE),]

R_C2 = subset(R_order, response == "C2") #C2 Only
R_C2 = R_C2[order(R_C2[,2], decreasing = FALSE),]

R_C3 = subset(R_order, response == "C3") #C3 Only
R_C3 = R_C3[order(R_C3[,3], decreasing = TRUE),]

R_samp = rownames(rbind(R_C3, R_C1, R_C2)) # Ravel Samples Ordered by probability of belonging to each class

#Srinivasan
S_order = readRDS("git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/VAL1-27.rds")
S_order = S_order$Predictions$Ravel_3C_8$data

S_C1 = subset(S_order, response == "C1") #C1 Only
S_C1 = S_C1[order(S_C1[,4],decreasing = FALSE),]

S_C2 =  subset(S_order, response == "C2") # C2 Only
S_C2 = S_C2[order(S_C2[,4], decreasing = FALSE),]

S_C3 =  subset(S_order, response == "C3") # C3 only
S_C3 = S_C3[order(S_C3[,5], decreasing = TRUE),]

S_samp = rownames(rbind(S_C3, S_C1, S_C2 )) # Srinivasan Samples Ordered by probability of belonging to each class
#####

# Ravel
# Phylum level
#####
Ravel.phy = readRDS("git/BVMetaGenomics/data/Phyloseqs/Ravel_phyloseq.rds") #Reset phyloseq
Ravel.phy = tax_glom(Ravel.phy, taxrank = "Rank7")
Ravel.phy <- filter_taxa(Ravel.phy, function(x) sum(x > 20) > (0.08*length(x)), TRUE)
Ravel.phy = relat.abun(Ravel.phy)
r = as.data.frame(Ravel.phy@tax_table)
r2 = unique(r$Rank2)


R_Phyla = list()
for (i in seq_along(r2)) {
  sub = subset_taxa(Ravel.phy, Rank2 %in% c(r2[i]))
  data = as.data.frame(sub@tax_table[,7])
  nmes = data$Rank7
  R_Phyla [[i]] = nmes
}
names(R_Phyla) =r2

tax = as.data.frame(Ravel.phy@tax_table)
otu = as.data.frame(Ravel.phy@otu_table)
identical(rownames(tax), rownames(otu))
tax = tax$Rank7
rownames(otu) = tax


ravel = readRDS("git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/Ravel_Clust3.rds")
ravel = ravel[28:29]

# Order samples by NS
ravel = ravel[order(ravel[,1],decreasing = TRUE),]

# Order samples by best model probs 
#ravel = as.data.frame(t(ravel))
#ravel = ravel[,R_samp]
#ravel = as.data.frame(t(ravel))

identical(rownames(ravel), colnames(otu))
R_Phyla   # cellCycle
otu =  as.matrix(otu[rownames(ravel)])# tt
ravel   # annot



# Firmicutes Phylum
ii = intersect(R_Phyla$p__Firmicutes, rownames(otu))
toPlot = otu[ii,]
p1 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = T,
                             top_annotation = HeatmapAnnotation(df = ravel, 
                                                                col = list(target = c('high' = '#440154FF', 'low' = '#FDE725FF', 'intermediate'='#21908CFF'),
                                                                           cluster = c('C1' = '#21908CFF', 'C2' = '#440154FF', 'C3'='#FDE725FF'))),
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#8DD3C7'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Fusobacteria Phylum
ii = intersect(R_Phyla$p__Fusobacteria, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p2 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FB8072'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Bacteroidetes Phylum
ii = intersect(R_Phyla$p__Bacteroidetes, rownames(otu))
toPlot = otu[ii,]
p3 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FDB462'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# p__Actinobacteria Phylum
ii = intersect(R_Phyla$p__Actinobacteria, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p4 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#BEBADA'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
ht_list = p1 %v% p2 %v% p3 %v% p4
phyla = substr(names(R_Phyla), 4, 50)
draw(ht_list, 
     use_raster = T,
     annotation_legend_list = list(Legend(labels = phyla, 
                                          title = 'Phy 20c 0.08s', 
                                          legend_gp = gpar(fill = c('#8DD3C7', '#FB8072','#FDB462', '#BEBADA')))))
#####
# Family level
#####
Ravel.phy = readRDS("git/BVMetaGenomics/data/Phyloseqs/Ravel_phyloseq.rds") #Reset phyloseq
Ravel.phy = tax_glom(Ravel.phy, taxrank = "Rank7")
Ravel.phy <- filter_taxa(Ravel.phy, function(x) sum(x > 20) > (0.08*length(x)), TRUE)
Ravel.phy = relat.abun(Ravel.phy)

r = as.data.frame(Ravel.phy@tax_table)
r2 = unique(r$Rank5)


R_Phyla = list()
for (i in seq_along(r2)) {
  sub = subset_taxa(Ravel.phy, Rank5 %in% c(r2[i]))
  data = as.data.frame(sub@tax_table[,7])
  nmes = data$Rank7
  R_Phyla [[i]] = nmes
}
names(R_Phyla) =r2

tax = as.data.frame(Ravel.phy@tax_table)
otu = as.data.frame(Ravel.phy@otu_table)
identical(rownames(tax), rownames(otu))
tax = tax$Rank7
rownames(otu) = tax


ravel = readRDS("git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/Ravel_Clust3.rds")
ravel = ravel[28:29]

# Order samples by NS
ravel = ravel[order(ravel[,1],decreasing = TRUE),]

# Order samples by best model probs 
#ravel = as.data.frame(t(ravel))
#ravel = ravel[,R_samp]
#ravel = as.data.frame(t(ravel))

identical(rownames(ravel), colnames(otu))
R_Phyla   # cellCycle
otu =  as.matrix(otu[rownames(ravel)])# tt
ravel   # annot



# Lactobacillaceae Family
ii = intersect(R_Phyla$f__Lactobacillaceae, rownames(otu))
toPlot = otu[ii,]
p1 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = T,
                             top_annotation = HeatmapAnnotation(df = ravel, 
                                                                col = list(target = c('high' = '#440154FF', 'low' = '#FDE725FF', 'intermediate'='#21908CFF'),
                                                                           cluster = c('C1' = '#21908CFF', 'C2' = '#440154FF', 'C3'='#FDE725FF'))),
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#8DD3C7'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Leptotrichiaceae Family
ii = intersect(R_Phyla$f__Leptotrichiaceae, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p2 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FB8072'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Veillonellaceae Family
ii = intersect(R_Phyla$f__Veillonellaceae, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p3 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FDB462'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Prevotellaceae Family
ii = intersect(R_Phyla$f__Prevotellaceae, rownames(otu))
toPlot = otu[ii,]
p4 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#BEBADA'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Peptoniphilaceae Family
ii = intersect(R_Phyla$f__Peptoniphilaceae, rownames(otu))
toPlot = otu[ii,]
p5 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FFFFB3'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Actinomycetaceae Family - hay que sustituir por Coriobacteriaceae -
ii = intersect(R_Phyla$f__Actinomycetaceae, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p6 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#80B1D3'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))

ht_list = p1 %v% p2 %v% p3 %v% p4 %v% p5 %v%p6
phyla = substr(names(R_Phyla), 4, 50)
draw(ht_list, 
     use_raster = T,
     annotation_legend_list = list(Legend(labels = phyla, 
                                          title = 'Fam 20c 0.08s', 
                                          legend_gp = gpar(fill = c('#8DD3C7', '#FB8072','#FDB462', '#BEBADA','#FFFFB3', '#80B1D3')))))
#####
# Genus Level
#####
Ravel.phy = readRDS("git/BVMetaGenomics/data/Phyloseqs/Ravel_phyloseq.rds") #Reset phyloseq
Ravel.phy = tax_glom(Ravel.phy, taxrank = "Rank7")
Ravel.phy <- filter_taxa(Ravel.phy, function(x) sum(x > 20) > (0.08*length(x)), TRUE)
Ravel.phy = relat.abun(Ravel.phy)

r = as.data.frame(Ravel.phy@tax_table)
r2 = unique(r$Rank6)


R_Phyla = list()
for (i in seq_along(r2)) {
  sub = subset_taxa(Ravel.phy, Rank6 %in% c(r2[i]))
  data = as.data.frame(sub@tax_table[,7])
  nmes = data$Rank7
  R_Phyla [[i]] = nmes
}
names(R_Phyla) =r2

tax = as.data.frame(Ravel.phy@tax_table)
otu = as.data.frame(Ravel.phy@otu_table)
identical(rownames(tax), rownames(otu))
tax = tax$Rank7
rownames(otu) = tax


ravel = readRDS("git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/Ravel_Clust3.rds")
ravel = ravel[28:29]

# Order samples by NS
ravel = ravel[order(ravel[,1],decreasing = TRUE),]

# Order samples by best model probs 
#ravel = as.data.frame(t(ravel))
#ravel = ravel[,R_samp]
#ravel = as.data.frame(t(ravel))

identical(rownames(ravel), colnames(otu))
R_Phyla   # cellCycle
otu =  as.matrix(otu[rownames(ravel)])# tt
ravel   # annot



# Lactobacillus Genus
ii = intersect(R_Phyla$g__Lactobacillus, rownames(otu))
toPlot = otu[ii,]
p1 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = T,
                             top_annotation = HeatmapAnnotation(df = ravel, 
                                                                col = list(target = c('high' = '#440154FF', 'low' = '#FDE725FF', 'intermediate'='#21908CFF'),
                                                                           cluster = c('C1' = '#21908CFF', 'C2' = '#440154FF', 'C3'='#FDE725FF'))),
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#8DD3C7'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Sneathia Genus
ii = intersect(R_Phyla$g__Sneathia, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p2 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FB8072'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Dialister Genus
ii = intersect(R_Phyla$g__Dialister, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p3 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FDB462'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Prevotella Genus
ii = intersect(R_Phyla$g__Prevotella, rownames(otu))
toPlot = otu[ii,]
p4 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#BEBADA'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Finegoldia Genus
ii = intersect(R_Phyla$g__Finegoldia, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p5 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FFFFB3'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Peptoniphilus Genus
ii = intersect(R_Phyla$g__Peptoniphilus, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p6 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#80B1D3'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Fannyhessea Genus - Hay que cambiar por Atopobium -
ii = intersect(R_Phyla$g__Fannyhessea, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p7 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(ravel),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#B3DE69'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))

ht_list = p1 %v% p2 %v% p3 %v% p4 %v% p5 %v% p6 %v% p7
phyla = substr(names(R_Phyla), 4, 50)
draw(ht_list, 
     use_raster = T,
     annotation_legend_list = list(Legend(labels = phyla, 
                                          title = 'Gen 20c 0.08s', 
                                          legend_gp = gpar(fill = c('#8DD3C7', '#FB8072','#FDB462', '#BEBADA','#FFFFB3', '#80B1D3','#B3DE69')))))
#####

# Srinivasan
# Phylum level
#####
Srin.phy = readRDS("git/BVMetaGenomics/data/Phyloseqs/Sriniv_phyloseq.rds")
Srin.phy = tax_glom(Srin.phy, taxrank = "Rank7")
Srin.phy <- filter_taxa(Srin.phy, function(x) sum(x > 20) > (0.08*length(x)), TRUE)
Srin.phy = relat.abun(Srin.phy)
s = as.data.frame(Srin.phy@tax_table)
s2 = unique(s$Rank2)

S_Phyla = list()
for (i in seq_along(s2)) {
  sub = subset_taxa(Srin.phy, Rank2 %in% c(s2[i]))
  data = as.data.frame(sub@tax_table[,7])
  nmes = data$Rank7
  S_Phyla [[i]] = nmes
}
names(S_Phyla) =s2

tax = as.data.frame(Srin.phy@tax_table)
otu = as.data.frame(Srin.phy@otu_table)
identical(rownames(tax), rownames(otu))
tax = tax$Rank7
rownames(otu) = tax


sriniv = readRDS("git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/Sriniv_Clust3.rds")
sriniv = sriniv[28:29]
sriniv = sriniv[, c("target", "cluster")]
# Order samples by NS
sriniv = sriniv[order(sriniv[,1],decreasing = TRUE),]

# Order samples by best model probs 
#sriniv = as.data.frame(t(sriniv))
#sriniv = sriniv[,S_samp]
#sriniv = as.data.frame(t(sriniv))

identical(rownames(sriniv), colnames(otu))
S_Phyla   # cellCycle
otu =  as.matrix(otu[rownames(sriniv)])# tt
sriniv   # annot



# Firmicutes Phylum
ii = intersect(S_Phyla$p__Firmicutes, rownames(otu))
toPlot = otu[ii,]
p1 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = T,
                             top_annotation = HeatmapAnnotation(df = sriniv, 
                                                                col = list(target = c('high' = '#440154FF', 'low' = '#FDE725FF', 'intermediate'='#21908CFF'),
                                                                           cluster = c('C1' = '#21908CFF', 'C2' = '#440154FF', 'C3'='#FDE725FF'))),
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#8DD3C7'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Fusobacteria Phylum
ii = intersect(S_Phyla$p__Fusobacteria, rownames(otu))
toPlot = otu[ii,]
p2 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FB8072'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Bacteroidetes Phylum
ii = intersect(S_Phyla$p__Bacteroidetes, rownames(otu))
toPlot = otu[ii,]
p3 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FDB462'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# p__Actinobacteria Phylum
ii = intersect(S_Phyla$p__Actinobacteria, rownames(otu))
toPlot = otu[ii,]
p4 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#BEBADA'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# BVAB2 Phylum
ii = intersect(S_Phyla$p__BVAB2, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p5 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FFFFB3'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))



ht_list = p1 %v% p2 %v% p3 %v% p4 %v% p5
phyla = substr(names(S_Phyla), 4, 50)
draw(ht_list, 
     use_raster = T,
     annotation_legend_list = list(Legend(labels = phyla, 
                                          title = 'Phy 20c 0.08s', 
                                          legend_gp = gpar(fill = c('#8DD3C7', '#FB8072','#FDB462', '#BEBADA', '#FFFFB3')))))
#####
# Family level
#####
Srin.phy = readRDS("git/BVMetaGenomics/data/Phyloseqs/Sriniv_phyloseq.rds")
Srin.phy = tax_glom(Srin.phy, taxrank = "Rank7")
Srin.phy <- filter_taxa(Srin.phy, function(x) sum(x > 20) > (0.08*length(x)), TRUE)
Srin.phy = relat.abun(Srin.phy)
s = as.data.frame(Srin.phy@tax_table)
s2 = unique(s$Rank5)

S_Phyla = list()
for (i in seq_along(s2)) {
  sub = subset_taxa(Srin.phy, Rank5 %in% c(s2[i]))
  data = as.data.frame(sub@tax_table[,7])
  nmes = data$Rank7
  S_Phyla [[i]] = nmes
}
names(S_Phyla) =s2

tax = as.data.frame(Srin.phy@tax_table)
otu = as.data.frame(Srin.phy@otu_table)
identical(rownames(tax), rownames(otu))
tax = tax$Rank7
rownames(otu) = tax


sriniv = readRDS("git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/Sriniv_Clust3.rds")
sriniv = sriniv[28:29]
sriniv = sriniv[, c("target", "cluster")]
# Order samples by NS
sriniv = sriniv[order(sriniv[,1],decreasing = TRUE),]

# Order samples by best model probs 
#sriniv = as.data.frame(t(sriniv))
#sriniv = sriniv[,S_samp]
#sriniv = as.data.frame(t(sriniv))

identical(rownames(sriniv), colnames(otu))
S_Phyla   # cellCycle
otu =  as.matrix(otu[rownames(sriniv)])# tt
sriniv   # annot


# Lactobacillaceae Family
ii = intersect(S_Phyla$f__Lactobacillaceae, rownames(otu))
toPlot = otu[ii,]
p1 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = T,
                             top_annotation = HeatmapAnnotation(df = sriniv, 
                                                                col = list(target = c('high' = '#440154FF', 'low' = '#FDE725FF', 'intermediate'='#21908CFF'),
                                                                           cluster = c('C1' = '#21908CFF', 'C2' = '#440154FF', 'C3'='#FDE725FF'))),
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#8DD3C7'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Leptotrichiaceae Family
ii = intersect(S_Phyla$f__Leptotrichiaceae, rownames(otu))
toPlot = otu[ii,]
p2 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FB8072'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Bifidobacteriaceae Family
ii = intersect(S_Phyla$f__Bifidobacteriaceae, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p3 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FDB462'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
#  Family Actinomycetaceae - Hay que cambiarla por Coriobacteriaceae
ii = intersect(S_Phyla$f__Actinomycetaceae, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p4 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#BEBADA'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Family Veillonellaceae
ii = intersect(S_Phyla$f__Veillonellaceae, rownames(otu))
toPlot = otu[ii,]
p5 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FFFFB3'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Lachnospiraceae Family
ii = intersect(S_Phyla$f__Lachnospiraceae, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p6 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#80B1D3'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Prevotellaceae Family
ii = intersect(S_Phyla$f__Prevotellaceae, rownames(otu))
toPlot = otu[ii,]
p7 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#B3DE69'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Eggerthellaceae Family
ii = intersect(S_Phyla$f__Eggerthellaceae, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p8 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FCCDE5'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# BVAB2 Family
ii = intersect(S_Phyla$f__BVAB2, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p9 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#D9D9D9'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Aerococcaceae Family
ii = intersect(S_Phyla$f__Aerococcaceae, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p10 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#BC80BD'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Peptoniphilaceae Family
ii = intersect(S_Phyla$f__Peptoniphilaceae, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p11 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "#CCEBC5"),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Porphyromonadaceae Family
ii = intersect(S_Phyla$f__Porphyromonadaceae, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p12 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "#FFED6F"),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))





ht_list = p1 %v% p2 %v% p3 %v% p4 %v% p5 %v% p6 %v% p7 %v% p8 %v% p9 %v% p10 %v% p11 %v% p12
phyla = substr(names(S_Phyla), 4, 50)
draw(ht_list, 
     use_raster = T,
     annotation_legend_list = list(Legend(labels = phyla, 
                                          title = 'Fam 20c 0.08s', 
                                          legend_gp = gpar(fill = c('#8DD3C7', '#FB8072','#FDB462', '#BEBADA','#FFFFB3', '#80B1D3',
                                                                    '#B3DE69', '#FCCDE5','#D9D9D9', '#BC80BD', "#CCEBC5", "#FFED6F")))))
#####
# Genus Level
#####
Srin.phy = readRDS("git/BVMetaGenomics/data/Phyloseqs/Sriniv_phyloseq.rds")
Srin.phy = tax_glom(Srin.phy, taxrank = "Rank7")
Srin.phy <- filter_taxa(Srin.phy, function(x) sum(x > 20) > (0.08*length(x)), TRUE)
Srin.phy = relat.abun(Srin.phy)
s = as.data.frame(Srin.phy@tax_table)
s2 = unique(s$Rank6)

S_Phyla = list()
for (i in seq_along(s2)) {
  sub = subset_taxa(Srin.phy, Rank6 %in% c(s2[i]))
  data = as.data.frame(sub@tax_table[,7])
  nmes = data$Rank7
  S_Phyla [[i]] = nmes
}
names(S_Phyla) =s2

tax = as.data.frame(Srin.phy@tax_table)
otu = as.data.frame(Srin.phy@otu_table)
identical(rownames(tax), rownames(otu))
tax = tax$Rank7
rownames(otu) = tax


sriniv = readRDS("git/BVMetaGenomics/data/SpeciesIntersect/Cluster/3C/Sriniv_Clust3.rds")
sriniv = sriniv[28:29]
sriniv = sriniv[, c("target", "cluster")]
# Order samples by NS
sriniv = sriniv[order(sriniv[,1],decreasing = TRUE),]

# Order samples by best model probs 
#sriniv = as.data.frame(t(sriniv))
#sriniv = sriniv[,S_samp]
#sriniv = as.data.frame(t(sriniv))

identical(rownames(sriniv), colnames(otu))
S_Phyla   # cellCycle
otu =  as.matrix(otu[rownames(sriniv)])# tt
sriniv   # annot


# Lactobacillus Genus
ii = intersect(S_Phyla$g__Lactobacillus, rownames(otu))
toPlot = otu[ii,]
p1 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = T,
                             top_annotation = HeatmapAnnotation(df = sriniv, 
                                                                col = list(target = c('high' = '#440154FF', 'low' = '#FDE725FF', 'intermediate'='#21908CFF'),
                                                                           cluster = c('C1' = '#21908CFF', 'C2' = '#440154FF', 'C3'='#FDE725FF'))),
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#8DD3C7'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Sneathia Genus
ii = intersect(S_Phyla$g__Sneathia, rownames(otu))
toPlot = otu[ii,]
p2 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FB8072'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Gardnerella Genus
ii = intersect(S_Phyla$g__Gardnerella, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p3 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FDB462'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
#  Fannyhessea Genus - Hay que cambiarla por Atopobium
ii = intersect(S_Phyla$g__Fannyhessea, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p4 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#BEBADA'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Megasphaera genus
ii = intersect(S_Phyla$g__Megasphaera, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p5 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FFFFB3'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Lachnospiraceae Family
ii = intersect(S_Phyla$g__, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p6 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#80B1D3'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Prevotella Genus
ii = intersect(S_Phyla$g__Prevotella, rownames(otu))
toPlot = otu[ii,]
p7 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#B3DE69'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Eggerthella Genus
ii = intersect(S_Phyla$g__Eggerthella, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p8 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#FCCDE5'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# BVAB2 Genus
ii = intersect(S_Phyla$g__BVAB2, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p9 = ComplexHeatmap::Heatmap(toPlot,
                             show_row_names = F, show_column_names = F,
                             column_order = rownames(sriniv),
                             cluster_rows = F,
                             show_heatmap_legend = F,
                             right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#D9D9D9'),
                                                                               labels_gp = gpar(col = 'white', fontsize = 5))))
# Aerococcus Genus
ii = intersect(S_Phyla$g__Aerococcus, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p10 = ComplexHeatmap::Heatmap(toPlot,
                              show_row_names = F, show_column_names = F,
                              column_order = rownames(sriniv),
                              cluster_rows = F,
                              show_heatmap_legend = F,
                              right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = '#BC80BD'),
                                                                                labels_gp = gpar(col = 'white', fontsize = 5))))
# Parvimonas Genus
ii = intersect(S_Phyla$g__Parvimonas, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p11 = ComplexHeatmap::Heatmap(toPlot,
                              show_row_names = F, show_column_names = F,
                              column_order = rownames(sriniv),
                              cluster_rows = F,
                              show_heatmap_legend = F,
                              right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "#CCEBC5"),
                                                                                labels_gp = gpar(col = 'white', fontsize = 5))))
# Dialister Genus
ii = intersect(S_Phyla$g__Dialister, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p12 = ComplexHeatmap::Heatmap(toPlot,
                              show_row_names = F, show_column_names = F,
                              column_order = rownames(sriniv),
                              cluster_rows = F,
                              show_heatmap_legend = F,
                              right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "#FFED6F"),
                                                                                labels_gp = gpar(col = 'white', fontsize = 5))))
# Porphyromonas Genus
ii = intersect(S_Phyla$g__Porphyromonas, rownames(otu))
toPlot = otu[ii,]
toPlot = as.data.frame(t(toPlot))
rownames(toPlot)= ii
toPlot = as.matrix(toPlot)
p13 = ComplexHeatmap::Heatmap(toPlot,
                              show_row_names = F, show_column_names = F,
                              column_order = rownames(sriniv),
                              cluster_rows = F,
                              show_heatmap_legend = F,
                              right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "#FFE5EE"),
                                                                                labels_gp = gpar(col = 'white', fontsize = 5))))




ht_list = p1 %v% p2 %v% p3 %v% p4 %v% p5 %v% p6 %v% p7 %v% p8 %v% p9 %v% p10 %v% p11 %v% p12 %v% p13
# We have to add genus name for BVAB1
names(S_Phyla)[6] = "g__BVAB1"
phyla = substr(names(S_Phyla), 4, 50)
draw(ht_list, 
     use_raster = T,
     annotation_legend_list = list(Legend(labels = phyla, 
                                          title = 'Gen 20c 0.08s', 
                                          legend_gp = gpar(fill = c('#8DD3C7', '#FB8072','#FDB462', '#BEBADA','#FFFFB3', '#80B1D3',
                                                                    '#B3DE69', '#FCCDE5','#D9D9D9', '#BC80BD', "#CCEBC5", "#FFED6F", "#FFE5EE")))))
#####