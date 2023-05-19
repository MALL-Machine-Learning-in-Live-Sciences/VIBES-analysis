# Circos plot
# ========

library(circlize)
library(RColorBrewer)
library(grDevices)
lista.betas <- readRDS(file = 'figures/data/betas.rds')
taxonomia <- readRDS(file = 'figures/data/tax_table.rds')
lista.especies <- setdiff(unique(c(rownames(lista.betas$N),rownames(lista.betas$IDN),rownames(lista.betas$IDD),rownames(lista.betas$D))), '(Intercept)')
df.raw <- taxonomia
df.raw$betas <- 0
i=4
df.betas <- as.data.frame(as.matrix(lista.betas[[i]]))
df.betas$especie <- rownames(df.betas)
colnames(df.betas)[1]<- 'betas'
df.raw$betas <- df.betas$betas[match(df.raw$Species,df.betas$especie)]
df.raw[is.na(df.raw)] <- 0
df.raw$betas <- abs(df.raw$betas)
df <- df.raw[,c('Species','Genus','betas')]
rownames(df)<-NULL
nm <- unique(c(df$Species,df$Genus))
colores <- colorRampPalette(brewer.pal(8, "Set2"))(length(nm))
grid.col <- data.frame(nm,colores)
col.pal <- tibble::deframe(grid.col)
group = structure(c(rep('Species',20),
                    rep('Genus',length(nm)-20)), names = nm)
# tiff(paste0('~/git/pruebasCodigo/cluster_',i,'.tif'),
#    width = 3.74,height = 3.74,units = 'in',res = 200)
svg(paste0('figures/plots/fig2/circos_cluster_',i,'.svg'),
     width = 12, 
     height = 12,
     pointsize = 7)
circos.clear()
circos.par(start.degree=90)
circos.par("canvas.xlim" = c(-2.5, 2.5))  # Ajustar el espacio en blanco horizontal
circos.par("canvas.ylim" = c(-2.5, 2.5))
chordDiagram(df,grid.col = col.pal,
             annotationTrack = c('grid'),
             annotationTrackHeight = c(0.03, 0.01),
             group = group,big.gap = 30, small.gap = 0.5
             #col = colorRampPalette(brewer.pal(8, "Set2"))(20)
)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  if(sector.name %in% c(unique(df$Genus),unique(df$Species))){
    circos.text(CELL_META$xcenter,
                ylim[1] + cm_h(2),
                sector.name,
                facing = "clockwise",
                niceFacing = TRUE,
                adj = c(0, 0.5),
                cex = 3,
                #col=grid.col[sector.name],
                #font = 2
    )
  }
}, bg.border = NA)
#circos.export(file = paste0('~/git/pruebasCodigo/cluster_',i,'.svg'), type = "svg")
dev.off()
