library(argparse)
library(org.Hs.eg.db)
library(velocyto.R)
library(Seurat)
library(gdata)
library(tidyverse)
library(Matrix)
library(clustree)
library(pheatmap)
library(foreach)
library(SeuratWrappers)
library(UpSetR)
library(cellAlign)
set.seed(100)

## args$seurat = '../results/seurat_intersect_velocity/20190624_seurat-object_MAST111_seu.rds'
## args$vel = '../results/seurat_intersect_velocity/20190624_seurat-object_MAST111_vel.rds'
## args$label = '20190624_seurat-object_MAST111_seu.rds'

velocity = Tool(object=vel, slot = "RunVelocity")
colortab = read.table('color_table.xls', sep='\t', header=T)
## args$clusterlabel = c(1, 2, 3, 1, 4, 5, 6, 3, 7)
## args$clusterlabel = 'MAST111'
metalabels = c("GROUND", "Hypoxia", "EMT",  "G1S", "UNASSIGNED", "G2M",  "MUSCLE", "INTERFERON", "PROLIF", "Histones")
metacolors = c(rgb(166, 166, 166, maxColorValue = 255),
               rgb(241, 149, 69, maxColorValue = 255),
               rgb(103, 35,  102, maxColorValue = 255),
               rgb(52, 101, 252, maxColorValue = 255),
               rgb(242, 242, 242, maxColorValue = 255),
               rgb(52, 101, 252, maxColorValue = 255),
               rgb(233, 63,  51, maxColorValue = 255),
               rgb(65, 129,  7, maxColorValue = 255),
               rgb(52, 101, 252, maxColorValue = 255),
               rgb(253, 247, 49, maxColorValue = 255))

cluster = as.character(colortab[,args$clusterlabel])
labels = na.omit(metalabels[match(cluster, metalabels)])
colors = na.omit(metacolors[match(cluster, metalabels)])

ident.colors = colors
print(ident.colors)
print(levels(obj@meta.data$RNA_snn_res.0.8))
names(ident.colors) <- levels(obj@meta.data$RNA_snn_res.0.8)
cell.colors <- ident.colors[obj@meta.data$RNA_snn_res.0.8]
names(cell.colors) <- colnames(vel)

embedvel = Embeddings(object=vel, reduction = "umap")
embed = Embeddings(object=obj, reduction = "umap")
rownames(embed) = rownames(embedvel)
## embed = Embeddings(object=vel, reduction = "umap")

pdf(paste0('../results/', name, '_velocity_tumoronly.pdf'), width=16, height=9)
## par(mfrow=c(1, 2), font=2, cex=1.3)
par(font=2, cex=1.3, xpd=T, mar=c(4, 4, 3, 15))
show.velocity.on.embedding.cor(emb = embed,
    			       vel = velocity,
    			       n = 200, scale = "sqrt", 
    			       cell.colors = ac(cell.colors, alpha = 0.5),
    			       cex=1, n.cores=6,
    			       arrow.scale = 1.2, 
    			       show.grid.flow = TRUE, 
    			       min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = F, cell.border.alpha = 0.3)
legend(max(embed[,1])*1.15, max(embed[,2])*1.0, legend=paste0(names(ident.colors), '_', labels), pch=19, col=ident.colors, bty = "n")
title(args$label)
## show.velocity.on.embedding.cor(emb = Embeddings(object=vel, reduction = "tsne"), n.cores=10,
##     			       vel = Tool(object=vel, slot = "RunVelocity"), 
##     			       n = 200, scale = "sqrt", 
##     			       cell.colors = ac(cell.colors, alpha = 0.5),
##     			       cex=1,
##     			       arrow.scale = 2, 
##     			       show.grid.flow = TRUE, 
##     			       min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
## legend('topright', legend=levels(obj@meta.data$RNA_snn_res.0.8), pch=19, col=ident.colors)
dev.off()
