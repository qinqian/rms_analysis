library(Seurat)

## rh = readRDS('../results/seurat_intersect_velocity/20190617_RH74-10cells_seu.rds')

## mast111 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds')
## mine = readRDS('../results/seurat/MAST111_seurat_obj_tumors.rds')
## ## mast111 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds')

## mast111 = FindClusters(mast111, resolution=0.8)
## mine = subset(mine, cells=rownames(mast111@meta.data))
## mast111 = subset(mast111, cells=rownames(mine@meta.data))

## table(mine$seurat_clusters, mast111$RNA_snn_res.0.8)
## Idents(mine) = mine@meta.data$seurat_clusters = mast111$RNA_snn_res.0.8

## pdf('test.pdf', width=18, height=10)
## p1 = DimPlot(mast111, reduction='tsne', label=T, group.by='RNA_snn_res.0.8')
## p2 = DimPlot(mast111, reduction='umap', label=T, group.by='RNA_snn_res.0.8')
## p3 = DimPlot(mine, reduction='tsne', label=T, group.by='seurat_clusters')
## p4 = DimPlot(mine, reduction='umap', label=T, group.by='seurat_clusters')
## p5 = DimPlot(mine, reduction='tsne', label=T, group.by='SCT_snn_res.0.45')
## p6 = DimPlot(mine, reduction='umap', label=T, group.by='SCT_snn_res.0.8')
## p7 = DimPlot(mine, reduction='tsne', label=T, group.by='SCT_snn_res.0.45')
## p8 = DimPlot(mine, reduction='tsne', label=T, group.by='SCT_snn_res.0.8')
## CombinePlots(plots=list(p1, p3, p5, p7, p2, p4, p6, p8), ncol=4)
## dev.off()


fishes = readRDS('../results/final_seurat_normalize_cluster_obj.RDS')
tumor24 = readRDS('../results/seurat/Tumor24_seurat_obj_tumors.rds')

orig = subset(fishes[[1]], ident=c(0, 3))


## tumor24 = readRDS('../results/seurat_intersect_velocity/Tumor24_seu.rds')

tumor24.orig = subset(fishes[[1]], cells=rownames(tumor24@meta.data))
tumor24 = subset(tumor24, cells=rownames(tumor24.orig@meta.data))

source('functions.R')
tumor24.orig = recluster.withtree(tumor24.orig)
tumor24.orig = FindClusters(tumor24.orig, resolution=0.8)

tumor24 = recluster.withtree(tumor24)
tumor24 = FindClusters(tumor24, resolution=0.8)

table(tumor24.orig@meta.data$seurat_clusters, tumor24@meta.data$seurat_clusters)
