library(argparse)
library(velocyto.R)
library(Seurat)
library(gdata)
#library(tidyverse)
library(Matrix)
#library(clustree)
#library(pheatmap)
#library(foreach)
library(SeuratWrappers)
#library(UpSetR)
set.seed(100)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--seuratobj1', dest='seurat1', metavar='N', type="character", nargs='+')
    parser$add_argument('--seuratobj2', dest='seurat2', default='')
    parser$add_argument('--label', dest='label', default='')

    parser$add_argument('--species', dest='species', type='character', default='human')
    args = parser$parse_args()
    args
}

args = get_args()
if (length(args$seurat1) == 0 || args$label == '') {
    cat('empty argument, exit..')
    q()
}

#args$seurat1 = '../results/seurat/20190418_MAST111_5Kcells_hg19_seurat_obj_tumors.rds'
#args$seurat2 = '/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds'
#args$label = 'MAST111'

seurat1 = readRDS(args$seurat1) ## mine

if (args$seurat2 != '') {
   seurat2 = readRDS(args$seurat2) ## sara
   print(length(intersect(rownames(seurat2@meta.data), rownames(seurat1@meta.data))))
}

if (args$seurat2 != '') {
    seurat1 = subset(seurat1, cells=rownames(seurat2@meta.data))
    seurat2 = subset(seurat2, cells=rownames(seurat1@meta.data))
    print(dim(seurat1))
    print(dim(seurat2))
}

seurat1.de = FindAllMarkers(seurat1, only.pos=T, min.pct=0.1,
                            test.use='MAST', min.diff.pct=0.1,
                            random.seed=100, logfc.threshold = 0.2)

seurat1.de$enrichment = seurat1.de$pct.1 - seurat1.de$pct.2
system('mkdir -p ../results/seurat_v6')
write.table(seurat1.de, file=paste0('../results/seurat_v6/', args$label, '_SCT_res0.8.xls'), sep='\t', quote=F, col.names=NA)

#pdf(paste0('../results/seurat_v6/', args$label, '_v6_comparison.pdf'), width=22, height=6)
#Idents(seurat1) = seurat1@meta.data$seurat_clusters = seurat2$RNA_snn_res.0.8
#p1 = DimPlot(seurat2, reduction='umap', label=T, group.by='RNA_snn_res.0.8') ## Sara umap with sara cluster
#p2 = DimPlot(seurat1, reduction='umap', label=T, group.by='seurat_clusters') ## Alvin umap with sara cluster
#p3 = DimPlot(seurat1, reduction='umap', label=T, group.by='SCT_snn_res.0.8') ## Alvin umap with Alvin cluster
#seurat1@reductions$umap = seurat2@reductions$umap
#p4 = DimPlot(seurat1, reduction='umap', label=T, group.by='SCT_snn_res.0.8') ## Sara umap with Alvin cluster
#p = CombinePlots(plots=list(p1, p4, p2, p3), ncol=4)
#print(p)
#dev.off()
#
#library(ggplot2)
#
#pdf(paste0('../results/seurat_v6/', args$label, '_v6_heatmap.pdf'), width=22, height=30)
#par(mar=c(3,5,8,3))
##heat = DoHeatmap(subset(seurat1, downsample = 100), features = seurat1.de$gene, size = 3) + theme(
##do not downsample cells
#heat = DoHeatmap(seurat1, features = seurat1.de$gene, size = 5) + theme(
#  panel.background = element_rect(fill = "white"),
#  plot.margin = margin(2, 3, 5, 1, "cm"),
#  plot.background = element_rect(
#    fill = "white",
#    colour = "white",
#    size = 1
#  )
#)
#print(heat)
#dev.off()

