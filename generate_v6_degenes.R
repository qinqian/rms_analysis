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
set.seed(100)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--seuratobj1', dest='seurat1', metavar='N', type="character", nargs='+')
    parser$add_argument('--seuratobj2', dest='seurat2', metavar='N', type="character", nargs='+')
    parser$add_argument('--label', dest='label', default='')

    parser$add_argument('--species', dest='species', type='character', default='human')
    args = parser$parse_args()
    args
}

args = get_args()
if (length(args$seurat) == 0 || args$vel == '' || args$label == '') {
    cat('empty argument, exit..')
    q()
}

## args$seurat1 = '../results/seurat/20190418_MAST111_5Kcells_hg19_seurat_obj_tumors.rds'
## args$seurat2 = '/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds'
## args$label = 'MAST111'

seurat1 = readRDS(args$seurat1) ## mine
seurat2 = readRDS(args$seurat2) ## sara

seurat1 = subset(seurat1, cells=rownames(seurat2@meta.data))
seurat2 = subset(seurat2, cells=rownames(seurat1@meta.data))

seurat1.de = FindAllMarkers(seurat1, only.pos=T, min.pct=0.1,
                            test.use='MAST',
                            random.seed=100, logfc.threshold = 0.2)

system('mkdir -p ../results/seurat_v6')
write.table(seurat1.de, file=paste0('../results/seurat_v6/', args$label, '_SCT_res0.8.xls'), sep='\t', quote=F)

pdf(paste0('../results/seurat_v6/', args$label, '_v6_comparison.pdf'), width=22, height=6)
Idents(seurat1) = seurat1@meta.data$seurat_clusters = seurat2$RNA_snn_res.0.8
p1 = DimPlot(seurat2, reduction='umap', label=T, group.by='RNA_snn_res.0.8') ## Sara umap with sara cluster
p2 = DimPlot(seurat1, reduction='umap', label=T, group.by='seurat_clusters') ## Alvin umap with sara cluster
p3 = DimPlot(seurat1, reduction='umap', label=T, group.by='SCT_snn_res.0.8') ## Alvin umap with Alvin cluster
seurat1@reductions$umap = seurat2@reductions$umap
p4 = DimPlot(seurat1, reduction='umap', label=T, group.by='SCT_snn_res.0.8') ## Sara umap with Alvin cluster
CombinePlots(plots=list(p1, p4, p2, p3), ncol=4)
dev.off()
