library(Seurat)
annotation <- read.delim('../final_annotations/primary_clusters.txt', sep='\t',
                         row.names=1)

## primary <- lapply(c('../results/seurat_sara/20082_hg19_premrna_seurat-object.rds',
##                     '../results/seurat_sara/20696_seurat-object.rds',
##                     '../results/seurat_sara/21202_hg19_premrna_seurat-object.rds',
##                     '../results/seurat_sara/29806_hg19_premrna_seurat-object.rds'), readRDS)

primary.tumors <- lapply(c('20082_recluster2_tumor_only.rds',
                           '../figures/20696_hg19_tumoronly_res0.8_umap.rds',
                           '../figures/21202_hg19_premrna_tumoronly_res0.8_umap.rds',
                           '../figures/29806_hg19_premrna_tumoronly_res0.8_umap.rds'), readRDS)

labels <- c("20082", "20696", "21202", "29806")

## primary <- lapply(c('../results/seurat/Tumor21_unfilter_seurat_obj_tumors.rds',
##                     '../results/seurat/Tumor22_unfilter_seurat_obj_tumors.rds',
##                     '../results/seurat/Tumor24_unfilter_seurat_obj_tumors.rds'), readRDS)
## primary.tumors <- lapply(c('../results/seurat_v6/Tumor21_recluster1.8.rds',
##                            '../results/seurat_v6/Tumor22_recluster1.8.rds',
##                            '../results/seurat_intersect_velocity/Tumor24_seu.rds'), readRDS)
## labels <- c("Tumor21", "Tumor22", "Tumor24")
pdxs = Sys.glob('../data/seurat_obj/*rds')[1:10]
cols = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1,
                  header=T,
                  check.names=F, stringsAsFactors=F)
pdxs = c(pdxs, '../results/seurat_sara/20191031_MSK74711_seurat-object.rds',
         '../results/seurat_sara/MAST118_seurat-object.rds',
         '../results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds',
         '../results/seurat_sara/MAST139_1cells_seurat-object.rds')
pdxs.objs = lapply(pdxs, readRDS)

labels2 = unlist(lapply(pdxs.objs, function(x) {
    levels(x$orig.ident[1])
}))
labels2[9] = 'MAST85-1'
labels2[11] = 'MSK74711'
labels2[13] = 'MSK72117'
labels2[14] = 'MAST139-1'
labels2[10] = 'RH74-10'

genes = unique(unlist(lapply(pdxs.objs, rownames)))
genes2 = unique(unlist(lapply(primary.tumors, rownames)))

genes = intersect(genes, genes2)

PI3K = scan("PI3K.txt", what='')

target = intersect(PI3K, genes)

library(ggplot2)
library(patchwork)

for (gene in target){
    try({
    for (i in seq_along(labels2)) {
        if (i==1)
            p = FeaturePlot(pdxs.objs[[i]], features=gene, reduction='umap')+ggtitle(labels2[i])
        else
            p = p+FeaturePlot(pdxs.objs[[i]], features=gene, reduction='umap')+ggtitle(labels2[i])
    }
    pdf(paste0(gene, '_pdx.pdf'), width=28, height=28)
    print(p)
    dev.off()
    for (i in seq_along(labels)) {
        if (i==1)
            p = FeaturePlot(primary.tumors[[i]], features=gene, reduction='umap')+ggtitle(labels[i])
        else
            p = p+FeaturePlot(primary.tumors[[i]], features=gene, reduction='umap')+ggtitle(labels[i])
    }
    pdf(paste0(gene, '_primary.pdf'), width=14, height=14)
    print(p)
    dev.off()
    })
}
