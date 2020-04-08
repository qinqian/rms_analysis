library(Seurat)

tumor21 = readRDS('../results/seurat_intersect_velocity/Tumor21_seu.rds')
tumor21 = FindClusters(tumor21, resolution=1.8)

tumor21.de = FindAllMarkers(tumor21, only.pos=T, min.pct=0.1,
                            test.use='MAST', min.diff.pct=0.1,
                            random.seed=100, logfc.threshold = 0.2)
tumor21.de$enrichment = tumor21.de$pct.1 - tumor21.de$pct.2

saveRDS(tumor21, paste0('../results/seurat_v6/', 'Tumor21_recluster1.8.rds'))

write.table(tumor21.de, file=paste0('../results/seurat_v6/', 'Tumor21_recluster1.8', '_SCT_res0.8.xls'),
            sep='\t', quote=F, col.names=NA)

system("Rscript annotate_celltypes.R --seuratobj ../results/seurat_v6/Tumor21_recluster1.8.rds --label Tumor21_recluster --species fish")

tumor22 = readRDS('../results/seurat_intersect_velocity/Tumor22_seu.rds')
tumor22 = FindClusters(tumor22, resolution=1.8)
tumor22.de = FindAllMarkers(tumor22, only.pos=T, min.pct=0.1,
                            test.use='MAST', min.diff.pct=0.1,
                            random.seed=100, logfc.threshold = 0.2)
tumor22.de$enrichment = tumor22.de$pct.1 - tumor22.de$pct.2

saveRDS(tumor22, paste0('../results/seurat_v6/', 'Tumor22_recluster1.8.rds'))

write.table(tumor22.de, file=paste0('../results/seurat_v6/', 'Tumor22_recluster1.8', '_SCT_res0.8.xls'),
            sep='\t', quote=F, col.names=NA)

system("Rscript annotate_celltypes.R --seuratobj ../results/seurat_v6/Tumor22_recluster1.8.rds --label Tumor22_recluster --species fish")
