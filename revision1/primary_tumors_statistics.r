library(Seurat)

res = c('/PHShome/qq06/langenau/projects/01_sc_rms/results/seurat_sara/20082_hg19_premrna_seurat-object.rds', 
        '/PHShome/qq06/langenau/projects/01_sc_rms/results/seurat_sara/20696_seurat-object.rds',
        '/PHShome/qq06/langenau/projects/01_sc_rms/results/seurat_sara/21202_hg19_premrna_seurat-object.rds',
        '/PHShome/qq06/langenau/projects/01_sc_rms/results/seurat_sara/29806_hg19_premrna_seurat-object.rds')

res = lapply(res, readRDS)

labels = gsub('_hg19_premrna', '', unlist(lapply(res, function(x) {
    levels(x$orig.ident[1])
})))

res.merge = Reduce(merge, res)

print(mean(res.merge@meta.data$nFeature_RNA))
print(sd(res.merge@meta.data$nFeature_RNA))
