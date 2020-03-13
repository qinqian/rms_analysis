library(Seurat)
library(cellassign)

primary1 = readRDS('../data/seurat_obj/all_primary_patient.RDS')

data(example_TME_markers)
library(SingleCellExperiment)

primary1.sce = as.SingleCellExperiment(primary1)

sizeFactors(primary1.sce) <- colSums(assay(primary1.sce))

marker_mat = marker_list_to_mat(example_TME_markers$symbol)

sce_marker <- primary1.sce[intersect(rownames(marker_mat), rownames(primary1.sce)),]

s = sizeFactors(sce_marker)
sce_marker = sce_marker[, s > 0]
s = s[s > 0]

cas <- cellassign(exprs_obj = sce_marker,
                  marker_gene_info = marker_mat[intersect(rownames(marker_mat), rownames(primary1.sce)),],
                  s = s)

saveRDS(cas, '../data/all_primary_train_cellassignmodel.rds')

primary1.objann = subset(primary1, cells=colnames(sce_marker))
primary1.objann$seurat_clusters = cas$cell_type

library(glue)
pdf(glue('../figures/all_primary_integrated_cellassign.pdf'), width=15, height=5)
p1 = DimPlot(primary1, reduction='umap', group.by='orig.ident')
p2 = DimPlot(primary1.objann, reduction='umap', group.by='seurat_clusters')
print(CombinePlots(plots=list(p1, p2)))
dev.off()
