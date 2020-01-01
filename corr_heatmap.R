library(Seurat)
library(RColorBrewer)
library(ggcorrplot)
library(pheatmap)

cell.num = list()
for (s in Sys.glob("../results/seurat_sara/*.rds")) {
    seu <- readRDS(s)
    index = order(seu@meta.data$RNA_snn_res.0.8)
    d = sub('.rds', '_SCT_res0.8.xls', s)
    l = sub('_seurat-object_SCT_res0.8.xls', '', d)
    cell.num[[basename(l)]] = data.frame(table(seu@meta.data$RNA_snn_res.0.8))
    f = sub('.rds', '.png', s)
    diffgene = read.table(d)
    seu.data <- as.matrix(seu$RNA@data)
    ## seu.data <- seu.data[seu@assays$RNA@var.features, ]
    seu.data <- seu.data[intersect(seu@assays$RNA@var.features, diffgene$gene), ]
    ## seu.data <- seu.data[diffgene$gene, ]
    seu.cor <- cor(seu.data, method='pearson')
    sum(colnames(seu.data) == rownames(seu@meta.data))
    plot_data = seu.cor[index, index]
    ## png('test.png', width=5, height=5, res=500, units='in')
    ## pdf('test.pdf', width=5, height=5)
    ## image(1:2000, 1:2000, plot_data[1:2000, 1:2000],
    ##       col=colorRampPalette(c("navy", "white", "firebrick3"))(50),
    ##       breaks=seq(-0.1, 1, length.out=51), axes=F)
    ## dev.off()
    annrow = data.frame(cluster=sort(seu@meta.data$RNA_snn_res.0.8))
    rownames(annrow) = rownames(seu.cor)[index]
    pheatmap(plot_data, annotation_row=annrow,annotation_col=annrow, 
             show_rownames = F, show_colnames = F, filename=f, width=8.5, height=8, res=600,
             breaks=seq(-0.1, 1, length.out=51),
             cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("blue", "white", "red"))(50))
}

cell.nums = Reduce(function(x, y) {merge(x, y, by='Var1',all=T)}, cell.num)

colnames(cell.nums) <- c("Cluster", names(cell.num))

write.table(cell.nums, file='../results/seurat_sara/V4_New_sample_cluster_cell_number.xls', sep='\t', quote=F, row.names=F)
