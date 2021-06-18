library(Seurat)
library(patchwork)

rd = readRDS('../results/seurat_sara/RD_seurat-object.rds')
rh = readRDS('../results/seurat_sara/RH41_seurat-object.rds')

pdf('RD_RH41_MYCMYCN.pdf', width=18, height=8.5)
p1.0=DimPlot(rd, group.by = 'RNA_snn_res.0.8')
p1=FeaturePlot(rd, 'MYC')
p1.1=FeaturePlot(rd, 'MYCN')
p2.0=DimPlot(rh, group.by = 'RNA_snn_res.0.8')
p2=FeaturePlot(rh, 'MYC')
p2.1=FeaturePlot(rh, 'MYCN')
print(p1.0+p1+p1.1+p2.0+p2+p2.1+plot_layout(ncol=3))
dev.off()

pdxs = Sys.glob('../data/seurat_obj/*rds')[1:10]
cols = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1,
                  header=T,
                  check.names=F, stringsAsFactors=F)
pdxs = c(pdxs, '../results/seurat_sara/20191031_MSK74711_seurat-object.rds',
         '../results/seurat_sara/MAST118_seurat-object.rds',
         '../results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds',
         '../results/seurat_sara/MAST139_1cells_seurat-object.rds')
pdxs.objs = lapply(pdxs, readRDS)
labels = unlist(lapply(pdxs.objs, function(x) {
    levels(x$orig.ident[1])
}))
labels[9] = 'MAST85-1'
labels[11] = 'MSK74711'
labels[13] = 'MSK72117'
labels[14] = 'MAST139-1'
labels[10] = 'RH74-10'

pdf('MAST111_MAST139_MYCMYCN.pdf', width=18, height=8.5)
p1.0=DimPlot(pdxs.objs[[1]], group.by = 'RNA_snn_res.0.8')
p1=FeaturePlot(pdxs.objs[[1]], 'MYC')
p1.1=FeaturePlot(pdxs.objs[[1]], 'MYCN')
p2.0=DimPlot(pdxs.objs[[2]], group.by = 'RNA_snn_res.0.8')
p2=FeaturePlot(pdxs.objs[[2]], 'MYC')
p2.1=FeaturePlot(pdxs.objs[[2]], 'MYCN')
print(p1.0+p1+p1.1+p2.0+p2+p2.1+plot_layout(ncol=3))
dev.off()


pdf('MAST85_MAST74711_MYCMYCN.pdf', width=18, height=8.5)
p1.0=DimPlot(pdxs.objs[[5]], group.by = 'RNA_snn_res.0.8')
p1=FeaturePlot(pdxs.objs[[5]], 'MYC')
p1.1=FeaturePlot(pdxs.objs[[5]], 'MYCN')
p2.0=DimPlot(pdxs.objs[[11]], group.by = 'RNA_snn_res.0.8')
p2=FeaturePlot(pdxs.objs[[11]], 'MYC')
p2.1=FeaturePlot(pdxs.objs[[11]], 'MYCN')
print(p1.0+p1+p1.1+p2.0+p2+p2.1+plot_layout(ncol=3))
dev.off()


library(ggplot2)
pdf('all_PDX_MYCMYCN.pdf', width=15, height=8.5*8)
n = 1
for (i in seq_along(pdxs.objs)) {
    p1.0=DimPlot(pdxs.objs[[i]], group.by = 'RNA_snn_res.0.8')+ggtitle(labels[i])
    p1=FeaturePlot(pdxs.objs[[i]], 'MYC')+ggtitle("MYC")
    p1.1=FeaturePlot(pdxs.objs[[i]], 'MYCN')+ggtitle("MYCN")
    if (n==1)
        p = p1.0+p1+p1.1
    else
        p = p+p1.0+p1+p1.1
    n <- n+1
}
print(p+plot_layout(ncol=3))
dev.off()
