library(Seurat)
target = readRDS('../results/seurat_sara/20191031_MSK74711_seurat-object.rds')
prior_6 = c("CD44", "CXCR4", "LRRN1", "TSPAN33", "THY1", "CHODL")

surface = scan('surface_protein.txt', what='')
muscle = scan('Muscle.txt', what='')
emt = scan('EMT.txt', what='')

muscle = intersect(surface, muscle)
emt = intersect(surface, emt)

cols = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1,
                  header=T,
                  check.names=F, stringsAsFactors=F)

ann = unlist(as.vector(cols['MSK74711',,drop=T]))
ann = ann[ann!='']
levels(target@meta.data$RNA_snn_res.0.8) = ann

pdf("Claudia_RED_Surface.pdf", width=16, height=11)
CombinePlots(plots=list(DimPlot(target, group.by='RNA_snn_res.0.8'),
                        FeaturePlot(target, features=prior_6)), ncol=2)
dev.off()

pdf("Claudia_MUSCLESurface_Surface.pdf", width=20, height=9)
CombinePlots(plots=list(DimPlot(target, group.by='RNA_snn_res.0.8'),
                        FeaturePlot(target, features=muscle)), ncol=2)
dev.off()

pdf("Claudia_EMTSurface_Surface.pdf", width=24, height=9)
CombinePlots(plots=list(DimPlot(target, group.by='RNA_snn_res.0.8'),
                        FeaturePlot(target, features=emt, ncol=3)), ncol=2)
dev.off()

qpcr = scan('qpcr.txt', what='')

pdf("Claudia_QPCR.pdf", width=28, height=28)
## CombinePlots(plots=list(DimPlot(target, group.by='RNA_snn_res.0.8'),
##                         FeaturePlot(target, features=qpcr, ncol=5)), ncol=2)
FeaturePlot(target, features=qpcr, ncol=5)
dev.off()


library(Seurat)
cols = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1,
                  header=T,
                  check.names=F, stringsAsFactors=F)

targets = list(readRDS('../results/seurat_sara/20191031_MSK74711_seurat-object.rds'),
               readRDS('../results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds'),
               readRDS('../results/seurat_sara/MAST118_seurat-object.rds'))

label = c("MSK74711", "MSK72117", "MAST118")
prior_6 = c("MKI67", "BNIP3", "MMP2", "LRRN1", "IFI6")
n = 1

library(patchwork)
pdf('three_tumor_biomarkers.pdf', height=7.5, width=22)
for (target in targets) {
    ann = unlist(as.vector(cols[label[n],,drop=T]))
    ann = ann[ann!='']
    levels(target@meta.data$RNA_snn_res.0.8) = ann
    print(DimPlot(target, group.by='RNA_snn_res.0.8')+FeaturePlot(target, slot='data', features=prior_6, ncol=3))
    n = n+1
}
dev.off()
