##1. remove blood cell clusters
library(tidyverse)
library(ggplot2)
library(foreach)
library(Seurat)
library(SeuratWrappers)
library(velocyto.R)
library(clustree)
library(fgsea)
library(ReactomePA)
library(clusterProfiler)
library(reactome.db)
source('functions.R')

## 2.recluster using intersection of velocity and seurat preprocessed datasets
samples <- readRDS('../results/final_seurat_normalize_cluster_obj.RDS')
sample1.obj <- samples[[1]]
sample2.merge.obj <- samples[[2]]
sample3.merge.obj <- samples[[3]]

sample1.vel <- as.Seurat(ReadVelocity(file = "/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color_v2/velocyto/Tumor24_zebrafish_with_orf_color_v2.loom"))
sample2_2.vel <- as.Seurat(ReadVelocity(file = "/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2/velocyto/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2.loom"))
sample2.vel <- as.Seurat(ReadVelocity(file = "/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish_with_orf_color_v2/velocyto/Tumor21_zebrafish_with_orf_color_v2.loom"))
sample3_2.vel <- as.Seurat(ReadVelocity(file = "/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2/velocyto/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2.loom"))
sample3.vel <- as.Seurat(ReadVelocity(file = "/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish_with_orf_color_v2/velocyto/Tumor22_zebrafish_with_orf_color_v2.loom"))

sample2.vel <- merge(sample2.vel, sample2_2.vel, add.cell.ids = c("sort", 'bulk'), project = 'RMS21_zebrafish')
sample3.vel <- merge(sample3.vel, sample3_2.vel, add.cell.ids = c("sort", 'bulk'), project = 'RMS22_zebrafish')

sample1.vel[["percent.mt"]] <- PercentageFeatureSet(sample1.vel, pattern='^mt-')
sample2.vel[["percent.mt"]] <- PercentageFeatureSet(sample2.vel, pattern='^mt-')
sample3.vel[["percent.mt"]] <- PercentageFeatureSet(sample3.vel, pattern='^mt-')

sample1.vel <- subset(sample1.vel, 
                      subset=nFeature_spliced > 1000 & nFeature_spliced <4000 & percent.mt<10)
sample2.vel <- subset(sample2.vel, 
                      subset=nFeature_spliced > 1000 & nFeature_spliced <4000 & percent.mt<10)
sample3.vel <- subset(sample3.vel, 
                      subset=nFeature_spliced > 1000 & nFeature_spliced <4000 & percent.mt<10)

sample1.vel <- process_standard(sample1.vel, assaytype='spliced', output='../results/tumor24_velocity_jackstraw.pdf', norm=F)
sample2.vel <- process_standard(sample2.vel, assaytype='spliced', output='../results/tumor21_velocity_jackstraw.pdf', norm=F)
sample3.vel <- process_standard(sample3.vel, assaytype='spliced', output='../results/tumor22_velocity_jackstraw.pdf', norm=F)

#saveRDS(list(sample1.vel, sample2.vel, sample3.vel), '../results/final_velocity_obj.RDS')

sample1.tumor.cells <- WhichCells(sample1.obj, ident=c(0, 3))
sample2.tumor.cells <- WhichCells(sample2.merge.obj, ident=c(0, 1))
sample3.tumor.cells <- WhichCells(sample3.merge.obj, ident=c(0, 1, 4))

## add prefix and suffix for veloctiy
sample1.tumor.cellsv <- paste0('Tumor24_zebrafish_with_orf_color_v2:', sample1.tumor.cells, 'x')

sample2.tumor.cellsv <- rep("", length(sample2.tumor.cells))

sample2.tumor.cellsv[grepl('sort_', sample2.tumor.cells)] = paste0('sort_Tumor21_zebrafish_with_orf_color_v2:', gsub('sort_', '', sample2.tumor.cells[grepl('sort_', sample2.tumor.cells)]), 'x')

sample2.tumor.cellsv[grepl('bulk_', sample2.tumor.cells)] = paste0('bulk_Tumor21_2ndlibrary_zebrafish_with_orf_color_v2:', gsub('bulk_', '', sample2.tumor.cells[grepl('bulk_', sample2.tumor.cells)]), 'x')

sample3.tumor.cellsv <- rep("", length(sample3.tumor.cells))

sample3.tumor.cellsv[grepl('sort_', sample3.tumor.cells)] = paste0('sort_Tumor22_zebrafish_with_orf_color_v2:', gsub('sort_', '', sample3.tumor.cells[grepl('sort_', sample3.tumor.cells)]), 'x')

sample3.tumor.cellsv[grepl('bulk_', sample3.tumor.cells)] = paste0('bulk_Tumor22_2ndlibrary_zebrafish_with_orf_color_v2:', gsub('bulk_', '', sample3.tumor.cells[grepl('bulk_', sample3.tumor.cells)]), 'x')

sample1.vel <- subset(sample1.vel, cells=sample1.tumor.cellsv)
sample2.vel  <- subset(sample2.vel, cells=sample2.tumor.cellsv)
sample3.vel  <- subset(sample3.vel, cells=sample3.tumor.cellsv)

sample1.tumor.cells <- sample1.tumor.cells[sample1.tumor.cellsv%in%rownames(sample1.vel@meta.data)]
sample2.tumor.cells <- sample2.tumor.cells[sample2.tumor.cellsv%in%rownames(sample2.vel@meta.data)]
sample3.tumor.cells <- sample3.tumor.cells[sample3.tumor.cellsv%in%rownames(sample3.vel@meta.data)]

sample1.obj <- subset(sample1.obj, cells=sample1.tumor.cells)
sample2.obj <- subset(sample2.merge.obj, cells=sample2.tumor.cells)
sample3.obj <- subset(sample3.merge.obj, cells=sample3.tumor.cells)

# avoid intersection of rna transcripts to get back the color sequences
# sample1.rna <- intersect(rownames(sample1.obj), rownames(sample1.vel))
# sample2.rna <- intersect(rownames(sample2.obj), rownames(sample2.vel))
# sample3.rna <- intersect(rownames(sample3.obj), rownames(sample3.vel))
# sample1.obj <- subset(sample1.obj, features=sample1.rna)
# sample1.vel <- subset(sample1.vel, features=sample1.rna)
# sample2.obj <- subset(sample2.obj, features=sample2.rna)
# sample2.vel <- subset(sample2.vel, features=sample2.rna)
# sample3.obj <- subset(sample3.obj, features=sample3.rna)
# sample3.vel <- subset(sample3.vel, features=sample3.rna)

sample1.vel <- recluster.withtree(sample1.vel, name='Tumor24')
sample1.obj <- recluster.withtree(sample1.obj, name='Tumor24')
sample2.vel <- recluster.withtree(sample2.vel, name='Tumor21')
sample2.obj <- recluster.withtree(sample2.obj, name='Tumor21')
sample3.vel <- recluster.withtree(sample3.vel, name='Tumor22')
sample3.obj <- recluster.withtree(sample3.obj, name='Tumor22')

pdf(paste0('../results/recluster_', 'Tumor24', '.pdf'), width=8, height=20)
clustree(sample1.obj, prefix='SCT_snn_res.') # node_colour='sc3_stability')
clustree(sample1.obj, prefix='SCT_snn_res.', node_colour='gfp-mandeline-orf-chr', node_colour_aggr='mean')
clustree(sample1.obj, prefix='SCT_snn_res.', node_colour='mCherry', node_colour_aggr='mean')
clustree(sample1.obj, prefix='SCT_snn_res.', node_colour='vangl2', node_colour_aggr='mean')
clustree(sample1.obj, prefix='SCT_snn_res.', node_colour='krasg12d-madeline', node_colour_aggr='mean')
dev.off()

pdf(paste0('../results/recluster_', 'Tumor21', '.pdf'), width=8, height=20)
clustree(sample2.obj, prefix='SCT_snn_res.') # node_colour='sc3_stability')
clustree(sample2.obj, prefix='SCT_snn_res.', node_colour='gfp-mandeline-orf-chr', node_colour_aggr='mean')
clustree(sample2.obj, prefix='SCT_snn_res.', node_colour='mCherry', node_colour_aggr='mean')
clustree(sample2.obj, prefix='SCT_snn_res.', node_colour='vangl2', node_colour_aggr='mean')
clustree(sample2.obj, prefix='SCT_snn_res.', node_colour='krasg12d-madeline', node_colour_aggr='mean')
dev.off()

pdf(paste0('../results/recluster_', 'Tumor22', '.pdf'), width=8, height=20)
clustree(sample3.obj, prefix='SCT_snn_res.') #node_colour='sc3_stability')
clustree(sample3.obj, prefix='SCT_snn_res.', node_colour='gfp-mandeline-orf-chr', node_colour_aggr='mean')
clustree(sample3.obj, prefix='SCT_snn_res.', node_colour='mCherry', node_colour_aggr='mean')
clustree(sample3.obj, prefix='SCT_snn_res.', node_colour='vangl2', node_colour_aggr='mean')
clustree(sample3.obj, prefix='SCT_snn_res.', node_colour='krasg12d-madeline', node_colour_aggr='mean')
dev.off()

fish = read.delim('fish_color_table.txt', sep='\t')

sample1.obj <- FindClusters(object=sample1.obj, resolution=0.8)
sample2.obj <- FindClusters(object=sample2.obj, resolution=0.8)
sample3.obj <- FindClusters(object=sample3.obj, resolution=0.8)

saveRDS(sample1.obj, '../results/seurat/Tumor24_orig_seurat.obj')
saveRDS(sample2.obj, '../results/seurat/Tumor21_orig_seurat.obj')
saveRDS(sample3.obj, '../results/seurat/Tumor22_orig_seurat.obj')

## double check fish cell labels
pdf('test_fish.pdf', width=30, height=8)
new.ids = paste0(seq(0, 9), '_', fish[, 'Tumor24'])
names(new.ids) = levels(sample1.obj)
sample1.obj <- RenameIdents(sample1.obj, new.ids)
p1=DimPlot(sample1.obj, reduction='umap', label=T, pt.size=0.5)
new.ids = paste0(seq(0, 8), '_', fish[, 'Tumor21'][fish[, 'Tumor21']!=''])
names(new.ids) = levels(sample2.obj)
sample2.obj <- RenameIdents(sample2.obj, new.ids)
p2=DimPlot(sample2.obj, reduction='umap', label=T, pt.size=0.5)
new.ids = paste0(seq(0, 8), '_', fish[, 'Tumor22'][fish[, 'Tumor22']!=''])
names(new.ids) = levels(sample3.obj)
sample3.obj <- RenameIdents(sample3.obj, new.ids)
p3=DimPlot(sample3.obj, reduction='umap', label=T, pt.size=0.5)
CombinePlots(plots=list(p1,p2,p3), ncol=3)
dev.off()

 ## when using scale.data, the fold change header is named as avg_diff
sample1.obj.markers = FindAllMarkers(sample1.obj, only.pos=T, min.pct=0.1, test.use='MAST', assay='SCT', slot='scale.data',
    				     random.seed=100, logfc.threshold = 0.1)
sample2.obj.markers = FindAllMarkers(sample2.obj, only.pos=T, min.pct=0.1, test.use='MAST', assay='SCT', slot='scale.data',
    				     random.seed=100, logfc.threshold = 0.1)
sample3.obj.markers = FindAllMarkers(sample3.obj, only.pos=T, min.pct=0.1, test.use='MAST', assay='SCT', slot='scale.data',
    				     random.seed=100, logfc.threshold = 0.1)

human_ortholog = read.table('~/alvin_singlecell/01_rms_projects//01_fish/data/ortholog_mapping/Beagle_fish_human_all_genes.txt', header=T, sep='\t', stringsAsFactors=F)
sample1.results <- cbind(sample1.obj.markers, human_ortholog[match(sample1.obj.markers$gene, human_ortholog$Gene), ])
sample1.results$diff.enrich  <- sample1.results$pct.1 / sample1.results$pct.2
sample2.results <- cbind(sample2.obj.markers, human_ortholog[match(sample2.obj.markers$gene, human_ortholog$Gene), ])
sample2.results$diff.enrich  <- sample2.results$pct.1 / sample2.results$pct.2
sample3.results <- cbind(sample3.obj.markers, human_ortholog[match(sample3.obj.markers$gene, human_ortholog$Gene), ])
sample3.results$diff.enrich  <- sample3.results$pct.1 / sample3.results$pct.2

## default using slot data, the fold change header is named as avg_logFC for log fold change
## sample1.obj.markers2 = FindAllMarkers(sample1.obj, only.pos=T, min.pct=0.1, assay='SCT', random.seed=100, logfc.threshold = 0.2)
## sample1.obj.results2 <- scDE.output(sample1.obj, sample1.obj.markers2) ## double check avg_logFC, it's kind of opposite......, it's different from sowmya's codes
write.table(sample1.results, sep='\t', quote=F, file='Tumor24_tumoronly_res0.8.xls')
write.table(sample2.results, sep='\t', quote=F, file='Tumor21_tumoronly_res0.8.xls')
write.table(sample3.results, sep='\t', quote=F, file='Tumor22_tumoronly_res0.8.xls')
 
markers <- c('myf5', 'krasg12d-madeline', 'mCherry', 'vangl2', 'gfp-mandeline-orf-chr')
allmarkers <- readRDS('../results/final_manual_markers.RDS')
for (category in names(allmarkers$rms_human_markers)) {
    pdf(paste0(paste0('../results/human_rms_', category, '_markers_on_zebrafish.pdf')), width=32, height=5)
    markers = na.omit(allmarkers$rms_human_markers[[category]][,1])
    markers.copy = intersect(markers, rownames(sample1.obj$SCT))
    p1 = DotPlot(sample1.obj, features = markers.copy) + RotatedAxis()
    p1h = DoHeatmap(subset(sample1.obj, downsample = 100), features=markers.copy, size = 3)
    markers.copy = intersect(markers, rownames(sample2.obj$SCT))
    p2 = DotPlot(sample2.obj, features = markers.copy) + RotatedAxis()
    p2h = DoHeatmap(subset(sample2.obj, downsample = 100), features=markers.copy, size = 3)
    markers.copy = intersect(markers, rownames(sample3.obj$SCT))
    p3 = DotPlot(sample3.obj, features = markers.copy) + RotatedAxis()
    p3h = DoHeatmap(subset(sample3.obj, downsample = 100), features=markers.copy, size = 3)
    print(CombinePlots(plots=list(p1, p2, p3), ncol=3))
    print(CombinePlots(plots=list(p1h, p2h, p3h), ncol=3))
    dev.off()
}


dave.zebrafish.markers11 <- lapply(Sys.glob('../data/markers/Zebrafish_Signatures/*'), scan, what='')
names(dave.zebrafish.markers11) <- gsub('\\.txt', '', basename(Sys.glob('../data/markers/Zebrafish_Signatures/*')))

for (category in 1:length(dave.zebrafish.markers11)) {
    pdf(paste0(paste0('../results/human_rms_', names(dave.zebrafish.markers11)[category], '_gsea_dave_markers_on_zebrafish.pdf')), width=32, height=5)
    markers = na.omit(dave.zebrafish.markers11[[category]])
    markers <- human_ortholog[match(markers, human_ortholog$Hsortholog), 2]
    markers.copy = intersect(markers, rownames(sample1.obj$SCT))
    if (names(dave.zebrafish.markers11)[category] == 'Hayes_List') {
         markers.copy <- c('krasg12d-madeline', 'gfp-mandeline-orf-chr', 'mCherry', markers.copy)
    }
    p1 = DotPlot(sample1.obj, features = markers.copy) + RotatedAxis()
    p1h = DoHeatmap(subset(sample1.obj, downsample = 100), features=markers.copy, size = 3)
    markers.copy = intersect(markers, rownames(sample2.obj$SCT))
    if (names(dave.zebrafish.markers11)[category] == 'Hayes_List') {
         markers.copy <- c('krasg12d-madeline', 'gfp-mandeline-orf-chr', 'mCherry', markers.copy)
    }
    p2 = DotPlot(sample2.obj, features = markers.copy) + RotatedAxis()
    p2h = DoHeatmap(subset(sample2.obj, downsample = 100), features=markers.copy, size = 3)
    markers.copy = intersect(markers, rownames(sample3.obj$SCT))
    if (names(dave.zebrafish.markers11)[category] == 'Hayes_List') {
         markers.copy <- c('krasg12d-madeline', 'gfp-mandeline-orf-chr', 'mCherry', markers.copy)
    }
    p3 = DotPlot(sample3.obj, features = markers.copy) + RotatedAxis()
    p3h = DoHeatmap(subset(sample3.obj, downsample = 100), features=markers.copy, size = 3)
    print(CombinePlots(plots=list(p1, p2, p3), ncol=3))
    print(CombinePlots(plots=list(p1h, p2h, p3h), ncol=3))
    dev.off()
#    break
}

#for (category in 4:length(dave.zebrafish.markers11)) {
for (category in 7:7) {
    markers = na.omit(dave.zebrafish.markers11[[category]])
    markers <- human_ortholog[match(markers, human_ortholog$Hsortholog), 2]
    markers.copy = intersect(markers, rownames(sample1.vel$SCT))
    markers.copy = intersect(markers.copy, rownames(sample2.vel$SCT))
    markers.copy = intersect(markers.copy, rownames(sample3.vel$SCT))
    markers.copy <- c('krasg12d-madeline', 'gfp-mandeline-orf-chr', 'mCherry', markers.copy)
    pdf(paste0('../results/recluster_tsne_umap_', names(dave.zebrafish.markers11)[category], '_gsea_dave_markers_threefish_tumors_res0.8', '.pdf'), width=18, height=8)
    for (m in markers.copy) {
        p1 <- DimPlot(sample1.obj, reduction='umap')
        p2 <- DimPlot(sample2.obj, reduction='umap')
        p3 <- DimPlot(sample3.obj, reduction='umap')
        p1.1 <- FeaturePlot(sample1.obj, reduction='umap', features=m)
        p2.1 <- FeaturePlot(sample2.obj, reduction='umap', features=m)
        p3.1 <- FeaturePlot(sample3.obj, reduction='umap', features=m)
        #p1 <- DimPlot(sample1.vel, reduction='umap')
        #p2 <- DimPlot(sample2.vel, reduction='umap')
        #p3 <- DimPlot(sample3.vel, reduction='umap')
        #p1.1 <- FeaturePlot(sample1.vel, reduction='umap', features=m)
        #p2.1 <- FeaturePlot(sample2.vel, reduction='umap', features=m)
        #p3.1 <- FeaturePlot(sample3.vel, reduction='umap', features=m)
        print(CombinePlots(plots=list(CombinePlots(plots=list(p1, p2, p3), ncol=3), CombinePlots(plots=list(p1.1, p2.1, p3.1), ncol=3)), ncol=1))
    }
    dev.off()
#    break
}

library(gridExtra)
library(grid)
library(ggpubr)
set.seed(100)
pdf(paste0('../results/recluster_tsne_umap', 'Tumor24_res0.8', '.pdf'), width=40, height=10)
p1 <- DimPlot(sample1.obj, reduction='tsne', group.by=c('orig.ident', 'seurat_clusters'))
p2 <- DimPlot(sample1.obj, reduction='umap', group.by=c('orig.ident', 'seurat_clusters'))
p3 <- FeaturePlot(sample1.obj, reduction='tsne', features=markers, ncol=5)
p4 <- FeaturePlot(sample1.obj, reduction='umap', features=markers, ncol=5)
#CombinePlots(plots=list(CombinePlots(plots=list(p1, p2), ncol=1), CombinePlots(plots=list(p3, p4), ncol=1)), ncol=2)
plots.all <- ggarrange(CombinePlots(plots=list(p1, p2), ncol=1), CombinePlots(plots=list(p3, p4), ncol=1), widths = c(2, 5), ncol=2)
plots.all
dev.off()

pdf(paste0('../results/recluster_tsne_umap', 'Tumor21_res0.8', '.pdf'), width=40, height=10)
p1 <- DimPlot(sample2.obj, reduction='tsne', group.by=c('orig.ident', 'seurat_clusters'))
p2 <- DimPlot(sample2.obj, reduction='umap', group.by=c('orig.ident', 'seurat_clusters'))
p3 <- FeaturePlot(sample2.obj,        reduction='tsne', features=markers, ncol=5)
p4 <- FeaturePlot(sample2.obj,        reduction='umap', features=markers, ncol=5)
plots.all <- ggarrange(CombinePlots(plots=list(p1, p2), ncol=1), CombinePlots(plots=list(p3, p4), ncol=1), widths = c(2, 5), ncol=2)
plots.all
dev.off()

pdf(paste0('../results/recluster_tsne_umap', 'Tumor22_res0.8', '.pdf'), width=40, height=10)
p1 <- DimPlot(sample3.obj, reduction='tsne', group.by=c('orig.ident', 'seurat_clusters'))
p2 <- DimPlot(sample3.obj, reduction='umap', group.by=c('orig.ident', 'seurat_clusters'))
p3 <- FeaturePlot(sample3.obj,        reduction='tsne', features=markers, ncol=5)
p4 <- FeaturePlot(sample3.obj,        reduction='umap', features=markers, ncol=5)
plots.all <- ggarrange(CombinePlots(plots=list(p1, p2), ncol=1), CombinePlots(plots=list(p3, p4), ncol=1), widths = c(2, 5), ncol=2)
plots.all
dev.off()

# Velocity analysis
sample1.vel <- RunVelocity(object=sample1.vel, deltaT = 1, kCells = 25, fit.quantile = 0.02)
sample2.vel <- RunVelocity(object=sample2.vel, deltaT = 1, kCells = 25, fit.quantile = 0.02)
sample3.vel <- RunVelocity(object=sample3.vel, deltaT = 1, kCells = 25, fit.quantile = 0.02)

# label cluster by seurat processed rna results
## a bit different in terms of the velocity spliced rna clusterings
#sample1.vel <- FindClusters(object=sample1.vel, resolution=0.8)
#sample2.vel <- FindClusters(object=sample2.vel, resolution=0.8)
#sample3.vel <- FindClusters(object=sample3.vel, resolution=0.8)

## label cluster by seurat processed rna results
sample1.vel.copy <- sample1.vel
sample1.vel.copy$seurat_clusters <- sample1.obj$seurat_clusters

sample2.vel.copy <- sample2.vel
sample2.vel.copy$seurat_clusters <- sample2.obj$seurat_clusters
sample3.vel.copy <- sample3.vel
sample3.vel.copy$seurat_clusters <- sample3.obj$seurat_clusters

pdf(paste0('../results/recluster_tsne_umap_velo', 'Tumor24_res0.8', '.pdf'), width=13, height=10)
p1 <- DimPlot(sample1.vel.copy, reduction='tsne', group.by=c('orig.ident', 'seurat_clusters'))
p2 <- DimPlot(sample1.vel.copy, reduction='umap', group.by=c('orig.ident', 'seurat_clusters'))
CombinePlots(plots=list(p1, p2), ncol=1)
dev.off()

pdf(paste0('../results/recluster_tsne_umap_velo', 'Tumor21_res0.8', '.pdf'), width=13, height=10)
p1 <- DimPlot(sample2.vel.copy, reduction='tsne', group.by=c('orig.ident', 'seurat_clusters'))
p2 <- DimPlot(sample2.vel.copy, reduction='umap', group.by=c('orig.ident', 'seurat_clusters'))
CombinePlots(plots=list(p1, p2), ncol=1)
dev.off()
pdf(paste0('../results/recluster_tsne_umap_velo', 'Tumor22_res0.8', '.pdf'), width=13, height=10)
p1 <- DimPlot(sample3.vel.copy, reduction='tsne', group.by=c('orig.ident', 'seurat_clusters'))
p2 <- DimPlot(sample3.vel.copy, reduction='umap', group.by=c('orig.ident', 'seurat_clusters'))
CombinePlots(plots=list(p1, p2), ncol=1)
dev.off()

ident.colors <- (scales::hue_pal())(length(levels(sample1.obj)))
names(ident.colors) <- levels(sample1.obj)
cell.colors <- ident.colors[Idents(sample1.obj)]
names(cell.colors) <- colnames(sample1.vel.copy)
## pdf('../results/tumor24_velocity_tumoronly.pdf', width=9, height=9)
pdf('../results/tumor24_velocity_tumoronly_redo.pdf', width=9, height=9)
v1 = show.velocity.on.embedding.cor(emb = Embeddings(object=sample1.vel.copy, reduction = "tsne"), 
			            vel = Tool(object=sample1.vel.copy, slot = "RunVelocity"), 
			            n = 200, scale = "sqrt", 
			            cell.colors = ac(cell.colors, alpha = 0.5), return.details=T,
			            cex=1, 
			            arrow.scale = 3.5, 
			            show.grid.flow = TRUE, 
			            min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
legend('topright', legend=levels(sample1.obj), pch=19, col=ident.colors)
## v2 = show.velocity.on.embedding.cor(emb = Embeddings(object=sample1.vel.copy, reduction = "umap"), 
## 			       vel = Tool(object=sample1.vel.copy, slot = "RunVelocity"), 
## 			       n = 200, scale = "sqrt", 
## 			       cell.colors = ac(cell.colors, alpha = 0.5), return.details=T,
## 			       cex=1, 
## 			       arrow.scale = 1.5, 
## 			       show.grid.flow = TRUE, 
## 			       min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
## legend('topright', legend=levels(sample1.obj), pch=19, col=ident.colors)
dev.off()

ident.colors <- (scales::hue_pal())(length(levels(sample2.obj)))
names(ident.colors) <- levels(sample2.obj)
cell.colors <- ident.colors[Idents(sample2.obj)]
names(cell.colors) <- colnames(sample2.vel.copy)
pdf('../results/tumor21_velocity_tumoronly.pdf', width=9, height=9)
v1 = show.velocity.on.embedding.cor(emb = Embeddings(object=sample2.vel.copy, reduction = "tsne"), 
			       vel = Tool(object=sample2.vel.copy, slot = "RunVelocity"), 
			       n = 200, scale = "sqrt", 
			       cell.colors = ac(cell.colors, alpha = 0.5), return.details=T,
			       cex=1, 
			       arrow.scale = 2, 
			       show.grid.flow = TRUE, 
			       min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
legend('topright', legend=levels(sample2.obj), pch=19, col=ident.colors)
v2= show.velocity.on.embedding.cor(emb = Embeddings(object=sample2.vel.copy, reduction = "umap"), 
			       vel = Tool(object=sample2.vel.copy, slot = "RunVelocity"), 
			       n = 200, scale = "sqrt", 
			       cell.colors = ac(cell.colors, alpha = 0.5),
			       cex=1, 
			       arrow.scale = 2, 
			       show.grid.flow = TRUE, 
			       min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
legend('topright', legend=levels(sample2.obj), pch=19, col=ident.colors)
dev.off()

ident.colors <- (scales::hue_pal())(length(levels(sample3.obj)))
names(ident.colors) <- levels(sample3.obj)
cell.colors <- ident.colors[Idents(sample3.obj)]
names(cell.colors) <- colnames(sample3.vel.copy)
pdf('../results/tumor22_velocity_tumoronly.pdf', width=9, height=9)
v1 = show.velocity.on.embedding.cor(emb = Embeddings(object=sample3.vel.copy, reduction = "tsne"), 
			       vel = Tool(object=sample3.vel.copy, slot = "RunVelocity"), 
			       n = 200, scale = "sqrt", 
			       cell.colors = ac(cell.colors, alpha = 0.5), return.details=T,
			       cex=1, 
			       arrow.scale = 3, 
			       show.grid.flow = TRUE, 
			       min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
legend('topright', legend=levels(sample3.obj), pch=19, col=ident.colors)
v2 = show.velocity.on.embedding.cor(emb = Embeddings(object=sample3.vel.copy, reduction = "umap"), 
			       vel = Tool(object=sample3.vel.copy, slot = "RunVelocity"), 
			       n = 300, scale = "sqrt", 
			       cell.colors = ac(cell.colors, alpha = 0.5),
			       cex=1, 
			       arrow.scale = 3, 
			       show.grid.flow = TRUE, 
			       min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = FALSE, cell.border.alpha = 0.1)
legend('topright', legend=levels(sample3.obj), pch=19, col=ident.colors)
dev.off()

## Run monocle3 on zebrafish datasets
## recluster for tumor 22 cluster 6 and cluster 7
sample3.obj <- FindClusters(sample3.obj, resolution=0.9)

pdf('../results/Tumor22_cluster6_recluster.pdf', width=16, height=7)
p1 <- DimPlot(sample3.obj, reduction='umap')
p2 <- DimPlot(sample3.obj, reduction='tsne')
CombinePlots(plots=list(p1, p2), ncol=2)
dev.off()

saveRDS(sample1.obj@assays[['SCT']]@scale.data, file='../results/tumor24_sct_scale_mat.rds')

## Combined cluster for DE analysis 
## https://bioinformatics.stackexchange.com/questions/4249/manually-define-clusters-in-seurat-and-determine-marker-genes
## Manual determination of the cluster number
sample1.obj <- FindClusters(object=sample1.obj, resolution=0.8)

levels(sample1.obj@active.ident)[levels(sample1.obj@active.ident) %in% c(0, 1, 2, 6)] = 10
levels(sample1.obj$seurat_clusters)[levels(sample1.obj$seurat_clusters) %in% c(0, 1, 2, 6)] = 10

sample1.obj.mergecluster.markers = FindMarkers(sample1.obj, only.pos=T, min.pct=0.1, test.use='MAST', assay='SCT', slot='scale.data',
					       ident.1 = 10,
					       random.seed=100)
