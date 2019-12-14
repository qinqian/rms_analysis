#rmote::start_rmote()
library(Seurat)
library(gdata)
library(tidyverse)
library(Matrix)
library(clustree)
library(pheatmap)
library(foreach)
library(SeuratWrappers)
library(garnett)
library(fgsea)

#library(cellassign)
library(UpSetR)
set.seed(100)

source("functions.R")
human_ortholog = read.table('~/alvin_singlecell/01_rms_projects//01_fish/data/ortholog_mapping/Beagle_fish_human_all_genes.txt', header=T, sep='\t', stringsAsFactors=F)
dave_top3_marker_names = Sys.glob("/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/markers/top3_markers/*")
dave_top3_marker <- lapply(dave_top3_marker_names, scan, what="")
names(dave_top3_marker) <- gsub('.txt', '', basename(dave_top3_marker_names))
sara_markers = read.xls('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/markers/raw/non_tumor_cell_markers.xlsx')
sara_markers_list  = list()
for (s in colnames(sara_markers)[1:6]) {
    gene_set = as.vector(sara_markers[[s]])
    sara_markers_list[[s]] = gene_set[gene_set!= ""]
}
individual_markers = c('vangl2', 'myf5', 'mCherry', 'GFP', 'myc_mus_musculus',
                       'krasg12d_madeline', 'dtomato_yan', 'kaede_yan', 'gfp_mandeline_orf_chr',
		       'KRASG12D_mandeline_second_orf_chr',
                       'cas9_ally', 'Cyan', 'dsRed', 'ZsYellow1')
sara_markers_list$tumor_markers = individual_markers

indrop_markers_name = Sys.glob("../data/markers/dave_jme_indrop/*")
indrop_markers = lapply(indrop_markers_name, function(x) {
    y = read.delim(x, na.strings="")
    y[, c(1, ncol(y))]
})
names(indrop_markers) = gsub('.txt', '', basename(indrop_markers_name))
smartseq_markers_name = Sys.glob("../data/markers/qin_jme_smart_seq//*")
smartseq_markers = lapply(smartseq_markers_name, function(x) {
    y = read.delim(x, na.strings="")
    y[, c(1, ncol(y))]
})
names(smartseq_markers) = gsub('.txt', '', basename(smartseq_markers_name))
dave_human_rms_markers_name = Sys.glob("../data/markers/dave_annotated_human_rms_genes/*")
dave_human_rms_markers = lapply(dave_human_rms_markers_name, function(x) {
    y=scan(x, what="")
    cbind(human_ortholog[match(y, human_ortholog[,3]), 2], y)
})
names(dave_human_rms_markers) = gsub('.txt', '', basename(dave_human_rms_markers_name))

all_markers_final = list(tumor_related=sara_markers_list, top3_markers=dave_top3_marker, rms_human_markers=dave_human_rms_markers,
          		 smart_seq=smartseq_markers, in_drop=indrop_markers)

all_markers_final$tumor_related$Macrophages <- c(all_markers_final$tumor_related$Macrophages, 'fcer1g')
all_markers_final$tumor_related$Macrophages <- c(all_markers_final$tumor_related$Macrophages, 'tyrobp')
all_markers_final$tumor_related$Macrophages <- c(all_markers_final$tumor_related$Macrophages, 'csf1r')
all_markers_final$tumor_related$oligodendrocytes <- c('mbpa', 'tfa', 'plp1a', 'mag', 'cldn11b')


#data_dir <- '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/'
#sample1 <- Read10X(paste0(data_dir, 'Tumor24_zebrafish_with_orf_color/outs/filtered_feature_bc_matrix'))
#sample2 <- Read10X(paste0(data_dir, 'Tumor21_zebrafish_with_orf_color/outs/filtered_feature_bc_matrix'))
#sample2_2 <- Read10X(paste0(data_dir, 'Tumor21_2ndlibrary_zebrafish_with_orf_color/outs/filtered_feature_bc_matrix'))
#sample3 <- Read10X(paste0(data_dir, 'Tumor22_zebrafish_with_orf_color/outs/filtered_feature_bc_matrix'))
#sample3_2 <- Read10X(paste0(data_dir, 'Tumor22_2ndlibrary_zebrafish_with_orf_color/outs/filtered_feature_bc_matrix'))

sample1 <- Read10X('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix/')
sample2 <- Read10X('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix/')
sample2_2 <- Read10X('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix/')
sample3 <- Read10X('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix/')
sample3_2 <- Read10X('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix/')
library_names <- c("Tumor24", "Tumor21_sort", "Tumor21_bulk", "Tumor22_sort", "Tumor22_bulk")

sample1.obj <- CreateSeuratObject(counts=sample1, project='RMS24', min.cells=3, min.features=10)
sample2.obj <- CreateSeuratObject(counts=sample2, project='RMS21sort', min.cells=3, min.features=10)
sample2_2.obj <- CreateSeuratObject(counts=sample2_2, project='RMS21bulk', min.cells=3, min.features=10)
sample3.obj <- CreateSeuratObject(counts=sample3, project='RMS22sort', min.cells=3, min.features=10)
sample3_2.obj <- CreateSeuratObject(counts=sample3_2, project='RMS22bulk', min.cells=3, min.features=10)

sample2.merge.obj <- merge(sample2.obj, sample2_2.obj, add.cell.ids = c("sort", 'bulk'), project = 'RMS21_zebrafish')
sample3.merge.obj <- merge(sample3.obj, sample3_2.obj, add.cell.ids = c("sort", 'bulk'), project = 'RMS22_zebrafish')

samples.list <- Reduce(rbind, lapply(list(sample1, sample2, sample2_2, sample3, sample3_2), function(x) {
  rowSums(x[individual_markers, ])	        
}))
rownames(samples.list) <- library_names
samples.list <- samples.list[, colSums(samples.list) != 0]
expressed.markers <- colnames(samples.list)
pheatmap(samples.list, display_numbers=T, cex=1.2, filename='../results/3fish_QC_markers.pdf', height=10, width=10)
sample1.obj[["percent.mt"]] <- PercentageFeatureSet(sample1.obj, pattern='^mt-')
sample2.merge.obj[["percent.mt"]] <- PercentageFeatureSet(sample2.merge.obj, pattern='^mt-')
sample3.merge.obj[["percent.mt"]] <- PercentageFeatureSet(sample3.merge.obj, pattern='^mt-')
#plot1 = VlnPlot(sample1.obj, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
#plot2 = VlnPlot(sample2.merge.obj, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
#plot3 = VlnPlot(sample3.merge.obj, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
#options(repr.plot.width=8, repr.plot.height=18)
#CombinePlots(plots=list(plot1, plot2, plot3), ncol=1)
#ggsave("../results/Cell_QC_violin_three_samples.pdf", width=8, height=18)

pdf("../results/Cell_QC_histgram_three_samples.pdf", width=18, height=12)
par(mfrow=c(2, 3), font=2, cex=1.5, mar=c(5, 5, 2, 0))
hist(sample2.merge.obj@meta.data$nFeature_RNA, n=50, xlab='Number of detected gene number', ylab="Cell barcode number",
     main='QC of tumor 21 ', col='red')
hist(sample3.merge.obj@meta.data$nFeature_RNA, xlab='Number of detected gene number', ylab="Cell barcode number",
     main='QC of tumor 22 ', n=50, col="blue")
hist(sample1.obj@meta.data$nFeature_RNA, xlab='Number of detected gene number', ylab="Cell barcode number",
     main='QC of tumor 24', n=50, col="blue")
hist(sample2.merge.obj@meta.data$percent.mt, n=50, xlab="Mitochondria ratio", ylab="Cell barcode number",
     main='QC of tumor 21 ', col='red')
hist(sample3.merge.obj@meta.data$percent.mt, xlab="Mitochondria ratio", ylab="Cell barcode number",
     main='QC of tumor 22 ', n=50, col="blue")
hist(sample1.obj@meta.data$percent.mt, xlab="Mitochondria ratio", ylab="Cell barcode number",
     main='QC of tumor 24', n=50, col="blue")
dev.off()

qc1 <- read.csv('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color_v2/outs/metrics_summary.csv')
qc2 <- read.csv('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish_with_orf_color_v2/outs/metrics_summary.csv')
qc2_2 <- read.csv('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2/outs/metrics_summary.csv')
qc3 <- read.csv('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish_with_orf_color_v2/outs/metrics_summary.csv')
qc3_2 <- read.csv('/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2/outs/metrics_summary.csv')
qc.default <- rbind(qc1, qc2, qc2_2, qc3, qc3_2)
rownames(qc.default) <- c("Tumor24", "Tumor21 sort", "Tumor21 bulk", "Tumor22 sort", "Tumor22 bulk")
write.table(qc.default, file="../results/cellranger_qc_report.xls", sep="\t", quote=F)

sample1.obj <- subset(sample1.obj, 
                      subset=nFeature_RNA>1000 & nFeature_RNA<4000 & percent.mt<10)
sample2.merge.obj <- subset(sample2.merge.obj, 
                      subset=nFeature_RNA>1000 & nFeature_RNA<4000 & percent.mt<10)
sample3.merge.obj <- subset(sample3.merge.obj, 
                      subset=nFeature_RNA>1000 & nFeature_RNA<4000 & percent.mt<10)

sample1.obj = process_standard(sample1.obj, norm=F, output='../results/sample1_pc_jackstraw.pdf')
pdf("../results/sample1_decide_resolution.pdf", width=6, height=10)
clustree(sample1.obj, prefix='SCT_snn_res.')
dev.off()

sample2.merge.obj = process_standard(sample2.merge.obj, norm=F, output='../results/sample2_pc_jackstraw.pdf')
pdf("../results/sample2_decide_resolution.pdf", width=6, height=10)
clustree(sample2.merge.obj, prefix='SCT_snn_res.')
dev.off()

sample3.merge.obj = process_standard(sample3.merge.obj, norm=F, output='../results/sample3_pc_jackstraw.pdf')
pdf("../results/sample3_decide_resolution.pdf", width=6, height=10)
clustree(sample3.merge.obj, prefix='SCT_snn_res.')
dev.off()

sample1.obj  <- FindClusters(sample1.obj, resolution=0.05)
sample2.merge.obj  <- FindClusters(sample2.merge.obj, resolution=0.1)
sample3.merge.obj  <- FindClusters(sample3.merge.obj, resolution=0.1)

saveRDS(sample1.obj@assays[['RNA']]@data, file='../results/tumor24_sct_mat_withblood.rds')
saveRDS(as.data.frame(sample1.obj$seurat_clusters), file='../results/tumor24_sct_cluster_withblood.rds')
saveRDS(sample2.merge.obj@assays[['RNA']]@data, file='../results/tumor21_sct_mat_withblood.rds')
saveRDS(as.data.frame(sample2.merge.obj$seurat_clusters), file='../results/tumor21_sct_cluster_withblood.rds')
saveRDS(sample3.merge.obj@assays[['RNA']]@data, file='../results/tumor22_sct_mat_withblood.rds')
saveRDS(as.data.frame(sample3.merge.obj$seurat_clusters), file='../results/tumor22_sct_cluster_withblood.rds')

pdf("../results/tumor_zebrafish_umap_vs_tsne.pdf", width=30, height=8)
p1 = DimPlot(object=sample1.obj, reduction='umap', group.by=c("orig.ident", "seurat_clusters"), ncol=1)
p2 = DimPlot(object=sample1.obj, reduction='tsne', group.by=c("orig.ident", "seurat_clusters"), ncol=1)
p3 = DimPlot(object=sample2.merge.obj, reduction='umap', group.by=c("orig.ident", "seurat_clusters"), ncol=1)
p4 = DimPlot(object=sample2.merge.obj, reduction='tsne', group.by=c("orig.ident", "seurat_clusters"), ncol=1)
p5 = DimPlot(object=sample3.merge.obj, reduction='umap', group.by=c("orig.ident", "seurat_clusters"), ncol=1)
p6 = DimPlot(object=sample3.merge.obj, reduction='tsne', group.by=c("orig.ident", "seurat_clusters"), ncol=1)
CombinePlots(plots=list(p1, p2, p3, p4, p5, p6), ncol=6)
dev.off()

sample1.markers = FindAllMarkers(sample1.obj, only.pos=T, min.pct=0.2, test.use='MAST', assay='SCT', slot='scale.data', # use scale data instead of data
                                 random.seed=100, logfc.threshold = 0.2)
sample2.markers = FindAllMarkers(sample2.merge.obj, only.pos=T, min.pct=0.2, test.use='MAST', assay='SCT', slot='scale.data',
                                 random.seed=100, logfc.threshold = 0.2)
sample3.markers = FindAllMarkers(sample3.merge.obj, only.pos=T, min.pct=0.2, test.use='MAST', assay='SCT', slot='scale.data',
                                 random.seed=100, logfc.threshold = 0.2)

plots <- CombinePlots(plots=list(p1, p2, p3, p4, p5, p6), ncol=6)


#saveRDS(list(sample1.obj, sample2.merge.obj, sample3.merge.obj), '../results/final_seurat_normalize_cluster_obj.RDS')
#saveRDS(list(sample1.markers, sample2.markers, sample3.markers), '../results/final_seurat_normalize_cluster_DEgenes.RDS')
#saveRDS(all_markers_final, '../results/final_manual_markers.RDS')

samples <- readRDS('../results/final_seurat_normalize_cluster_obj.RDS')
sample1.obj <- samples[[1]]
sample2.merge.obj <- samples[[2]]
sample3.merge.obj <- samples[[3]]

samples.marker <- readRDS('../results/final_seurat_normalize_cluster_DEgenes.RDS')
sample1.marker <- samples.marker[[1]]
sample2.marker <- samples.marker[[2]]
sample3.marker <- samples.marker[[3]]

expressed.markers <- gsub('_', '-', expressed.markers)
expressed.markers  <- intersect(rownames(sample1.obj), expressed.markers)

library(ggpubr)
pdf("../results/tumor_markers_zebrafish.pdf", width=30, height=10)
for (m in expressed.markers) {
    p7 = FeaturePlot(sample1.obj,        reduction='umap', features=m)
    p8 = FeaturePlot(sample1.obj,        reduction='tsne', features=m)
    p9 = FeaturePlot(sample2.merge.obj,  reduction='umap', features=m)
    p10 = FeaturePlot(sample2.merge.obj, reduction='tsne', features=m)
    p11 = FeaturePlot(sample3.merge.obj, reduction='umap', features=m)
    p12 = FeaturePlot(sample3.merge.obj, reduction='tsne', features=m)
    plot.markers <- CombinePlots(list(p7, p8, p9, p10, p11, p12), ncol=6)
    plots.all <- ggarrange(plots, plot.markers, heights = c(2, 1), ncol=1)
    print(plots.all)
}
dev.off()

library(gridExtra)
library(grid)
pdf("../results/nontumor_markers_zebrafish.pdf", width=30, height=10)
top3markers <- unlist(all_markers_final$top3_markers)
top3markers <- top3markers[top3markers %in% rownames(sample1.obj)]
top3markers <- top3markers[top3markers %in% rownames(sample2.merge.obj)]
top3markers <- top3markers[top3markers %in% rownames(sample3.merge.obj)]
for (cell.type in 1:length(top3markers)) {
	j = top3markers[cell.type]
	cell.type = names(top3markers)[cell.type]
        p7 = FeaturePlot(sample1.obj,        reduction='umap', features=j)
        p8 = FeaturePlot(sample1.obj,        reduction='tsne', features=j)
        p9 = FeaturePlot(sample2.merge.obj,  reduction='umap', features=j)
        p10 = FeaturePlot(sample2.merge.obj, reduction='tsne', features=j)
        p11 = FeaturePlot(sample3.merge.obj, reduction='umap', features=j)
        p12 = FeaturePlot(sample3.merge.obj, reduction='tsne', features=j)
        plot.markers <- CombinePlots(list(p7, p8, p9, p10, p11, p12), ncol=6)
        grid.arrange(plots, plot.markers,
                     ncol=1,
                     top=paste0(cell.type, " ", j),
                     heights = c(2, 1),
                     bottom = textGrob(paste0(cell.type, " ", j), gp = gpar(fontface = 2, fontsize = 14),
                     hjust = 1,
                     x = 1))
}
dev.off()

all_markers_final <- readRDS('../results/final_manual_markers.RDS')
results1 <- scDE.output(sample1.obj, sample1.marker)
results2 <- scDE.output(sample2.merge.obj, sample2.marker)
results3 <- scDE.output(sample3.merge.obj, sample3.marker)
for (i in 1:length(results1)) {
    write.table(results1[[i]], file=paste0('../results/tumor24_', i-1, '_cluster_tree_res0.05.xls'), quote=F, sep='\t')
}
for (i in 1:length(results2)) {
    write.table(results2[[i]], file=paste0('../results/tumor21_', i-1, '_cluster_tree_res0.1.xls'), quote=F, sep='\t')
}
for (i in 1:length(results3)) {
    write.table(results3[[i]], file=paste0('../results/tumor22_', i-1, '_cluster_tree_res0.1.xls'), quote=F, sep='\t')
}

cluster.top3markers  <- tolower(unique(c(unlist(lapply(results1, function(x) as.vector(x$gene))), 
         unlist(lapply(results2, function(x) as.vector(x$gene))), 
         unlist(lapply(results3, function(x) as.vector(x$gene))),
	 unique(unlist(all_markers_final$top3_markers)), 
	 all_markers_final$tumor_related$tumor_markers)))
cluster.top3markers <- union(cluster.top3markers, rownames(sample1.obj[['SCT']]@data))
cluster.top3markers <- union(cluster.top3markers, rownames(sample2.merge.obj[['SCT']]@data))
cluster.top3markers <- sort(tolower(union(cluster.top3markers, rownames(sample3.merge.obj[['SCT']]@data))))
cluster.top3markers.mat <- matrix(0, nrow=length(cluster.top3markers), ncol=length(all_markers_final$top3_markers)+length(results1)+1+length(results2)+length(results3)) # add one tumor markers from FACS

rownames(cluster.top3markers.mat) <- cluster.top3markers
colnames(cluster.top3markers.mat) <- c(names(all_markers_final$top3_markers), 
				       'Tumor_FACS_marker',
                                       paste0("Tumor24_cluster", 1:length(results1)-1),
                                       paste0("Tumor21_cluster", 1:length(results2)-1),
                                       paste0("Tumor22_cluster", 1:length(results3)-1))
results.diff  <- list()
results.diff[['Tumor_FACS_marker']] <- tolower(all_markers_final$tumor_related$tumor_markers)
for (n in 1:length(results1)) {
   results.diff[[paste0('Tumor24_cluster',n-1)]] <- results1[[n]]$gene
}
for (n in 1:length(results2)) {
   results.diff[[paste0('Tumor21_cluster',n-1)]] <- results2[[n]]$gene
}
for (n in 1:length(results3)) {
   results.diff[[paste0('Tumor22_cluster',n-1)]] <- results3[[n]]$gene
}
for (n in names(all_markers_final$top3_markers)) {
    cluster.top3markers.mat[tolower(all_markers_final$top3_markers[[n]]), n] <- 1
}
for (n in names(results.diff)) {
    cluster.top3markers.mat[tolower(results.diff[[n]]), n] <- 1
}
cluster.top3markers.mat <- cluster.top3markers.mat[rowSums(cluster.top3markers.mat)>0, ]

pdf("../results/cluster3samples_consensus_24.pdf")
upset(as.data.frame(cluster.top3markers.mat), sets=c("Macrophages", "Neutrophils", names(results.diff)[1:5]), empty.intersections = "on")
dev.off()

pdf("../results/cluster3samples_consensus_21.pdf")
upset(as.data.frame(cluster.top3markers.mat), sets=c("Macrophages", "Neutrophils", names(results.diff)[6:9]), empty.intersections = "on")
dev.off()

pdf("../results/cluster3samples_consensus_22.pdf")
upset(as.data.frame(cluster.top3markers.mat), sets=c("Macrophages", "Neutrophils", names(results.diff)[10:13]), empty.intersections = "on")
dev.off()

cluster.top3markers  <- tolower(unique(c(all_markers_final$tumor_related$tumor_markers,
	 unique(unlist(all_markers_final$top3_markers)))))
cluster.top3markers = gsub('_', '-', cluster.top3markers)
cluster.top3markers <- intersect(cluster.top3markers, rownames(sample1.obj[['RNA']]@data))
cluster.top3markers.mat <- matrix(0, nrow=length(cluster.top3markers), ncol=length(all_markers_final$top3_markers)+1) # add one tumor markers from FACS
rownames(cluster.top3markers.mat) <- cluster.top3markers
colnames(cluster.top3markers.mat) <- c(names(all_markers_final$top3_markers), 'Tumor_FACS_marker')
for (n in names(all_markers_final$top3_markers)) {
    cluster.top3markers.mat[intersect(rownames(sample1.obj[['RNA']]@data), tolower(all_markers_final$top3_markers[[n]])), n] <- 1
}
cluster.top3markers.mat[intersect(rownames(sample1.obj[['RNA']]@data), gsub('_', '-', tolower(all_markers_final$tumor_related$tumor_markers))), 'Tumor_FACS_marker'] <- 1
cluster.top3markers.mat = cluster.top3markers.mat[, colSums(cluster.top3markers.mat)!=0]
library(SingleCellExperiment)
sample1.sce = subset(sample1.obj, features=rownames(cluster.top3markers.mat))
sample1.sce = as.SingleCellExperiment(sample1.sce)
sizeFactors(sample1.sce) <- colSums(assay(sample1.sce))
sample1.sce.size = sizeFactors(sample1.sce)
sample1.sce = sample1.sce[, sample1.sce.size > 0]
sample1.sce.size = sample1.sce.size[sample1.sce.size > 0]
fit <- cellassign(exprs_obj = sample1.sce,
                  marker_gene_info = cluster.top3markers.mat[rownames(sample1.sce), ],
                  s = sample1.sce.size,
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = FALSE)
print(head(celltypes(fit)))

pdf("../results/tumor24_cellassign.pdf")
pheatmap(cellprobs(fit))
dev.off()

cluster.top3markers  <- tolower(unique(c(all_markers_final$tumor_related$tumor_markers,
	 unique(unlist(all_markers_final$top3_markers)))))
cluster.top3markers = gsub('_', '-', cluster.top3markers)
cluster.top3markers <- intersect(cluster.top3markers, rownames(sample2.merge.obj[['RNA']]@data))
cluster.top3markers.mat <- matrix(0, nrow=length(cluster.top3markers), ncol=length(all_markers_final$top3_markers)+1) # add one tumor markers from FACS
rownames(cluster.top3markers.mat) <- cluster.top3markers
colnames(cluster.top3markers.mat) <- c(names(all_markers_final$top3_markers), 'Tumor_FACS_marker')
for (n in names(all_markers_final$top3_markers)) {
    cluster.top3markers.mat[intersect(rownames(sample2.merge.obj[['RNA']]@data), tolower(all_markers_final$top3_markers[[n]])), n] <- 1
}
cluster.top3markers.mat[intersect(rownames(sample2.merge.obj[['RNA']]@data), gsub('_', '-', tolower(all_markers_final$tumor_related$tumor_markers))), 'Tumor_FACS_marker'] <- 1
cluster.top3markers.mat = cluster.top3markers.mat[, colSums(cluster.top3markers.mat)!=0]
library(SingleCellExperiment)
sample2.sce = subset(sample2.merge.obj, features=rownames(cluster.top3markers.mat))
sample2.sce = as.SingleCellExperiment(sample2.sce)
sizeFactors(sample2.sce) <- colSums(assay(sample2.sce))
sample2.sce.size = sizeFactors(sample2.sce)
sample2.sce = sample2.sce[, sample2.sce.size > 0]
sample2.sce.size = sample2.sce.size[sample2.sce.size > 0]

fit <- cellassign(exprs_obj = sample2.sce,
                  marker_gene_info = cluster.top3markers.mat[rownames(sample2.sce), ],
                  s = sample2.sce.size,
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = FALSE)
print(head(celltypes(fit)))

pdf("../results/tumor21_cellassign.pdf")
pheatmap(cellprobs(fit))
dev.off()

cluster.top3markers  <- tolower(unique(c(all_markers_final$tumor_related$tumor_markers,
	 unique(unlist(all_markers_final$top3_markers)))))
cluster.top3markers = gsub('_', '-', cluster.top3markers)
cluster.top3markers <- intersect(cluster.top3markers, rownames(sample3.merge.obj[['RNA']]@data))
cluster.top3markers.mat <- matrix(0, nrow=length(cluster.top3markers), ncol=length(all_markers_final$top3_markers)+1) # add one tumor markers from FACS
rownames(cluster.top3markers.mat) <- cluster.top3markers
colnames(cluster.top3markers.mat) <- c(names(all_markers_final$top3_markers), 'Tumor_FACS_marker')
for (n in names(all_markers_final$top3_markers)) {
    cluster.top3markers.mat[intersect(rownames(sample3.merge.obj[['RNA']]@data), tolower(all_markers_final$top3_markers[[n]])), n] <- 1
}
cluster.top3markers.mat[intersect(rownames(sample3.merge.obj[['RNA']]@data), gsub('_', '-', tolower(all_markers_final$tumor_related$tumor_markers))), 'Tumor_FACS_marker'] <- 1
cluster.top3markers.mat = cluster.top3markers.mat[, colSums(cluster.top3markers.mat)!=0]

library(SingleCellExperiment)
sample3.sce = subset(sample3.merge.obj, features=rownames(cluster.top3markers.mat))
sample3.sce = as.SingleCellExperiment(sample3.sce)
sizeFactors(sample3.sce) <- colSums(assay(sample3.sce))
sample3.sce.size = sizeFactors(sample3.sce)
sample3.sce = sample3.sce[, sample3.sce.size > 0]
sample3.sce.size = sample3.sce.size[sample3.sce.size > 0]

fit <- cellassign(exprs_obj = sample3.sce,
                  marker_gene_info = cluster.top3markers.mat[rownames(sample3.sce), ],
                  s = sample3.sce.size,
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = FALSE)
print(head(celltypes(fit)))

pdf("../results/tumor22_cellassign.pdf")
pheatmap(cellprobs(fit))
dev.off()
