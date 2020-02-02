library(Seurat) 
library(cowplot)
library(data.table)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
library(doubletFinder)
library(Matrix)
library(foreach)

projectPathEris <- c("/data/langenau/human_rms_pdxs/")

names <- c("RH74", "MAST39", "MAST95", "MAST35", "MAST71", "MAST85", "MAST111", "MAST139", "MSK72117", "MSK82489")

type <- c("ERMS","ERMS","ARSM-P7","ARMS-P3","ERMS","ERMS","ERMS","ERMS","ARMS-P7","ARMS-P3") 

samples <- c("20190618_PDX1_RH74_5Kcells_hg19","20190618_PDX2_MAST39_5Kcells_hg19","20190618_PDX3_MAST95_10Kcells_hg19", 
             "20190418_MAST35_5Kcells_hg19","20190418_MAST71_5Kcells_hg19","20190418_MAST85_5Kcells_hg19", 
             "20190418_MAST111_5Kcells_hg19","20190418_MAST139_5Kcells_hg19", 
             "20190617_MSK72117_5Kcells_hg19","20190617_MSK82489_5Kcells_hg19")

files.mapstats <- c("20190524_mapstats_PDX1_RH74_5Kcells_hg19_mm10_counts_uniq.txt",
                    "20190524_mapstats_PDX2_MAST39_5Kcells_hg19_mm10_counts_uniq.txt",
                    "20190524_mapstats_PDX3_MAST95_10Kcells_hg19_mm10_counts_uniq.txt",
                    "20190524_mapstats_MAST35_5Kcells_hg19_mm10_counts_uniq.txt",
                    "20190620_mapstats_MAST71_5Kcells_hg19_mm10_counts_uniq.txt",
                    "20190620_mapstats_MAST85_5Kcells_hg19_mm10_counts_uniq.txt",
                    "20190524_mapstats_MAST111_5Kcells_hg19_mm10_counts_uniq.txt",
                    "20190620_mapstats_MAST139_5Kcells_hg19_mm10_counts_uniq.txt",
                    "20190619_mapstats_MSK72117_5Kcells_hg19_mm10_counts_uniq.txt",
                    "20190619_mapstats_MSK82489_5Kcells_hg19_mm10_counts_uniq.txt")

lib.size <- list()
nr.genes <- list()
perc.mouse <- list()
perc.mito <- list()
perc.ribo <- list()

for (i in 1:length(samples)){
  mat <- Read10X(data.dir = paste0(projectPathEris, samples[i], "/outs/filtered_feature_bc_matrix/"))

  seurat.obj <- CreateSeuratObject(counts = mat, project = names[i], min.cells = 3, min.features = 1) 
  
  # add info to metadata 
  # 1) percentage of mitochondrial genes
  seurat.obj[["percent.mito"]] <- PercentageFeatureSet(object = seurat.obj, pattern = "^MT-")
  # 2) percentage of ribosomal genes
  seurat.obj[["percent.ribo"]] <- PercentageFeatureSet(object = seurat.obj, pattern = "^RP[SL][[:digit:]]")
  # 3) percentage of human and mouse reads
  map.counts <- read.delim(paste0(projectPathEris, "mapstats/", files.mapstats[i]), header = F, sep = " ", stringsAsFactors = F)
  df <- map.counts %>%
      mutate("cell" = sapply(strsplit(sapply(strsplit(V1,"-"),'[',i=1),":"),'[',i=3),
             "map" = substr(V1, 25, 25),
             "count" = V2) %>% select(-V1, -V2) %>% filter(cell %in% rownames(seurat.obj@meta.data)) 
  totals <- df %>% group_by(cell) %>% summarise("total" = sum(count))
  df2 <- inner_join(df, totals, by = "cell") %>% mutate("fraction" = count/total)
  df.filt <- df2 %>% filter(map=="h") %>% select(cell, count, fraction)
  vec <- as.vector(df.filt$count)
  names(vec) <- df.filt$cell
  seurat.obj$counts.human <- vec
  vec <- as.vector(df.filt$fraction)
  names(vec) <- df.filt$cell
  seurat.obj$fraction.human <- vec
  df.filt <- df2 %>% filter(map=="m") %>% select(cell, count, fraction) 
  vec <- as.vector(df.filt$count)
  names(vec) <- df.filt$cell
  seurat.obj$counts.mouse <- vec
  vec <- as.vector(df.filt$fraction)
  names(vec) <- df.filt$cell
  seurat.obj$fraction.mouse <- vec
  # 4) library size for filtering
  seurat.obj[["lib.size.10k"]] <- log10(seurat.obj[["nCount_RNA"]]/1e4+1)

  # QC plots for sample-specific filtering
  lib.size[[i]] <- log10(seurat.obj[["nCount_RNA"]][,1]/1e4+1)
  nr.genes[[i]] <- seurat.obj[["nFeature_RNA"]][,1]
  perc.mouse[[i]] <- 100*seurat.obj[["fraction.mouse"]][,1]
  perc.mito[[i]] <- seurat.obj[["percent.mito"]][,1]
  perc.ribo[[i]] <- seurat.obj[["percent.ribo"]][,1]
  # mean library size
  mean.lib.size <- mean(lib.size[[i]])
  # median library size
  med.lib.size <- median(lib.size[[i]])
  # mean number of genes
  mean.nr.genes <- mean(nr.genes[[i]])
  # median number of genes
  med.nr.genes <- median(nr.genes[[i]])
  # mean number of mouse genes
  mean.mous.genes <- mean(perc.mouse[[i]])
  # median number of mouse genes
  med.mous.genes <- median(perc.mouse[[i]])
  # mean number of mitochondrial genes
  mean.mito.genes <- mean(perc.mito[[i]])
  # median number of mitochondrial genes
  med.mito.genes <- median(perc.mito[[i]])
  # mean number of ribosomal genes
  mean.ribo.genes <- mean(perc.ribo[[i]])
  # median number of ribosomal genes
  med.ribo.genes <- median(perc.ribo[[i]])
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_QC_plot_1_", names[i], ".pdf"), width = 12, height = 8)
  par(mfrow=c(2,3))
  hist(lib.size[[i]], xlab="Library size (10^4)", main=paste0("mean ", round(mean.lib.size,2)),
       breaks=20, col="grey80", ylab="Number of cells")
  abline(v=mean.lib.size, col="red", lwd=2, lty=1)
  abline(v=med.lib.size, col="red", lwd=2, lty=2)
  hist(nr.genes[[i]], xlab="Number of expressed genes", main=paste0("mean ", round(mean.nr.genes,2)),
       breaks=20, col="grey80", ylab="Number of cells")
  abline(v=mean.nr.genes, col="red", lwd=2, lty=1)
  abline(v=med.nr.genes, col="red", lwd=2, lty=2)
  hist(perc.mouse[[i]], xlab="% Mouse reads", main=paste0("mean ", round(mean.mous.genes,2)),
       breaks=20, col="grey80", ylab="Number of cells")
  abline(v=mean.mous.genes, col="red", lwd=2, lty=1)
  abline(v=med.mous.genes, col="red", lwd=2, lty=2)
  hist(perc.mito[[i]], xlab="% Mitochondrial reads", main=paste0("mean ", round(mean.mito.genes,2)),
       breaks=20, col="grey80", ylab="Number of cells")
  abline(v=mean.mito.genes, col="red", lwd=2, lty=1)
  abline(v=med.mito.genes, col="red", lwd=2, lty=2)
  hist(perc.ribo[[i]], xlab="% Ribosomal reads", main=paste0("mean ", round(mean.ribo.genes,2)),
       breaks=20, col="grey80", ylab="Number of cells")
  abline(v=mean.ribo.genes, col="red", lwd=2, lty=1)
  abline(v=med.ribo.genes, col="red", lwd=2, lty=2)
  dev.off()

  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_QC_plot_2_", names[i], ".pdf"), width = 16, height = 6)
    VlnPlot(object = seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo","fraction.human","fraction.mouse"), 
            pt.size = 0.1, ncol = 6)
  dev.off()

  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_QC_plot_3_", names[i], ".pdf"), width = 12, height = 8)
    plot1 <- FeatureScatter(object = seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    plot2 <- FeatureScatter(object = seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mito")
    plot3 <- FeatureScatter(object = seurat.obj, feature1 = "nFeature_RNA", feature2 = "percent.mito")
    plot4 <- FeatureScatter(object = seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.ribo")
    CombinePlots(plots = list(plot1, plot4, plot2, plot3))
  dev.off()
  
  #  sample-specific thresholds
  ## to filter cells with log-library size > 3 mads below the median log-library size
  # median absolute deviation of the library size
  mad.lib.size <- mad(lib.size[[i]])
  # threshold for minimum library size 
  min.lib.size <- med.lib.size-3*mad.lib.size
  
  # to filter cells with log-transformed number of genes 3 mads below the median
  # median absolute deviation of the number of genes
  mad.nr.genes <- mad(nr.genes[[i]])
  # threshold for minimum number of genes
  min.nr.genes <- max(med.nr.genes-3*mad.nr.genes, 1000)

  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_QC_plot_1_", names[i], ".pdf"), width = 12, height = 8)
  par(mfrow=c(2,3))
  hist(lib.size[[i]], xlab="Library size (10^4)", main=paste0("mean ", round(mean.lib.size,2)),
       breaks=20, col="grey80", ylab="Number of cells")
  abline(v=mean.lib.size, col="red", lwd=2, lty=1)
  abline(v=med.lib.size, col="red", lwd=2, lty=2)
  hist(nr.genes[[i]], xlab="Number of expressed genes", main=paste0("mean ", round(mean.nr.genes,2)),
       breaks=20, col="grey80", ylab="Number of cells")
  abline(v=mean.nr.genes, col="red", lwd=2, lty=1)
  abline(v=med.nr.genes, col="red", lwd=2, lty=2)
  hist(perc.mouse[[i]], xlab="% Mouse reads", main=paste0("mean ", round(mean.mous.genes,2)),
       breaks=20, col="grey80", ylab="Number of cells")
  abline(v=mean.mous.genes, col="red", lwd=2, lty=1)
  abline(v=med.mous.genes, col="red", lwd=2, lty=2)
  hist(perc.mito[[i]], xlab="% Mitochondrial reads", main=paste0("mean ", round(mean.mito.genes,2)),
       breaks=20, col="grey80", ylab="Number of cells")
  abline(v=mean.mito.genes, col="red", lwd=2, lty=1)
  abline(v=med.mito.genes, col="red", lwd=2, lty=2)
  hist(perc.ribo[[i]], xlab="% Ribosomal reads", main=paste0("mean ", round(mean.ribo.genes,2)),
       breaks=20, col="grey80", ylab="Number of cells")
  abline(v=mean.ribo.genes, col="red", lwd=2, lty=1)
  abline(v=med.ribo.genes, col="red", lwd=2, lty=2)
  dev.off()
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_QC_plot_2_", names[i], ".pdf"), width = 16, height = 6)
  VlnPlot(object = seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo","fraction.human","fraction.mouse"), 
          pt.size = 0.1, ncol = 6)
  dev.off()
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_QC_plot_3_", names[i], ".pdf"), width = 12, height = 8)
  plot1 <- FeatureScatter(object = seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot2 <- FeatureScatter(object = seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.mito")
  plot3 <- FeatureScatter(object = seurat.obj, feature1 = "nFeature_RNA", feature2 = "percent.mito")
  plot4 <- FeatureScatter(object = seurat.obj, feature1 = "nCount_RNA", feature2 = "percent.ribo")
  CombinePlots(plots = list(plot1, plot4, plot2, plot3))
  dev.off()
  
  seurat.obj <- subset(x = seurat.obj, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mito < 20 & fraction.mouse< 0.05 ) 
  
   seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
  seurat.obj <- FindVariableFeatures(object = seurat.obj, selection.method = "vst",
                                     mean.function = ExpMean, dispersion.function = LogVMR,
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)

  highvar.genes <- head(VariableFeatures(object = seurat.obj), 1000)
  top10.highvar.genes <- head(VariableFeatures(object = seurat.obj), 10)
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_highvar_features_plot_", names[i], ".pdf"), width = 12, height = 6)
  plot1 <- VariableFeaturePlot(object = seurat.obj)
  plot2 <- LabelPoints(plot = plot1, points = top10.highvar.genes, repel = TRUE)
  CombinePlots(plots = list(plot1, plot2))
  dev.off()
  
  seurat.obj <- ScaleData(object = seurat.obj, genes.use = all.genes, 
                          vars.to.regress = c("nUMI", "nGene", "percent.mito", "percent.ribo", "fraction.mouse"), 
                          model.use = "linear", use.umi = FALSE) 
  
  seurat.obj <- RunPCA(object = seurat.obj, features = highvar.genes,
                       npcs = 50, ndims.print = 1:5, nfeatures.print = 1:5,
                       reduction.name = "pca", reduction.key = "PC_", seed.use = 123)
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_PCA_plot_", names[i], ".pdf"), width = 6, height = 6)
  DimPlot(object = seurat.obj, reduction = "pca")
  dev.off()
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_PCA_loadings_plot_", names[i], ".pdf"), width = 10, height = 12)
  VizDimLoadings(object = seurat.obj, dims = 1:4, reduction = "pca")
  dev.off()
  
  seurat.obj <- JackStraw(object = seurat.obj, num.replicate = 100)
  seurat.obj <- ScoreJackStraw(object = seurat.obj, dims = 1:20)
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_PCA_jackstraw_plot_", names[i], ".pdf"), width = 6, height = 6)
  JackStrawPlot(object = seurat.obj, dims = 1:20)
  dev.off()
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_PCA_elbow_plot_", names[i], ".pdf"), width = 6, height = 6)
  ElbowPlot(object = seurat.obj, ndims = 50) 
  dev.off()

  seurat.obj <- FindNeighbors(object = seurat.obj, k.param = 20, dims = 1:20, reduction = "pca")

  seurat.obj <- FindClusters(object = seurat.obj, reduction.type = "pca", dims.use = 1:20, 
                             algorithm = 1, 
                             resolution = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
                             random.seed = 123)
                     
  seurat.obj <- RunUMAP(seurat.obj, reduction.use = "pca", dims = 1:20, reduction.name = "umap")

  seurat.obj <- RunTSNE(seurat.obj, reduction.use = "pca", dims.use = 1:20, reduction.name = "tsne")
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_tSNE_plot_20PCs_color-nUMI_", names[i], ".pdf"), width = 6, height = 6)
  FeaturePlot(seurat.obj, features = "nCount_RNA", reduction = "tsne", cols = c("lightgrey", "red"))
  dev.off()
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_tSNE_plot_20PCs_color-nFeature_", names[i], ".pdf"), width = 6, height = 6)
  FeaturePlot(seurat.obj, features = "nFeature_RNA", reduction = "tsne", cols = c("lightgrey", "red"))
  dev.off()
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_tSNE_plot_20PCs_color-mito_", names[i], ".pdf"), width = 6, height = 6)
  FeaturePlot(seurat.obj, features = "percent.mito", reduction = "tsne", cols = c("lightgrey", "red"))
  dev.off()
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_tSNE_plot_20PCs_color-ribo_", names[i], ".pdf"), width = 6, height = 6)
  FeaturePlot(seurat.obj, features = "percent.ribo", reduction = "tsne", cols = c("lightgrey", "red"))
  dev.off()
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_tSNE_plot_20PCs_color-mouse_", names[i], ".pdf"), width = 6, height = 6)
  FeaturePlot(seurat.obj, features = "fraction.mouse", reduction = "tsne", cols = c("lightgrey", "red"))
  dev.off()
  
  for (j in c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)){
    pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_UMAP_plot_20PCs_res", j, "_", names[i], ".pdf"), width = 6, height = 6)
    UMAPPlot(seurat.obj, reduction = "umap", group.by = paste0("RNA_snn_res.", j))
    dev.off()
    
    pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_tSNE_plot_20PCs_res", j, "_", names[i], ".pdf"), width = 6, height = 6)
    TSNEPlot(seurat.obj, reduction.use = "tsne", group.by = paste0("RNA_snn_res.", j))
    dev.off()
  }
  
  pdf(paste0(projectPathEris, "plots/", format(Sys.time(), "%Y%m%d"), "_human_mouse_genes_", names[i], ".pdf"), width = 18, height = 12)
  FeaturePlot(object = seurat.obj, reduction = "tsne", cols = c("lightgrey", "red"), order = TRUE, 
              features = c("CTSS", "CD68", "HEXB", "MAFB", "LGALS3", "AIF1", "SRGN", "TYROBP", "SLC15A3", "BCL2A1", "GNGT2"))
  dev.off()
  
  saveRDS(seurat.obj, file = paste0(projectPathEris, "seurat_objects/", format(Sys.time(), "%Y%m%d"), "_seurat-object_", names[i], ".rds"))
}  

