## library(cowplot)
## library(data.table)
## library(tidyverse)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(Matrix)
library(foreach)
library(readr)

system('mkdir -p results/seurat_sara')

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--mixtureobj', dest='mixture', metavar='N', type="character", nargs='+')
    parser$add_argument('--seuratobj', dest='seurat', metavar='N', type="character", nargs='+')
    parser$add_argument('--label', dest='label', default='')
    parser$add_argument('--doublet', dest='doublet', default='')
    args = parser$parse_args()
    args
}

args = get_args()
if (length(args$seurat) == 0 || args$label == '') {
    cat('empty argument, exit..')
    q()
}

names <- args$label
type <- c("ERMS")
samples <- args$seurat
mixture.samples <- args$mixture
doublet <- args$doublet

## names <- 'mast139_muscle_plus'
## type  <- 'ERMS'
## samples <- 'MAST139_Muscle_plus/outs/filtered_feature_bc_matrix'
## mixture.samples <-'MAST139_Muscle_plus_mixture/outs/filtered_feature_bc_matrix'

## names <- 'mast139_muscle_minus'
## samples <- 'MAST139_Muscle_minus/outs/filtered_feature_bc_matrix'
## mixture.samples <-'MAST139_Muscle_minus_mixture/outs/filtered_feature_bc_matrix'

## names <- '20082'
## samples <- '20082_hg19_premrna/outs/filtered_feature_bc_matrix'
## mixture.samples <-'20082_hg19_premrna/outs/filtered_feature_bc_matrix'
## doublet <- '../results/20082_hg19_premrna_doublet_doublet.csv'

## names <- 'C12SC1'
## samples <- '/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/C12SC1_hg19/outs/filtered_feature_bc_matrix'
## mixture.samples <- '/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/C12SC1_mixture/outs/filtered_feature_bc_matrix'
## doublet <- '../results/C12SC1_hg19_doublet_doublet.csv'

mat <- Read10X(samples[1])

colnames(mat) <- gsub('-1', '', colnames(mat))


print(doublet)
if (doublet!="") {
    doublet = read.csv(doublet, stringsAsFactors=F)
    cells = gsub('-1', '', doublet[,2])
    non_doublet <- rep(T, nrow(doublet))
    non_doublet[doublet$prediction=='True'] = F
    cells = cells[non_doublet]
    mat = mat[, colnames(mat)%in%cells]
}

seurat.obj <- CreateSeuratObject(counts = mat, project=names[1], min.cells = 3, min.features = 1)

seurat.obj[["percent.mito"]] <- PercentageFeatureSet(object = seurat.obj, pattern = "^MT-")
## 2) percentage of ribosomal genes
seurat.obj[["percent.ribo"]] <- PercentageFeatureSet(object = seurat.obj, pattern = "^RP[SL][[:digit:]]")

print(111)
## calculate mouse ratio by sara
## map.counts <- read.delim(paste0(projectPathEris, "mapstats/", files.mapstats[i]), header = F, sep = " ", stringsAsFactors = F)
## df <- map.counts %>%
##   mutate("cell" = sapply(strsplit(sapply(strsplit(V1,"-"),'[',i=1),":"),'[',i=3),
##          "map" = substr(V1, 25, 25),
##          "count" = V2) %>% select(-V1, -V2) %>% filter(cell %in% rownames(seurat.obj@meta.data)) 
## totals <- df %>% group_by(cell) %>% summarise("total" = sum(count))
## df2 <- inner_join(df, totals, by = "cell") %>% mutate("fraction" = count/total)
## df.filt <- df2 %>% filter(map=="h") %>% select(cell, count, fraction)
## vec <- as.vector(df.filt$count)
## names(vec) <- df.filt$cell
## seurat.obj$counts.human <- vec
## vec <- as.vector(df.filt$fraction)
## names(vec) <- df.filt$cell
## seurat.obj$fraction.human <- vec
## df.filt <- df2 %>% filter(map=="m") %>% select(cell, count, fraction) 
## vec <- as.vector(df.filt$count)
## names(vec) <- df.filt$cell
## seurat.obj$counts.mouse <- vec
## vec <- as.vector(df.filt$fraction)
## names(vec) <- df.filt$cell
## seurat.obj$fraction.mouse <- vec

## replace with matrix operation 
seurat.obj[["lib.size.10k"]] <- log10(seurat.obj[["nCount_RNA"]]/1e4+1)

if (!is.na(mixture.samples)) {
    mix.mat <- Read10X(mixture.samples[1])
    colnames(mix.mat) <- gsub('-1', '', colnames(mix.mat))

    fraction.mouse = Matrix::colSums(mix.mat[grepl('mm10', rownames(mix.mat)), ])/Matrix::colSums(mix.mat)
    ## 4) library size for filtering
    seurat.obj[['fraction.mouse']] = fraction.mouse[match(colnames(seurat.obj), names(fraction.mouse))]
}
    
pdf(paste0('results/seurat_sara/', names[1], "_QC_plot", ".pdf"), width = 18, height = 8)
par(mfrow=c(2,3))
if (!is.na(mixture.samples)) {
VlnPlot(object = seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo","fraction.mouse"), 
        pt.size = 0.1, ncol = 6)
} else {
VlnPlot(object = seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo"), 
        pt.size = 0.1, ncol = 6)
}
hist(seurat.obj@meta.data$lib.size.10k, xlab="Library size (10^4)", main=paste0("mean ", round(mean(seurat.obj@meta.data$lib.size.10k),2)),
     breaks=20, col="grey80", ylab="Number of cells")
abline(v=mean(seurat.obj@meta.data$lib.size.10k), col="red", lwd=2, lty=1)
abline(v=median(seurat.obj@meta.data$lib.size.10k), col="red", lwd=2, lty=2)
hist(seurat.obj@meta.data$nFeature_RNA, xlab="Number of expressed genes", main=paste0("mean ", round(mean(seurat.obj@meta.data$nFeature_RNA), 2)),
     breaks=20, col="grey80", ylab="Number of cells")
abline(v=mean(seurat.obj@meta.data$nFeature_RNA), col="red", lwd=2, lty=1)
abline(v=median(seurat.obj@meta.data$nFeature_RNA), col="red", lwd=2, lty=2)
if (!is.na(mixture.samples)) {
    hist(seurat.obj@meta.data$fraction.mouse, xlab="% Mouse reads", main=paste0("mean ", round(mean(seurat.obj@meta.data$fraction.mouse),2)),
         breaks=20, col="grey80", ylab="Number of cells")
    abline(v=mean(seurat.obj@meta.data$fraction.mouse), col="red", lwd=2, lty=1)
    abline(v=median(seurat.obj@meta.data$fraction.mouse), col="red", lwd=2, lty=2)
}
hist(seurat.obj@meta.data$percent.mito, xlab="% Mitochondrial reads", main=paste0("mean ", round(mean(seurat.obj@meta.data$percent.mito),2)),
     breaks=20, col="grey80", ylab="Number of cells")
abline(v=mean(seurat.obj@meta.data$percent.mito), col="red", lwd=2, lty=1)
abline(v=median(seurat.obj@meta.data$percent.mito), col="red", lwd=2, lty=2)
hist(seurat.obj$percent.ribo, xlab="% Ribosomal reads", main=paste0("mean ", round(mean(seurat.obj$percent.ribo),2)),
     breaks=20, col="grey80", ylab="Number of cells")
abline(v=mean(seurat.obj$percent.ribo), col="red", lwd=2, lty=1)
abline(v=median(seurat.obj$percent.ribo), col="red", lwd=2, lty=2)
dev.off()

if (!is.na(mixture.samples)) {
    seurat.obj <- subset(x = seurat.obj, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mito < 20)

    print(max(seurat.obj@meta.data$fraction.mouse))
    print(summary(seurat.obj@meta.data$fraction.mouse))

    seurat.obj2 <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat.obj2 <- FindVariableFeatures(object = seurat.obj2, selection.method = "vst",
                                       mean.function = ExpMean, dispersion.function = LogVMR,
                                       x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
    if (!is.na(mixture.samples)) {
        seurat.obj2 <- ScaleData(object = seurat.obj2, genes.use = rownames(seurat.obj2), 
                                vars.to.regress = c("nUMI", "nGene", "percent.mito", "percent.ribo", "fraction.mouse"),
                                model.use = "linear", use.umi = FALSE) 
    } else {
        seurat.obj2 <- ScaleData(object = seurat.obj2, genes.use = rownames(seurat.obj2), 
                                vars.to.regress = c("nUMI", "nGene", "percent.mito", "percent.ribo"), 
                                model.use = "linear", use.umi = FALSE) 
    }
    highvar.genes <- head(VariableFeatures(object = seurat.obj2), 1000)
    seurat.obj2 <- RunPCA(object = seurat.obj2, features = highvar.genes,
                         npcs = 50, ndims.print = 1:5, nfeatures.print = 1:5,
                         reduction.name = "pca", reduction.key = "PC_", seed.use = 123)
    seurat.obj2 <- JackStraw(object = seurat.obj2, num.replicate = 100)
    seurat.obj2 <- ScoreJackStraw(object = seurat.obj2, dims = 1:20)
    seurat.obj2 <- FindNeighbors(object = seurat.obj2, k.param = 20, dims = 1:20, reduction = "pca")
    seurat.obj2 <- FindClusters(object = seurat.obj2, reduction.type = "pca", dims.use = 1:20, 
                               algorithm = 1, 
                               resolution = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
                               random.seed = 123)
    seurat.obj2 <- RunUMAP(seurat.obj2, reduction.use = "pca", dims = 1:20, reduction.name = "umap")
    seurat.obj2@meta.data$mouse <- seurat.obj2@meta.data$fraction.mouse >= 0.05

    pdf(paste0('results/seurat_sara/', names[1], "_UMAP_with_mouse", ".pdf"), width = 8.2, height = 5)
    p1 = DimPlot(seurat.obj2, reduction = "umap", group.by='mouse')
    print(p1)
    dev.off()

    seurat.obj <- subset(seurat.obj, cells=names(which(fraction.mouse < 0.05)))
} else {
    seurat.obj <- subset(x = seurat.obj, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mito < 20)
}

seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)

seurat.obj <- FindVariableFeatures(object = seurat.obj, selection.method = "vst",
                                   mean.function = ExpMean, dispersion.function = LogVMR,
                                   x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
if (!is.na(mixture.samples)) {
seurat.obj <- ScaleData(object = seurat.obj, genes.use = rownames(seurat.obj), 
                        vars.to.regress = c("nUMI", "nGene", "percent.mito", "percent.ribo", "fraction.mouse"),
                        model.use = "linear", use.umi = FALSE) 
} else {
seurat.obj <- ScaleData(object = seurat.obj, genes.use = rownames(seurat.obj), 
                        vars.to.regress = c("nUMI", "nGene", "percent.mito", "percent.ribo"), 
                        model.use = "linear", use.umi = FALSE) 
}

highvar.genes <- head(VariableFeatures(object = seurat.obj), 1000)

seurat.obj <- RunPCA(object = seurat.obj, features = highvar.genes,
                     npcs = 50, ndims.print = 1:5, nfeatures.print = 1:5,
                     reduction.name = "pca", reduction.key = "PC_", seed.use = 123)

seurat.obj <- JackStraw(object = seurat.obj, num.replicate = 100)
seurat.obj <- ScoreJackStraw(object = seurat.obj, dims = 1:20)

seurat.obj <- FindNeighbors(object = seurat.obj, k.param = 20, dims = 1:20, reduction = "pca")

seurat.obj <- FindClusters(object = seurat.obj, reduction.type = "pca", dims.use = 1:20, 
                           algorithm = 1, 
                           resolution = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
                           random.seed = 123)
                     
seurat.obj <- RunUMAP(seurat.obj, reduction.use = "pca", dims = 1:20, reduction.name = "umap")
seurat.obj <- RunTSNE(seurat.obj, reduction.use = "pca", dims.use = 1:20, reduction.name = "tsne")

saveRDS(seurat.obj, file = paste0('results/seurat_sara/', names[1], "_seurat-object", ".rds"))

pdf(paste0('results/seurat_sara/', names[1], "_seurat-object_QC", ".pdf"), width = 6, height = 6)
FeaturePlot(seurat.obj, features = "nCount_RNA", reduction = "tsne", cols = c("lightgrey", "red"))
FeaturePlot(seurat.obj, features = "nFeature_RNA", reduction = "tsne", cols = c("lightgrey", "red"))
FeaturePlot(seurat.obj, features = "percent.mito", reduction = "tsne", cols = c("lightgrey", "red"))
FeaturePlot(seurat.obj, features = "percent.ribo", reduction = "tsne", cols = c("lightgrey", "red"))
if (!is.na(mixture.samples)) {
FeaturePlot(seurat.obj, features = "fraction.mouse", reduction = "tsne", cols = c("lightgrey", "red"))
}
dev.off()

for (j in c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2)) {
    pdf(paste0('results/seurat_sara/', names[1], "_UMAP_plot_20PCs_res", j, ".pdf"), width = 12, height = 6)
    p1 = DimPlot(seurat.obj, reduction = "umap", group.by = paste0("RNA_snn_res.", j))
    p2 = DimPlot(seurat.obj, reduction = "tsne", group.by = paste0("RNA_snn_res.", j))
    print(CombinePlots(plots=list(p1, p2)))
    dev.off()
}

pdf(paste0('results/seurat_sara/', names[1], "_human_mouse_genes_", ".pdf"), width = 18, height = 12)
FeaturePlot(object = seurat.obj, reduction = "tsne", cols = c("lightgrey", "red"), order = TRUE, 
            features = c("CTSS", "CD68", "HEXB", "MAFB", "LGALS3", "AIF1", "SRGN", "TYROBP", "SLC15A3", "BCL2A1", "GNGT2"))
dev.off()
