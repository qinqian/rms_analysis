library(Seurat)
library(glue)

tumor = readRDS('../figures/20082_hg19_premrna_tumoronly_res0.8_umap.rds')
table(tumor@meta.data$RNA_snn_res.0.8)

Idents(tumor) = tumor@meta.data$RNA_snn_res.0.8

tumor = subset(tumor, idents=c(11, 9), invert=T)

table(tumor@meta.data$RNA_snn_res.0.8)

pdf(glue('20082_refine_with_clusterlabel.pdf'), width=6, height=4)
tumor$seurat_clusters = tumor$RNA_snn_res.0.8
DimPlot(tumor, reduction='umap', group.by='seurat_clusters', label = TRUE, legend = "none")
dev.off()

tumor$seurat_clusters = droplevels(tumor$seurat_clusters)

library(tidyverse)
source("DEGs_seurat3_sara.R")
library(foreach)
FindDE <- function(x) {
    allcluster <- names((x@meta.data$seurat_clusters) %>% table())
    clusterde <- list()
    for (i in allcluster) {
        print(i)
        print(allcluster[-(as.integer(i)+1)])
        de.up <- get_upregulated_genes_in_clusters(x, i, allcluster[-(as.integer(i)+1)])
        if (nrow(de.up) > 0) {
            de.up$cluster <- i
        }
        clusterde[[i]] <- de.up
    }
    seurat1.de <- do.call('rbind', clusterde) ## still too many genes
    seurat1.de
}

tumor.de <- FindDE(tumor)

tumor.de <- subset(tumor.de, p.adjusted <= 0.01 & (fg_fraction >= 0.1 | bg_fraction >= 0.1))

write.table(tumor.de, file='20082_refine_tumoronly_res0.8.xls', sep='\t', quote=F, col.names=NA)

rerun_cluster <- function(obj) {
    obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(object = obj, selection.method = "vst",
                                mean.function = ExpMean, dispersion.function = LogVMR,
                                x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
    obj <- ScaleData(object = obj, genes.use = rownames(obj),
                     vars.to.regress = c("nUMI", "nGene", "percent.mito", "percent.ribo", "fraction.mouse"),
                     model.use = "linear", use.umi = FALSE) 
    highvar.genes <- head(VariableFeatures(object = obj), 1000)
    obj <- RunPCA(object = obj, features = highvar.genes,
                  npcs = 50, ndims.print = 1:5, nfeatures.print = 1:5,
                  reduction.name = "pca", reduction.key = "PC_", seed.use = 123)
    obj <- FindNeighbors(object = obj, k.param = 20, dims = 1:20, reduction = "pca")
    obj <- FindClusters(obj, resolution=0.8)
    obj <- RunUMAP(obj, reduction.use = "pca", dims = 1:20, reduction.name = "umap")
    obj
}

tumor = rerun_cluster(tumor)

tumor.de <- FindDE(tumor)

saveRDS(tumor, '20082_recluster2_tumor_only.rds')

tumor.de <- subset(tumor.de, p.adjusted <= 0.01 & (fg_fraction >= 0.1 | bg_fraction >= 0.1))
write.table(tumor.de, file='20082_recluster2_tumoronly_res0.8.xls', sep='\t', quote=F, col.names=NA)

pdf(glue('20082_recluster2_with_clusterlabel.pdf'), width=6, height=4)
DimPlot(tumor, reduction='umap', group.by='seurat_clusters', label = TRUE, legend = "none")
dev.off()
