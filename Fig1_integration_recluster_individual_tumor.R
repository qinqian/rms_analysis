library(Seurat)

primary.obj <- Reduce("merge", list(
                                readRDS('../results/seurat_sara/20082_hg19_premrna_seurat-object.rds'),
                                readRDS('../results/seurat_sara/21202_hg19_premrna_seurat-object.rds'),
                                readRDS('../results/seurat_sara/20696_seurat-object.rds'),
                                readRDS('../results/seurat_sara/29806_hg19_premrna_seurat-object.rds'))
                      )


p20082 <- readRDS('../results/seurat_sara/20082_hg19_premrna_seurat-object.rds')
Idents(p20082) <- p20082@meta.data$RNA_snn_res.0.8
p20082.tumor = subset(p20082, idents=c(0, 10, 11, 13, 14, 15, 2, 3, 4, 6, 7))

p21202 <- readRDS('../results/seurat_sara/21202_hg19_premrna_seurat-object.rds')
Idents(p21202) <- p21202@meta.data$RNA_snn_res.0.8
p21202.tumor = subset(p21202, idents=c(10, 12, 13, 15, 9), invert=T)

p20696 <- readRDS('../results/seurat_sara/20696_seurat-object.rds')
Idents(p20696) <- p20696@meta.data$RNA_snn_res.0.8
p20696.tumor <- subset(p20696, idents=c(0, 1, 10, 2, 4, 7, 8))

p29806 <- readRDS('../results/seurat_sara/29806_hg19_premrna_seurat-object.rds')
Idents(p29806) <- p29806@meta.data$RNA_snn_res.0.8
p29806.tumor <- subset(p29806, idents=c(12, 13, 14, 9), invert=T)

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

for (tumor in list(p20082.tumor, p21202.tumor, p20696.tumor, p29806.tumor)) {
    tumor <- rerun_cluster(tumor)
    ## TODO:
    saveRDS(tumor, paste0('../figures/', tumor@meta.data$orig.ident[1], '_tumoronly_res0.8_umap.rds'))
    ## pdf(paste0('../figures/', tumor@meta.data$orig.ident[1], '_tumoronly_res0.8_umap.pdf'))
    ## print(DimPlot(tumor, reduction='umap', label=T, group.by='seurat_clusters'))
    ## dev.off()
    ## tumor.de <- FindDE(tumor)
    ## tumor.de <- subset(tumor.de, p.adjusted <= 0.01 & (fg_fraction >= 0.1 | bg_fraction >= 0.1))
    ## write.table(tumor.de, file=paste0('../results/seurat_sara/', tumor@meta.data$orig.ident[1], '_tumoronly_res0.8.xls'), sep='\t', quote=F, col.names=NA)
}

library(ComplexHeatmap)
gene.modules <- Sys.glob('lisa_script_modules/*symbols')
gene.list <- lapply(gene.modules, scan, what='')
names(gene.list) <- basename(gsub('.symbols', '', gene.modules))
names(gene.list) <- c("EMT", "G1S", "G2M", "Histone", "Hypoxia", "INTERFERON", "MUSCLE", "TNFA")

list2df <- function(x) {
    ylist <- list()
    for (y in names(x)) {
        ylist[[y]] = data.frame(module=y, gene=x[[y]])
    }
    do.call('rbind', ylist)
}
gene.list.copy = list2df(gene.list)

library(ggplot2)

for (label in c('20082', '21202', '20696', '29806')) {
    tumor = readRDS(Sys.glob(paste0('../figures/', label, '*_tumoronly_res0.8_umap.rds')))
    gene.list = gene.list.copy
    gene.list[, 2] <- as.vector(gene.list[, 2])
    gene.list <- gene.list[gene.list[, 2] %in% rownames(tumor$RNA@data), ]
    df = read.table(paste0('../results/seurat_sara/', tumor@meta.data$orig.ident[1], '_tumoronly_res0.8.xls'))
    ## gene.list <- gene.list[gene.list[, 2] %in% rownames(tumor$RNA@counts), ]
    ## gene.list <- gene.list[gene.list[, 2] %in% rownames(tumor$RNA@data), ]
    sortcells <- sort(tumor$seurat_clusters)
    heatdata  <- as.matrix(tumor$RNA@data)[gene.list[, 2], order(tumor$seurat_clusters)]
    heatdata.cor <- as.matrix(tumor$RNA@data[df$genename, order(tumor$seurat_clusters)])
    annrow <- gene.list[, 1, drop=F]
    rownames(annrow) <- gene.list[, 2]
    anncol <- data.frame(cluster=sortcells)
    metacolors <- c(rgb(166, 166, 166, maxColorValue = 255),
                    rgb(241, 149, 69, maxColorValue = 255),
                    rgb(103, 35,  102, maxColorValue = 255),
                    rgb(52, 101, 252, maxColorValue = 255),
                    rgb(242, 242, 242, maxColorValue = 255),
                    rgb(52, 101, 252, maxColorValue = 255),
                    rgb(233, 63,  51, maxColorValue = 255),
                    rgb(65, 129,  7, maxColorValue = 255),
                    rgb(52, 101, 252, maxColorValue = 255),
                    rgb(253, 247, 49, maxColorValue = 255),
                    'purple')
    metalabels <- c("GROUND", "Hypoxia", "EMT",  "G1S", "UNASSIGNED",
                    "G2M",  "MUSCLE", "INTERFERON", "PROLIF",
                    "Histone", "TNFA")
    names(metacolors) <- metalabels
    ann_colors <- list(
        module = metacolors
    )
    ## selection <- apply(heatdata, 1, function(x) sum(x>0)/length(x) > 0.1)
    ## heatdata <- heatdata[selection, ,drop=F]
    get_correlated_variable_genes = function(mat, n = nrow(mat), cor_cutoff = 0, n_cutoff = 0) {
        ind = order(apply(mat, 1, function(x) {
            q = quantile(x, c(0.05, 0.95))
            x = x[x < q[1] & x > q[2]]
            var(x)/mean(x)
        }), decreasing = TRUE)[1:n]
        mat2 = mat[ind, , drop = FALSE]
        dt = cor(t(mat2), method = "spearman")
        diag(dt) = 0
        dt[abs(dt) < cor_cutoff] = 0
        dt[dt < 0] = -1
        dt[dt > 0] = 1
        i = colSums(abs(dt)) > n_cutoff
        mat3 = mat2[i, ,drop = FALSE]
        return(mat3)
    }
    ## mat = get_correlated_variable_genes(heatdata, cor_cutoff = 0.1, n_cutoff = 0)
    mat <- heatdata
    mat2 = t(apply(mat, 1, function(x) {
        q10 <- quantile(x, 0.001)
        q90 <- quantile(x, 0.999)
        x[x < q10] <- q10
        x[x > q90] <- q90
        scale(x)
    }))
    colnames(mat2) = colnames(mat)
    base_mean = rowMeans(t(mat))
    mat.cor = t(apply(heatdata.cor, 1, function(x) {
        q10 <- quantile(x, 0.001)
        q90 <- quantile(x, 0.999)
        x[x < q10] <- q10
        x[x > q90] <- q90
        scale(x)
    }))
    library(GetoptLong)
    library(RColorBrewer)
    library(ComplexHeatmap)
    library(circlize)
    if(length(unique(anncol[,1]))<=12) {
        ha = structure(brewer.pal(length(unique(anncol[, 1])), "Set3"),
                       names=levels(anncol[, 1]))
    } else {
        ha = structure(c(brewer.pal(name="Dark2", n = round(length(unique(anncol[, 1]))/2)), brewer.pal(name="Paired", n = length(unique(anncol[, 1])) - round(length(unique(anncol[, 1]))/2))),
                       names=levels(anncol[, 1]))
    }
    ha2 = HeatmapAnnotation(cluster=anncol[, 1],
                            col=list(cluster=ha), show_legend=F)
    ha4 = HeatmapAnnotation(modules=as.vector(gene.list[,1]), # as.vector(gene.list[selection, 1]),
                            col=list(modules=ann_colors$module), show_legend=T)
    ht_list = Heatmap(t(mat2), col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                      name = "scaled_expr", column_title = qq("relative expression for @{ncol(mat)} cells"),
                      show_column_dend = FALSE, cluster_rows=F, cluster_columns=F,
                      show_column_names = FALSE, width = unit(8, "cm"), show_row_dend = FALSE,
                      heatmap_legend_param = list(title = "Scaled expr"), top=ha4) +
        Heatmap(anncol[,1], col=ha, width = unit(0.7, "cm"), name='clusters') +
        Heatmap(cor(mat.cor), name = "cor",
                col = colorRamp2(c(-0.3, 0, 0.3), c("blue", "white", "red")),
                show_row_names = FALSE, show_column_names = FALSE, show_row_dend = FALSE, #row_dend_side = "right",
                show_column_dend = FALSE, column_title = "pairwise correlation between cells",
                heatmap_legend_param = list(title = "Correlation"), cluster_rows=F, cluster_columns=F,
                top=ha2)
    pdf(paste0(label, '_heatmap.pdf'), width=16, height=9)
    ht_list = draw(ht_list, main_heatmap = "cor")
    dev.off()
}

cell_markers = readRDS(file='../data/cancersea/cellmarkers.rds')
cancermarker = readRDS("../data/cancersea/cancersea.rds")

pemt_geneset = unique(readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene")$gene)
gene.modules <- Sys.glob('lisa_script_modules/*symbols')
gene.list <- lapply(gene.modules, scan, what='')
names(gene.list) <- basename(gsub('.symbols', '', gene.modules))
names(gene.list) <- c("EMT", "G1S", "G2M", "Histone", "Hypoxia", "INTERFERON", "MUSCLE", "TNFA")
gene.list$CAF = c("FAP", "PDPN", "THY1", "MMP2", "MMP11", "PDGFRA", "PDGFRL", "TGFB3", "CTGF")
gene.list$MyoFib = c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA")
gene.list$pEMT <- pemt_geneset
list2df <- function(x) {
    ylist <- list()
    for (y in names(x)) {
        ylist[[y]] = data.frame(module=y, gene=x[[y]])
    }
    do.call('rbind', ylist)
}
gene.list <- list2df(gene.list)

library(glue)
library(tidyverse)
library(vroom)
library(clusterProfiler)

for (label in c('20082', '21202', '20696', '29806')) {
    df = read.table(Sys.glob(glue('../results/seurat_sara/{label}*_tumoronly_res0.8.xls')))
    cluster_cancer = list()
    cluster_internal_cancer = list()
    cluster_normal = list()
    for (cluster in unique(df$cluster)) {
        y1 <- enricher(as.character(df[df$cluster==cluster, ]$genename),
                       TERM2GENE=cell_markers, minGSSize=1)
        print('aaaa')
        ## if (nrow(y1) > 0) {
        if (!is.null(y1)  && nrow(y1)!=0) {
            cluster_normal[[as.character(cluster)]] <- cbind(cluster, head(y1[order(y1$p.adjust), ], 5))
        }
        y2 <- enricher(as.character(df[df$cluster == cluster, ]$genename),
                       TERM2GENE = cancermarker, minGSSize = 1)
        ## if (nrow(y2) > 0) {
        if (!is.null(y2) && nrow(y2)!=0) {
            cluster_cancer[[as.character(cluster)]] <- cbind(cluster, head(y2[order(y2$p.adjust), ], 5))
        }
        y3 <- enricher(as.character(df[df$cluster == cluster, ]$genename),
                       TERM2GENE = gene.list, minGSSize = 1)
        ## if (nrow(y3) > 0) {
        if (!is.null(y3)  && nrow(y3)!=0) {
            cluster_internal_cancer[[as.character(cluster)]] <- cbind(cluster, head(y3[order(y3$p.adjust), ], 5))
        }
    }
    cluster_cancer <- do.call(rbind, cluster_cancer)
    cluster_internal_cancer <- do.call(rbind, cluster_internal_cancer)
    cluster_normal <- do.call(rbind, cluster_normal)
    normal.clusters <- subset(cluster_normal[,c(1, 7)], p.adjust<=0.05)
    cancer.clusters <- rbind(cbind(subset(cluster_cancer[,c(1,2,7)], p.adjust<=0.05), Resource='CancerSEA'),
                             cbind(subset(cluster_internal_cancer[,c(1,2,7)], p.adjust<=0.05), Resource='RMS'))
    cancer.clusters <- cancer.clusters[order(cancer.clusters$cluster), ]
    rownames(cancer.clusters) <- paste0(cancer.clusters$cluster, '.', cancer.clusters$ID, '.', cancer.clusters$Resource)
    cancer.clusters <- cancer.clusters[, c(1, 3)]
    rownames(normal.clusters) = paste0(rownames(normal.clusters), '.CellMarker')
    allann <- rbind(cancer.clusters, normal.clusters)
    allann <- allann[order(allann$cluster), ]
    allann$ID <- rownames(allann)
    readr::write_csv(allann, path=glue('{label}_tumor_only_annotation.csv'))
}
