## library(nichenetr)
library(Seurat)
library(cellassign)
library(tidyverse)
library(vroom)
library(clusterProfiler)

primary1_obj <- readRDS('../results/seurat_sara/20696_seurat-object.rds')
Idents(primary1_obj) <- primary1_obj$RNA_snn_res.0.8
pdf('20696_label_on_umap.pdf')
DimPlot(primary1_obj, label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

primary1_obj <- readRDS('../results/seurat_sara/20082_hg19_premrna_seurat-object.rds')
Idents(primary1_obj) <- primary1_obj$RNA_snn_res.0.8
pdf('20082_label_on_umap.pdf')
DimPlot(primary1_obj, label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

levels(Idents(primary1_obj)) <- c("Tumor", "Tumor", "Tumor", "Monocytes",
                                  "Tumor", "Alvelolar Epithelial Cells (Met site)", "Alvelolar Epithelial Cells (Met site)",
                                  "Tumor", "Tumor", "Smooth Muscle of Lung\nPericytes\nMSCs\nFibroblasts",
                                  "Tumor", "T-cells", "Vasculature\nEndothelial Cells", "Lymphocytes",
                                  "Smooth Muscle of Lung\nPericytes\nMSCs\nFibroblasts",  "Vasculature\nEndothelial Cells")

pdf("20696_daveannotation_on_umap.pdf")
DimPlot(primary1_obj, label.size=3.5, label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

library(ComplexHeatmap)
primary1_tumor <- subset(primary1_obj, ident='Tumor')

cluster_ids <- c("Ground", "Cell cycle", "Muscle", "UNK",
                 "EMT", "UNK", "Ribosome (stress)")
cells <- droplevels(primary1_tumor@meta.data$RNA_snn_res.0.8)
levels(cells) <- cluster_ids

primary1_tumor@meta.data$label <- cells

primary1_tumor <- RunPCA(object = primary1_tumor, npcs=30)
primary1_tumor <- FindNeighbors(object = primary1_tumor, dims = 1:20)
primary1_tumor <- FindClusters(primary1_tumor, resolution=0.8)
primary1_tumor <- RunUMAP(primary1_tumor, dims = 1:20, reduction = "pca")

allcluster <- names((primary1_tumor@meta.data$seurat_clusters) %>% table())

library(foreach)
source('DEGs_seurat3_sara.R')
clusterde <- list()
for (i in allcluster) {
    print(i)
    print(allcluster[-(as.integer(i)+1)])
    de.up <- get_upregulated_genes_in_clusters(primary1_tumor, i, allcluster[-(as.integer(i)+1)])
    de.up$cluster <- i
    clusterde[[i]] <- de.up
}
seurat1.de <- do.call('rbind', clusterde) ## still too many genes

## again, filter by adjusted p value and fraction of cells expressing the genes
seurat1.de <- subset(seurat1.de, p.adjusted <= 0.01 & (fg_fraction >= 0.1 | bg_fraction >= 0.1))

write.table(seurat1.de, file="20696_tumor_only_cluster_differential_genes.xls", quote=F, sep='\t', col.names=NA)

pdf('20696_tumor_recluster_on_umap.pdf', width=16, height=6)
p1 <- DimPlot(primary1_tumor, label.size=3.5, label = TRUE, pt.size = 0.5) + NoLegend()
Idents(primary1_tumor) <- primary1_tumor@meta.data$label
p2 <- DimPlot(primary1_tumor, label.size=3.5, label = TRUE, pt.size = 0.5) + NoLegend()
CombinePlots(plots=list(p1, p2))
dev.off()

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
gene.list = list2df(gene.list)

library(ggplot2)
gene.list[, 2] <- as.vector(gene.list[, 2])
gene.list <- gene.list[gene.list[, 2] %in% rownames(primary1_tumor$RNA@data), ]
sortcells <- sort(primary1_tumor$seurat_clusters)
heatdata  <- as.matrix(primary1_tumor$RNA@data)[gene.list[, 2], order(primary1_tumor$seurat_clusters)]

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

selection <- apply(heatdata, 1, function(x) sum(x>0)/length(x) > 0.1)
heatdata <- heatdata[selection, ,drop=F]

get_correlated_variable_genes = function(mat, n = nrow(mat), cor_cutoff = 0, n_cutoff = 0) {
    ind = order(apply(mat, 1, function(x) {
            q = quantile(x, c(0.1, 0.9))
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
mat = heatdata

mat2 = t(apply(mat, 1, function(x) {
    q10 = quantile(x, 0.05)
    q90 = quantile(x, 0.95)
    x[x < q10] = q10
    x[x > q90] = q90
    scale(x)
}))
colnames(mat2) = colnames(mat)
base_mean = rowMeans(t(mat))

library(GetoptLong)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

ha = structure(brewer.pal(length(unique(anncol[, 1])), "Set3"),
               names=levels(anncol[, 1]))
ha2 = HeatmapAnnotation(cluster=anncol[, 1],
                       col=list(cluster=ha), show_legend=F)

ha4 = HeatmapAnnotation(modules=as.vector(gene.list[selection, 1]),
                        col=list(modules=ann_colors$module), show_legend=T)

ht_list = Heatmap(t(mat2), col = colorRamp2(c(-2.5, 0, 3.5), c("blue", "white", "red")),
    name = "scaled_expr", column_title = qq("relative expression for @{ncol(mat)} cells"),
    show_column_dend = FALSE, cluster_rows=F, cluster_columns=F,
    show_column_names = FALSE, width = unit(8, "cm"), show_row_dend = FALSE,
    heatmap_legend_param = list(title = "Scaled expr"), top=ha4) +
    Heatmap(anncol[,1], col=ha, width = unit(0.7, "cm"), name='clusters') +
    Heatmap(cor(mat2), name = "cor",
        col = colorRamp2(c(-0.3, 0, 0.3), c("blue", "white", "red")),
        show_row_names = FALSE, show_column_names = FALSE, show_row_dend = FALSE, #row_dend_side = "right",
        show_column_dend = FALSE, column_title = "pairwise correlation between cells",
        heatmap_legend_param = list(title = "Correlation"), cluster_rows=F, cluster_columns=F,
        top=ha2)

pdf('20696_heatmap.pdf', width=16, height=9)
ht_list = draw(ht_list, main_heatmap = "cor")
dev.off()


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

cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
   ## tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
   tidyr::unite("cellMarker", cellName, sep=", ") %>%
   dplyr::select(cellMarker, geneSymbol) %>%
   dplyr::mutate(geneID = strsplit(geneSymbol, ', '))
saveRDS(cell_markers, file='../data/cancersea/cellmarkers.rds')

cancersea <- lapply(Sys.glob("../data/cancersea/*txt"), read.table, header=T)
names(cancersea) <- gsub('.txt', '', basename(Sys.glob("../data/cancersea/*txt")))
cancermarker = list()
for (i in names(cancersea)) {
    cancermarker[[i]] <- cbind(cancermarker=i, GeneID=as.character(cancersea[[i]]$GeneName))
}
cancermarker = do.call(rbind, cancermarker)
saveRDS(cancermarker, file="../data/cancersea/cancersea.rds")

allmarkers <- read.table('20696_tumor_only_cluster_differential_genes.xls')

cluster_cancer = list()
cluster_internal_cancer = list()
cluster_normal = list()
for (cluster in unique(allmarkers$cluster)) {
    y1 <- enricher(as.character(allmarkers[allmarkers$cluster==cluster, ]$genename),
                   TERM2GENE=cell_markers, minGSSize=1)
    if (nrow(y1) > 0) {
        cluster_normal[[as.character(cluster)]] <- cbind(cluster, head(y1[order(y1$p.adjust), ], 3))
    }
    y2 <- enricher(as.character(allmarkers[allmarkers$cluster == cluster, ]$genename),
                  TERM2GENE = cancermarker, minGSSize = 1)
    if (nrow(y2) > 0) {
        cluster_cancer[[as.character(cluster)]] <- cbind(cluster, head(y2[order(y2$p.adjust), ], 3))
    }
    y3 <- enricher(as.character(allmarkers[allmarkers$cluster == cluster, ]$genename),
                  TERM2GENE = gene.list, minGSSize = 1)
    if (nrow(y3) > 0) {
        cluster_internal_cancer[[as.character(cluster)]] <- cbind(cluster, head(y3[order(y3$p.adjust), ], 3))
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

readr::write_csv(allann, path='20696_tumor_only_annotation.csv')

Idents(primary1_tumor) <- primary1_tumor$seurat_clusters

ids <- structure(c("GROUND", "G1S", "Muscle", "GROUND", "GROUND", "GROUND", "GROUND",
   "EMT", "G2M"), names=0:8)

primary1_tumor = RenameIdents(primary1_tumor, ids)

pdf('20696_tumor_daveann_on_umap.pdf', width=8, height=6)
p1 <- DimPlot(primary1_tumor, label.size=3.5, label = TRUE, pt.size = 0.5) + NoLegend()
print(p1)
dev.off()

write.table(as.data.frame(Idents(primary1_tumor)), file='20696_primary1_tumor.xls',
            quote=F, sep='\t', col.names=NA)

library(CytoTRACE)
cyto <- CytoTRACE(as.matrix(primary1_tumor@assays$RNA@counts))

pheno <- as.character(Idents(primary1_tumor))
names(pheno) <- colnames(primary1_tumor)

sum(names(cyto$CytoTRACE) == names(pheno))

plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOD1")

cyto_diff = cyto$CytoTRACE

primary1_tumor@meta.data[names(cyto_diff), 'cyto'] = cyto_diff

out_box = data.frame(state=Idents(primary1_tumor), diff=cyto_diff)

library(ggpubr)
ggplot(out_box, aes(state, diff) ) + geom_boxplot() + ggpubr::theme_pubclean()
ggsave('cyto_box.pdf')
