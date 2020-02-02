library(nichenetr)
library(Seurat)
library(tidyverse)
library(vroom)
library(clusterProfiler)

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
gene.list <- list2df(gene.list)

cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
   tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
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

primary1 <- read.table(Sys.glob('../results/seurat_sara/20696*SCT*'))
cluster_cancer = list()
cluster_internal_cancer = list()
cluster_normal = list()
for (cluster in unique(primary1$cluster)) {
    y1 <- enricher(as.character(primary1[primary1$cluster==cluster, ]$genename),
                   TERM2GENE=cell_markers, minGSSize=1)
    if (nrow(y1) > 0) {
        cluster_normal[[as.character(cluster)]] <- cbind(cluster, head(y1[order(y1$p.adjust), ], 5))
    }
    y2 <- enricher(as.character(primary1[primary1$cluster == cluster, ]$genename),
                  TERM2GENE = cancermarker, minGSSize = 1)
    if (nrow(y2) > 0) {
        cluster_cancer[[as.character(cluster)]] <- cbind(cluster, head(y2[order(y2$p.adjust), ], 5))
    }
    y3 <- enricher(as.character(primary1[primary1$cluster == cluster, ]$genename),
                  TERM2GENE = gene.list, minGSSize = 1)
    if (nrow(y3) > 0) {
        cluster_internal_cancer[[as.character(cluster)]] <- cbind(cluster, head(y3[order(y3$p.adjust), ], 5))
    }
}

cluster_cancer <- do.call(rbind, cluster_cancer)
cluster_internal_cancer <- do.call(rbind, cluster_internal_cancer)
cluster_normal <- do.call(rbind, cluster_normal)

## subset(cluster_normal[,c(1,2,7)], p.adjust<=0.05)
cancer.clusters <- rbind(cbind(subset(cluster_cancer[,c(1,2,7)], p.adjust<=0.05), Resource='CellMarker'),
                         cbind(subset(cluster_internal_cancer[,c(1,2,7)], p.adjust<=0.05), Resource='RMS'))
cancer.clusters <- cancer.clusters[order(cancer.clusters$cluster), ]

primary1_obj = readRDS('../results/seurat_sara/20696_seurat-object.rds')

DimPlot(primary1_obj, reduction='umap')
dev.off()

## hnscc_expression <- readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
## saveRDS(hnscc_expression, "hnscc_expression.rds")
## expression <- hnscc_expression$expression
## sample_info <- hnscc_expression$sample_info
## ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
## saveRDS(ligand_target_matrix, "ligand_target_matrix.rds")
## weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
## saveRDS(weighted_networks, "weighted_networks.rds")

## weighted_networks <- readRDS('weighted_networks.rds')
## ligand_target_matrix <- readRDS('ligand_target_matrix.rds')

## weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network %>%
##                                                                 distinct(from, to), by = c("from", "to"))

## mast111 <- readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds')

## colortab <- read.table('color_table.xls', sep="\t", header=T)
## col <- colortab[,'MAST111']
## col <- col[col!=""]

## mast111@meta.data$celltype <- mast111@meta.data$RNA_snn_res.0.8
## levels(mast111@meta.data$celltype) <- col

## mast111@meta.data %>% head()
## mast111@meta.data$celltype %>% table()
