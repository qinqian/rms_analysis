## library(nichenetr)
library(Seurat)
library(cellassign)
library(tidyverse)
library(vroom)
library(clusterProfiler)
library(glue)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--seuratobj', dest='seurat', metavar='N', type="character", nargs='+')
    parser$add_argument('--label', dest='label', default='')
    parser$add_argument('--species', dest='species', type='character', default='human')
    args = parser$parse_args()
    args
}

args = get_args()

if (length(args$seurat) == 0 || args$label == '') {
    cat('empty argument, exit..')
    q()
}

args$seurat <- '../results/seurat_sara/20696_seurat-object.rds'
args$label <- '20696'

## args$seurat <- '../results/seurat_sara/21202_hg19_premrna_seurat-object.rds'
## args$label <- '21202'

#args$seurat <- '../results/seurat_sara/20082_hg19_premrna_seurat-object.rds'
#args$label <- '20082'

## args$seurat <- '../results/seurat_sara/29806_hg19_premrna_seurat-object.rds'
## args$label <- '29806'
## print(args)

## args$seurat <- '../results/seurat_sara/MAST139_1cells_seurat-object.rds'
## args$label <- 'MAST139_1cells'
## print(args)

#args$seurat <- '../results/seurat_intersect_velocity/Tumor24_seu.rds'
#args$label <- 'Tumor24'
#args$species <- 'fish'
#print(args)

## args$seurat <- '../figures/20082_hg19_premrna_tumoronly_res0.8_umap.rds'
## args$label <- '20082'
## args$species <- 'human'
## print(args)

primary1_obj <- readRDS(args$seurat)

erms.topsign = c("ITM2A", "FABP5", "EMILIN1", "RRBP1", "CCND1", "EMP3", "EIF4EBP1", "TSTA3", "ADA", "HMGA2")
arms.topsign = c("HSPB2", "MYL4", "PIPOX", "TNNT2", "MYOG", "ENO3", "NRN1", "GYPC", "TSPAN3", "TFF3")
human.features = c("MYOD1", "MYF5", "MYOG", "DES", "MYC", "PAX3", "PAX7", "KRAS",
                   "NRAS", "MYLPF",
                   "CDK4", "CDK6", "WEE1")
if (args$species != 'human') {
    human_ortholog = read.table('~/langenau//01_rms_projects//01_fish/data/ortholog_mapping/Beagle_fish_human_all_genes.txt', header=T, sep='\t', stringsAsFactors=F)
    erms.topsign = human_ortholog[human_ortholog$Hsortholog %in% erms.topsign, 'Gene']
    arms.topsign = human_ortholog[human_ortholog$Hsortholog %in% arms.topsign, 'Gene']
    human.features = human_ortholog[human_ortholog$Hsortholog %in% human.features, 'Gene']
}

primary1_obj.copy <- readRDS(args$seurat)

pdf(glue('{args$label}_with_clusterlabel.pdf'), width=6, height=4)
if (args$species != 'fish') {
    primary1_obj.copy$seurat_clusters = primary1_obj.copy$RNA_snn_res.0.8
}
DimPlot(primary1_obj.copy, reduction='umap', group.by='seurat_clusters', label = TRUE) #, legend = "none")
dev.off()

pdf(glue('{args$label}_erms_markers.pdf'), width=30, height=10)
if (args$species != 'fish') {
    p1 <- DotPlot(primary1_obj, features=erms.topsign, cols = c("lightgrey", "red"), group.by='RNA_snn_res.0.8')
} else {
    p1 <- DotPlot(primary1_obj, features=erms.topsign, cols = c("lightgrey", "red"), group.by='seurat_clusters')
}
p2 <- FeaturePlot(primary1_obj, features=erms.topsign)
## p3 <- DoHeatmap(primary1_obj, features=erms.topsign, slots='raw.data')
CombinePlots(plots=list(p1, p2))
dev.off()

pdf(glue('{args$label}_arms_markers.pdf'), width=30, height=10)
if (args$species != 'fish') {
    p1 <- DotPlot(primary1_obj, features=arms.topsign, cols = c("lightgrey", "red"), group.by='RNA_snn_res.0.8')
} else {
    p1 <- DotPlot(primary1_obj, features=arms.topsign, cols = c("lightgrey", "red"), group.by='seurat_clusters')
}
p2 <- FeaturePlot(primary1_obj, features=arms.topsign)
## p3 <- DoHeatmap(primary1_obj, features=erms.topsign, slots='raw.data')
CombinePlots(plots=list(p1, p2))
dev.off()

pdf(glue('{args$label}_myodmyogdesmycpax_markers.pdf'), width=16, height=10)
print(FeaturePlot(primary1_obj, features=human.features))
dev.off()

ligand_target_matrix <- readRDS('../results/nichenet/ligand_target_matrix.rds')
pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.

## gene.modules <- Sys.glob('lisa_script_modules/*symbols')
gene.modules <- Sys.glob('../final_annotations/gene_modules/*txt')
gene.list <- lapply(gene.modules, scan, what='')
names(gene.list) <- basename(gsub('.txt', '', gene.modules))

gene.list$CAF = c("FAP", "PDPN", "THY1", "MMP2", "MMP11", "PDGFRA", "PDGFRL", "TGFB3", "CTGF")
gene.list$MyoFib = c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA")
gene.list$pEMT <- pemt_geneset

list2df <- function(x) {
    ylist <- list()
    for (y in names(x)) {
        ylist[[y]] = data.frame(module=y, gene=x[[y]])
    }
    do.call('rbind', ylist)
    ## ylist
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

print(glue('../results/seurat_sara/{args$label}*xls'))
if (args$species == 'human') {
    #allmarkers <- read.table(Sys.glob(glue('../results/seurat_sara/{args$label}*xls')))
    allmarkers <- read.table(Sys.glob(glue('../results/seurat_sara/{args$label}*SCT*xls')))
} else {
    allmarkers <- read.table(Sys.glob(glue('../results/seurat_v6/{args$label}*SCT*xls')), sep='\t', header=T)
    allmarkers <- cbind(allmarkers, human_ortholog[match(allmarkers$gene, human_ortholog$Gene), ])
    write.table(allmarkers, glue('{args$label}_humanortholog.xls'), quote=F, sep='\t', row.names=F)
}

cluster_cancer = list()
cluster_internal_cancer = list()
cluster_normal = list()

for (cluster in unique(allmarkers$cluster)) {
    if (args$species == 'human') {
        genes = as.character(allmarkers[allmarkers$cluster==cluster, ]$genename)
    } else {
        genes = as.character(allmarkers[allmarkers$cluster==cluster, ]$gene)
        genes = human_ortholog[human_ortholog$Gene %in% genes, 'Hsortholog']
    }
    y1 <- enricher(genes,
                   TERM2GENE=cell_markers, minGSSize=1)
    ## if (nrow(y1) > 0) {
    if (!is.null(y1)  && nrow(y1)!=0) {
        cluster_normal[[as.character(cluster)]] <- cbind(cluster, head(y1[order(y1$p.adjust), ], 5))
    }
    y2 <- enricher(genes,
                   TERM2GENE = cancermarker, minGSSize = 1)
    ## if (nrow(y2) > 0) {
    if (!is.null(y2) && nrow(y2)!=0) {
        cluster_cancer[[as.character(cluster)]] <- cbind(cluster, head(y2[order(y2$p.adjust), ], 5))
    }
    y3 <- enricher(genes,
                  TERM2GENE = gene.list, minGSSize = 1)
    if (!is.null(y3)  && nrow(y3)!=0) {
        cluster_internal_cancer[[as.character(cluster)]] <- cbind(cluster, head(y3[order(y3$p.adjust), ], 5))
    }
}

cluster_cancer <- do.call(rbind, cluster_cancer)
cluster_internal_cancer <- do.call(rbind, cluster_internal_cancer)
cluster_normal <- do.call(rbind, cluster_normal)

normal.clusters <- subset(cluster_normal[,c(1, 3, 4, 5, 7)], p.adjust<=0.05)

cancer.clusters <- rbind(cbind(subset(cluster_cancer[,c(1,2,4,5,7)], p.adjust<=0.05), Resource='CancerSEA'),
                         cbind(subset(cluster_internal_cancer[,c(1,2,4,5,7)], p.adjust<=0.05), Resource='RMS'))
cancer.clusters <- cancer.clusters[order(cancer.clusters$cluster), ]

rownames(cancer.clusters) <- paste0(cancer.clusters$cluster, '.', cancer.clusters$ID, '.', cancer.clusters$Resource)
cancer.clusters <- cancer.clusters[, c(1, 3, 4, 5)]

rownames(normal.clusters) = paste0(normal.clusters[,1], '.', normal.clusters[,2], '.CellMarker')
normal.clusters <- normal.clusters[, c(1, 3, 4, 5)]

allann <- rbind(cancer.clusters, normal.clusters)
allann <- allann[order(allann$cluster), ]
allann$ID <- rownames(allann)
readr::write_csv(allann, path=glue('{args$label}_all_annotation.csv'))

#quit()
library(data.table)

DT <- data.table(allann)
gsea.ann <- DT[, lapply(.SD, function(x) paste(head(x[order(p.adjust)], 5), collapse='|  ')), by=cluster, .SDcols=3]

data(example_TME_markers)
library(SingleCellExperiment)

primary1.sce = as.SingleCellExperiment(primary1_obj)
sizeFactors(primary1.sce) <- colSums(assay(primary1.sce))

if (args$species == 'fish') {
    example_TME_markers$symbol = lapply(example_TME_markers$symbol, function(x) {
        human_ortholog[human_ortholog$Hsortholog %in% x, 'Gene']
    })
}

marker_mat = marker_list_to_mat(example_TME_markers$symbol)

sce_marker <- primary1.sce[intersect(rownames(marker_mat), rownames(primary1.sce)),]

s = sizeFactors(sce_marker)
sce_marker = sce_marker[, s > 0]
s = s[s > 0]

cas <- cellassign(exprs_obj = sce_marker,
                  marker_gene_info = marker_mat[intersect(rownames(marker_mat), rownames(primary1.sce)),],
                  s = s)

colnames(sce_marker)

primary1.objann = subset(primary1_obj, cells=colnames(sce_marker))
primary1.objann$seurat_clusters = cas$cell_type
pdf(glue('{args$label}_cellassign.pdf'), width=6, height=4)
DimPlot(primary1.objann, reduction='umap', group.by='seurat_clusters')
dev.off()

#if (args$species == 'human') {
#    levels(primary1_obj$RNA_snn_res.0.8)[match(gsea.ann$cluster, levels(primary1_obj$RNA_snn_res.0.8))] <- gsea.ann$ID
#} else {
#    levels(primary1_obj$seurat_clusters)[match(gsea.ann$cluster, levels(primary1_obj$seurat_clusters))] <- gsea.ann$ID
#}
#
#pdf(glue('{args$label}_clusterprofiler.pdf'), width=20, height=9)
#if (args$species == 'human') {
#    DimPlot(primary1_obj, reduction='umap', group.by='RNA_snn_res.0.8')
#} else {
#    DimPlot(primary1_obj, reduction='umap', group.by='seurat_clusters')
#}
#dev.off()
#
