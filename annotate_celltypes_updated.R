library(Seurat)
library(tidyverse)
library(vroom)
library(clusterProfiler)
library(glue)
library(DOSE)

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

args$seurat <- '../../01_sc_rms/results/seurat_intersect_velocity/Tumor24_seu.rds'
args$species <- 'fish'
args$label <- 'Tumor24'

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
CombinePlots(plots=list(p1, p2))
dev.off()

pdf(glue('{args$label}_arms_markers.pdf'), width=30, height=10)
if (args$species != 'fish') {
    p1 <- DotPlot(primary1_obj, features=arms.topsign, cols = c("lightgrey", "red"), group.by='RNA_snn_res.0.8')
} else {
    p1 <- DotPlot(primary1_obj, features=arms.topsign, cols = c("lightgrey", "red"), group.by='seurat_clusters')
}
p2 <- FeaturePlot(primary1_obj, features=arms.topsign)
CombinePlots(plots=list(p1, p2))
dev.off()

pdf(glue('{args$label}_myodmyogdesmycpax_markers.pdf'), width=16, height=10)
print(FeaturePlot(primary1_obj, features=human.features))
dev.off()

ligand_target_matrix <- readRDS('../../01_sc_rms/results/nichenet/ligand_target_matrix.rds')
pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)]

gene.modules <- Sys.glob('../../01_sc_rms/final_annotations/gene_modules/*txt')
gene.list <- lapply(gene.modules, scan, what='')
names(gene.list) <- paste0('RMS.', basename(gsub('.txt', '', gene.modules)))

gene.list$CAF = c("FAP", "PDPN", "THY1", "MMP2", "MMP11", "PDGFRA", "PDGFRL", "TGFB3", "CTGF")
gene.list$MyoFib = c("ACTA2", "MCAM", "MYLK", "MYL9", "IL6", "PDGFA")
gene.list$pEMT <- pemt_geneset

list2df <- function(x) {
    ylist <- list()
    for (y in names(x)) {
        ylist[[y]] = data.frame(module=y, gene=paste(x[[y]], collapse=', '))
    }
    do.call('rbind', ylist)
}

gene.list <- list2df(gene.list)
colnames(gene.list) = c("cellMarker", "geneSymbol")

## cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
##    tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
##    dplyr::select(cellMarker, geneID) %>%
##    dplyr::mutate(geneID = strsplit(geneID, ', '))

cell_markers <- (vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
   tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
   dplyr::select(cellMarker, geneSymbol) %>%
   dplyr::mutate(geneID = lapply(strsplit(geneSymbol, ', '), function(x) {gsub('\\]', '', gsub('\\[', '', x))})))

gene.list = tidyr::as_tibble(gene.list) %>% dplyr::mutate(geneID = strsplit(as.character(geneSymbol), ', '))

cancersea <- lapply(Sys.glob("../../01_sc_rms/data/cancersea/*txt"), read.table, header=T)
names(cancersea) <- gsub('.txt', '', basename(Sys.glob("../../01_sc_rms/data/cancersea/*txt")))

cancermarker = list()
for (i in names(cancersea)) {
    cancermarker[[i]] <- cbind(cancermarker=i, GeneID=paste(as.character(cancersea[[i]]$GeneName), collapse=', '))
}
cancermarker = do.call(rbind, cancermarker)

colnames(cancermarker) = c("cellMarker", "geneSymbol")

cancermarker = tidyr::as_tibble(cancermarker) %>% dplyr::mutate(geneID = strsplit(as.character(geneSymbol), ', '))

allgenes = bind_rows(list(gene.list, cancermarker, cell_markers))

if (args$species == 'human') {
    #allmarkers <- read.table(Sys.glob(glue('../results/seurat_sara/{args$label}*xls')))
    allmarkers <- read.table(Sys.glob(glue('../results/seurat_sara/{args$label}*SCT*xls')))
} else {
    allmarkers <- read.table(Sys.glob(glue('{args$label}*SCT*xls')), sep='\t', header=T)
    allmarkers <- cbind(allmarkers, human_ortholog[match(allmarkers$gene, human_ortholog$Gene), ])
    write.table(allmarkers, glue('{args$label}_humanortholog.xls'), quote=F, sep='\t', row.names=F)
}

## filter by cutoffs
allmarkers <- subset(allmarkers, p_val_adj <= 0.01 & (pct.1 >= 0.1 | pct.2 >= 0.1) & (enrichment >= 0.1))

clusters = list()
for (i in unique(allmarkers$cluster)) {
    genes = as.character(allmarkers[allmarkers$cluster==i, ]$gene)
    genes = human_ortholog[human_ortholog$Gene %in% genes, 'Hsortholog']
    genes = genes[(!is.na(genes)) & (genes!="")]
    cat(i, "\t", length(genes), "\n")
    cat(genes, sep='\n', file=glue('lisa/{args$label}_cluster{i}.symbols'))
    y1 = enricher(genes, TERM2GENE=allgenes[,-2], maxGSSize=1000, minGSSize = 5)
    y1 = y1@result
    rownames(y1) = paste0(i, ".", rownames(y1))
    y1$clusters = i
    clusters[[as.character(i)]] <- y1
}

clusters = do.call('rbind', clusters)

clusters = clusters %>% mutate(FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
readr::write_csv(clusters, path=glue('{args$label}_all_annotation.csv'))

clusters = subset(clusters, p.adjust <= 0.01 & qvalue <= 0.01 & Count >= 3 & FoldEnrichment >= 1.2)
readr::write_csv(clusters, path=glue('{args$label}_filter_annotation.csv'))
