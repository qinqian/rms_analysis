metacolors <- c(rgb(119, 62, 20, maxColorValue = 255),
                rgb(236, 133, 40, maxColorValue = 255),
                rgb(59, 22, 115, maxColorValue = 255),
                rgb(52, 101, 252, maxColorValue = 255),
                rgb(242, 242, 242, maxColorValue = 255),
                rgb(52, 101, 252, maxColorValue = 255),
                rgb(225, 39, 39, maxColorValue = 255),
                rgb(72, 159,  75, maxColorValue = 255),
                rgb(20, 64, 216, maxColorValue = 255),
                rgb(226, 75, 143, maxColorValue = 255),
                rgb(158, 60, 200, maxColorValue = 255),
                rgb(241, 250, 100, maxColorValue = 255),
                rgb(166, 166, 166, maxColorValue = 255),
                rgb(166, 166, 166, maxColorValue = 255))
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'UPR', 'Unique#7', 'Unique#8')
names(metacolors) <- metalabels

library(Seurat)
library(tidyverse)
library(vroom)
library(clusterProfiler)
library(glue)
library(DOSE)
library(glue)
library(ComplexHeatmap)
library(GetoptLong)
library(ggplot2)
library(circlize)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--seuratobj', dest='seurat', metavar='N', type="character", nargs='+')
    parser$add_argument('--label', dest='label', default='')
    parser$add_argument('--de', dest='de', default='')
    parser$add_argument('--species', dest='species', type='character', default='human')
    args = parser$parse_args()
    args
}
args = get_args()
if (length(args$seurat) == 0 || args$label == '') {
    cat('empty argument, exit..')
    q()
}

#args$seurat <- '/PHShome/qq06/langenau/projects/01_sc_rms/phaseA_explore_rms/20082_recluster2_tumor_only.rds'
#args$species <- 'human'
#args$label <- '20082'
#args$de <- '/PHShome/qq06/langenau/projects/01_sc_rms/results/seurat_sara/20082_recluster2_tumoronly_res0.8.xls'

primary1_obj <- readRDS(args$seurat)

get_internal_geneset = function() {
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
    ligand_target_matrix <- readRDS('~/langenau/projects/01_sc_rms/results/nichenet/ligand_target_matrix.rds')
    pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)]

    gene.modules <- Sys.glob('~/langenau/projects/01_sc_rms/final_annotations/gene_modules/*txt')
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
    tidyr::as_tibble(gene.list) %>% dplyr::mutate(geneID = strsplit(as.character(geneSymbol), ', '))
}

get_cell_marker = function() {
    cell_markers <- (vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
                     tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>%
                     dplyr::select(cellMarker, geneSymbol) %>%
                     dplyr::mutate(geneID = lapply(strsplit(geneSymbol, ', '), function(x) {gsub('\\]', '', gsub('\\[', '', x))})))
    cell_markers
}

get_cancer_marker = function() {
    cancersea <- lapply(Sys.glob("~/langenau/projects/01_sc_rms/data/cancersea/*txt"), read.table, header=T)
    names(cancersea) <- gsub('.txt', '', basename(Sys.glob("~/langenau/projects/01_sc_rms/data/cancersea/*txt")))
    cancermarker = list()
    for (i in names(cancersea)) {
        cancermarker[[i]] <- cbind(cancermarker=i, GeneID=paste(as.character(cancersea[[i]]$GeneName), collapse=', '))
    }
    cancermarker = do.call(rbind, cancermarker)
    colnames(cancermarker) = c("cellMarker", "geneSymbol")
    tidyr::as_tibble(cancermarker) %>% dplyr::mutate(geneID = strsplit(as.character(geneSymbol), ', '))
}

get_msigdb = function(category='c2') {
    c2 = read.gmt(Sys.glob(glue("/PHShome/qq06/langenau/01_rms_projects/02_human/data/msigdb/{category}.all.v7.0.symbols.gmt")))
    ylist = list()
    for (i in unique(c2[,1])) {
        ylist[[i]] = data.frame(module=i, gene=paste(subset(c2, ont==i)$gene, collapse=', '))
    }
    ylist = do.call('rbind', ylist)
    colnames(ylist) = c("cellMarker", "geneSymbol")
    ylist = tidyr::as_tibble(ylist) %>% dplyr::mutate(geneID = strsplit(as.character(geneSymbol), ', '))
}

## get gene set
gene.list = get_internal_geneset()
#cell_markers = get_cell_marker()
#cancermarker = get_cancer_marker()
#c2 = get_msigdb()
#c5 = get_msigdb('c5')
#h = get_msigdb('h')
#allgenes = bind_rows(list(gene.list, cancermarker, cell_markers, c2, c5, h))

## all universe genes number
#print(length(na.omit(unique(unlist(allgenes$geneID)))))

print('------------')
print(args$de)
if (args$species == 'human') {
    if (args$de != '') {
        if (args$label == 'MSK72117_res0.9') {
            allmarkers <- read.csv(args$de)
            # resolution 0.9
            primary1_obj = FindClusters(primary1_obj, resolution=0.9)
        } else {
            allmarkers <- read.table(args$de)
            primary1_obj$seurat_clusters = primary1_obj$RNA_snn_res.0.8
        }
    } 
}

library(gdata)
if (args$label == 'C12SC2') {
    primary1_obj$seurat_clusters = reorder(primary1_obj$seurat_clusters, new.order=c('0', '1', '8', '9', '12', '11', '2', '6', '13', '3', '4', '5', '7', '10'))
} else if (args$label == 'MSK72117_res0.9') {
    primary1_obj$seurat_clusters = reorder(primary1_obj$seurat_clusters, new.order=c('0', '1', '2', '4', '3', '6', '7', '8', '9', '10', '11', '5', '12', '13', '14', '15'))
}

write.table(table(primary1_obj$seurat_clusters), file=glue('{args$label}_cells_number.xls'), sep='\t', quote=F)

## filter by cutoffs
##Seurat FindAllMarker
if (args$de != '') {
    ## allmarkers <- subset(allmarkers, p_val_adj <= 0.01 & (pct.1 >= 0.1 | pct.2 >= 0.1) & (enrichment >= 0.1))
    allmarkers <- subset(allmarkers, p.adjusted <= 0.01 & (fg_fraction >= 0.1 | bg_fraction >= 0.1) & (diff_in_enrichment >= 0.1))
} else {
    ##Sowmya 
    allmarkers <- subset(allmarkers, p.adjusted <= 0.01 & (fg_fraction >= 0.1 | bg_fraction >= 0.1) & (diff_in_enrichment >= 0.1))
}

## all overlapped DE genes
# print(length(unique(allmarkers[allmarkers$cluster == 3, 'genename'])))
# print(length(intersect(unique(allmarkers[allmarkers$cluster == 3, 'genename']),
#                        na.omit(unique(unlist(allgenes$geneID))))))

pdf(glue('{args$label}_with_umap.pdf'), width=6, height=4)
p1 = DimPlot(primary1_obj, reduction='umap', group.by='seurat_clusters', label = TRUE)
print(p1)
dev.off()

correct <- function(x) { return (1e6*x/sum(x))}
cpm = as.data.frame(apply(as.matrix(primary1_obj$RNA@counts), 2, correct))
## cpm = primary1_obj$RNA@scale.data # top 2000 variable genes only

plot.genes = list()
for (i in gene.list$cellMarker) {
    plot.genes[[i]] = data.frame(module=i, gene=gene.list[gene.list$cellMarker == i, 'geneID']$geneID[[1]])
}
plot.genes = do.call(rbind, plot.genes)

cpm.sub = cpm[as.character(plot.genes[,2]), ]


sortcells = order(primary1_obj$seurat_clusters)
heatdata  <- cpm.sub[, sortcells]

clusters <- primary1_obj$seurat_clusters[sortcells]
annrow <- plot.genes[, 1, drop=F]
rownames(annrow) <- paste0(plot.genes[, 1], '.', plot.genes[, 2])
anncol <- data.frame(cluster=primary1_obj$seurat_clusters[sortcells])

mat2 = t(apply(heatdata, 1, function(x) {
    q10 <- quantile(x, 0.1, na.rm=T)
    q90 <- quantile(x, 0.9, na.rm=T)
    x[x < q10] <- q10
    x[x > q90] <- q90
    ## x = (x - mean(x)) / sd(x) ^ as.logical(sd(x))
    ## x = log2(x+1)
    scale(x)
}))

selection = complete.cases(mat2)
mat2 = mat2[selection, ]
rownames(mat2) = rownames(annrow)[selection]
colnames(mat2) = rownames(anncol)

pdf(glue("{args$label}_heatmap.pdf"), width=22, height=12)
set.seed(99)
topha = HeatmapAnnotation(states=as.vector(anncol[,1]),
                          show_legend=T,
                          gp = gpar(col = NA),
                          border = c(states=F))
leftha = rowAnnotation(modules=annrow[selection, 1],
                       gp = gpar(col = NA))
ha = Heatmap(mat2, name = "Scaled Expression",
             use_raster = TRUE, raster_quality = 2,
             show_column_dend = FALSE, cluster_rows=F, cluster_columns=F,
             show_row_names=F,
             show_column_names = FALSE,
             column_title = qq(glue("{args$label} relative expression for @{ncol(heatdata)} cells")),
             top_annotation=topha,
             col = colorRamp2(c(min(mat2), -0.6, -0.4, 0, max(mat2)),
                              c("blue",
                                rgb(97, 233, 234, maxColorValue = 255),
                                "white",
                                "white", "red")),
             left_annotation=leftha)
draw(ha)
dev.off()

#clusters = list()
#for (i in unique(allmarkers$cluster)) {
#    genes = as.character(allmarkers[allmarkers$cluster==i, ]$gene)
#    if (args$species != 'human') {
#        genes = human_ortholog[human_ortholog$Gene %in% genes, 'Hsortholog']
#    }
#    genes = genes[(!is.na(genes)) & (genes!="")]
#    cat(i, "\t", length(genes), "\n")
#    ### cat(genes, sep='\n', file=glue('lisa/{args$label}_cluster{i}.symbols'))
#    y1 = enricher(genes, TERM2GENE=allgenes[,-2], maxGSSize=1000, minGSSize = 5)
#    y1 = y1@result
#    rownames(y1) = paste0(i, ".", rownames(y1))
#    y1$clusters = i
#    clusters[[as.character(i)]] <- y1
#}
#
#clusters = do.call('rbind', clusters)
#clusters = clusters %>% mutate(FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
#readr::write_csv(clusters, path=glue('{args$label}_all_annotation.csv'))
#
#clusters = subset(clusters, p.adjust <= 0.01 & qvalue <= 0.01 & Count >= 3 & FoldEnrichment >= 1.2)
#readr::write_csv(clusters, path=glue('{args$label}_filter_annotation.csv'))
