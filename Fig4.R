library(Seurat)
library(glue)
library(ComplexHeatmap)
library(GetoptLong)

source('DEGs_seurat3_sara.R')

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

annotation <- read.delim('../final_annotations/primary_clusters.txt', sep='\t',
                         row.names=1)

primary <- lapply(c('../results/seurat_sara/20082_hg19_premrna_seurat-object.rds',
                    '../results/seurat_sara/20696_seurat-object.rds',
                    '../results/seurat_sara/21202_hg19_premrna_seurat-object.rds',
                    '../results/seurat_sara/29806_hg19_premrna_seurat-object.rds'), readRDS)
primary.tumors <- lapply(c('20082_recluster2_tumor_only.rds', # '../figures/20082_hg19_premrna_tumoronly_res0.8_umap.rds',
                           '../figures/20696_hg19_tumoronly_res0.8_umap.rds',
                           '../figures/21202_hg19_premrna_tumoronly_res0.8_umap.rds',
                           '../figures/29806_hg19_premrna_tumoronly_res0.8_umap.rds'), readRDS)

labels <- c("20082", "20696", "21202", "29806")
erms.topsign = c("ITM2A", "FABP5", "EMILIN1", "RRBP1", "CCND1", "EMP3", "EIF4EBP1", "TSTA3", "ADA", "HMGA2")
arms.topsign = c("HSPB2", "MYL4", "PIPOX", "TNNT2", "MYOG", "ENO3", "NRN1", "GYPC", "TSPAN3", "TFF3")

cat(erms.topsign, sep='\n', file='../final_annotations/gene_modules/ERMS_core.txt')
cat(arms.topsign, sep='\n', file='../final_annotations/gene_modules/ARMS_core.txt')

library(colorBrewer)
library(ggplot2)
library(circlize)

props = list()
for (i in 1:length(labels)) {
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
    label = labels[i]
    primary_case = primary[[i]]
    tumor = primary.tumors[[i]]
    primary_case$seurat_clusters = rep('Non-tumor', length(primary_case$seurat_clusters))
    primary_case$seurat_clusters[names(primary_case$seurat_clusters) %in% names(Idents(tumor))] = 'Tumor'
    primary_case$seurat_clusters = as.factor(primary_case$seurat_clusters)
    states = unlist(as.vector(annotation[label, ]))
    states = as.vector(states[states!=''])
    tumor$RNA_snn_res.0.8 = tumor$seurat_clusters
    levels(tumor$seurat_clusters) = states
    props[[label]] = table(tumor$seurat_clusters)
    pdf(glue("Fig4A_{label}.pdf"), width=12.5, height=4)
    p1 = DimPlot(primary_case, group.by='seurat_clusters')
    p2 = DimPlot(tumor, group.by='seurat_clusters', cols=metacolors)
    print(CombinePlots(plots=list(p1, p2), ncol=2))
    dev.off()
    pdf(glue("Fig4A_{label}_label.pdf"), width=12.5, height=4)
    p1 = DimPlot(primary_case, group.by='RNA_snn_res.0.8')
    p2 = DimPlot(tumor, group.by='RNA_snn_res.0.8')
    print(CombinePlots(plots=list(p1, p2), ncol=2))
    dev.off()
    pdf(glue("Supp2_{label}_dot.pdf"), width=19, height=4)
    p1 <- DotPlot(tumor, features=erms.topsign,
                  cols = c("lightgrey", "red"), group.by='seurat_clusters')
    p2 <- DotPlot(tumor, features=arms.topsign,
                  cols = c("lightgrey", "red"), group.by='seurat_clusters')
    print(CombinePlots(plots=list(p1, p2)), ncol=2)
    dev.off()
    pdf(glue("Supp2_{label}_marker.pdf"), width=8.5, height=5)
    p2 <- FeaturePlot(tumor, features=c("MYOD1", "MYF5", "MYOG", "MYLPF"))
    print(p2)
    dev.off()
    tumor$seurat_clusters = reorder(tumor$seurat_clusters,
                                    new.order=names(sort(table(tumor$seurat_clusters))))
    gene.modules <- Sys.glob('../final_annotations/gene_modules/*txt')
    gene.list <- lapply(gene.modules, scan, what='')
    names(gene.list) <- basename(gsub('.txt', '', gene.modules))
    list2df <- function(x) {
        ylist <- list()
        for (y in names(x)) {
            ylist[[y]] = data.frame(module=y, gene=x[[y]])
        }
        do.call('rbind', ylist)
    }
    gene.list = list2df(gene.list)
    cores = gene.list[grepl("core", gene.list[,1]), ]
    gene.list <- gene.list[as.vector(gene.list[,1])%in%c(levels(tumor$seurat_clusters), 'ARMS_core', 'ERMS_core'), ]
    geneord = names(sort(table(tumor$seurat_clusters)))
    gene.list[,1] = reorder(droplevels(gene.list[,1]),
                            new.order=c(geneord[geneord!='Ground'],
                                        'ARMS_core', 'ERMS_core'))
    head(gene.list)
    cpm = as.data.frame(apply(as.matrix(tumor$RNA@counts), 2, correct))
    sortgenes = order(gene.list[, 1])
    gene.list = gene.list[sortgenes,]
    sortcells = order(tumor$seurat_clusters, tumor$RNA_snn_res.0.8)
    heatdata  <- cpm[as.character(gene.list[,2]), sortcells]
    clusters <- tumor$RNA_snn_res.0.8[sortcells]
    annrow <- gene.list[, 1, drop=F]
    rownames(annrow) <- paste0(gene.list[, 1], '.', gene.list[, 2])
    ## gene.list[grepl("core", gene.list[,1]), 2] = paste0(gene.list[grepl("core", gene.list[,1]), 2], '.1')
    ## rownames(annrow) <- gene.list[, 2]
    anncol <- data.frame(cluster=tumor$seurat_clusters[sortcells])
    ## mat2 = t(apply(heatdata, 1, function(x) {
    ##     q10 <- quantile(x, 0.1)
    ##     q90 <- quantile(x, 0.9)
    ##     x[x < q10] <- q10
    ##     x[x > q90] <- q90
    ##     scale(x, scale=T)
    ## }))
    ## selection = complete.cases(mat2)
    ## mat2 = mat2[selection, ]
    ## annrow = annrow[selection, ,drop=F]
    ## anncol = anncol[selection, ,drop=F]
    mat2= t(scale(t(heatdata)))
    (a<-levels(clusters))
    test=cbind(as.vector(tumor$RNA_snn_res.0.8)[sortcells], as.vector(tumor$seurat_clusters)[sortcells])
    ## test=table(test[,2], test[,1])
    test=table(test[,1], test[,2])
    test=t(test[,names(sort(table(tumor$seurat_clusters)))])
    cluster_cols=c()
    nn=c()
    for (i in seq_along(rownames(test))) {
        if (length(cluster_cols) > 0 && cluster_cols[length(cluster_cols)]=='gray')
            cluster_cols = c(cluster_cols, rep(c('white', 'gray'), sum(test[i, ]>0))[seq(1, sum(test[i, ]>0))])
        else
            cluster_cols = c(cluster_cols, rep(c('gray', 'white'), sum(test[i, ]>0))[seq(1, sum(test[i, ]>0))])
        nn = c(nn, colnames(test)[test[i, ]>0])
    }
    names(cluster_cols)= nn
    metacolors = c(metacolors, 'ARMS_core'='coral', 'ERMS_core'='chartreuse')
    tiff(glue('Fig4C_{label}_heatmap.tiff'), units="in", width=12, height=5, res=320)
    topha = HeatmapAnnotation(states=as.vector(anncol[,1]),
                              clusters=as.vector(clusters),
                              col=list(states=metacolors,
                                       clusters=cluster_cols),
                              show_legend=T,
                              gp = gpar(col = NA),
                              border = c(states=F, clusters=T))
    leftha = rowAnnotation(modules=annrow[, 1],
                           col=list(modules=metacolors),
                           gp = gpar(col = NA))
    ha = Heatmap(mat2, name = "Scaled Expression",
                 use_raster = TRUE, raster_quality = 2,
                 show_column_dend = FALSE, cluster_rows=F, cluster_columns=F,
                 show_row_names=F,
                 show_column_names = FALSE, #width = unit(8, "cm"),
                 column_title = qq(glue("{label} relative expression for @{ncol(heatdata)} cells")),
                 col = colorRamp2(c(-1, -0.6, -0.4, 0, 1),
                                  c("blue", rgb(97, 233, 234, maxColorValue = 255),
                                    "white",
                                    ## rgb(97, 233, 234, maxColorValue = 255),
                                    "white", "red")),
                 top_annotation=topha,
                 left_annotation=leftha)
    draw(ha)
    dev.off()
}

cols = unique(unlist(lapply(props, names)))
cell.prop = matrix(0, ncol=length(cols), nrow=length(labels))
colnames(cell.prop) = cols
rownames(cell.prop) = labels
for (i in seq_along(labels)) {
    cell.prop[labels[i], names(props[[labels[i]]])] = props[[labels[i]]]
}
write.table(cell.prop, file='primary_tumor_all_cells.xls', quote=F, sep='\t', col.names=NA)

annotation <- read.delim('../final_annotations/fish_clusters.txt', sep='\t', check.names=F,
                         row.names=1)

primary <- lapply(c('../results/seurat/Tumor21_unfilter_seurat_obj_tumors.rds',
                    '../results/seurat/Tumor22_unfilter_seurat_obj_tumors.rds',
                    '../results/seurat/Tumor24_unfilter_seurat_obj_tumors.rds'), readRDS)
primary.tumors <- lapply(c('../results/seurat_v6/Tumor21_recluster1.8.rds',
                           '../results/seurat_v6/Tumor22_recluster1.8.rds',
                           '../results/seurat_intersect_velocity/Tumor24_seu.rds'), readRDS)
labels <- c("Tumor21", "Tumor22", "Tumor24")

human_ortholog = read.table('~/langenau/01_rms_projects/01_fish/data/ortholog_mapping/Beagle_fish_human_all_genes.txt', header=T, sep='\t', stringsAsFactors=F)

props = list()
for (i in 1:length(labels)) {
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
    label = labels[i]
    primary_case = primary[[i]]
    tumor = primary.tumors[[i]]
    primary_case$seurat_clusters = rep('Non-tumor', length(primary_case$seurat_clusters))
    primary_case$seurat_clusters[names(primary_case$seurat_clusters) %in% names(Idents(tumor))] = 'Tumor'
    primary_case$seurat_clusters = as.factor(primary_case$seurat_clusters)
    states = unlist(as.vector(annotation[label, ]))
    states = as.vector(states[states!=''])
    levels(tumor$seurat_clusters) = states
    props[[label]] = table(tumor$seurat_clusters)
    pdf(glue("Fig4A_{label}.pdf"), width=12.5, height=4)
    p1 = DimPlot(primary_case, group.by='seurat_clusters')
    p2 = DimPlot(tumor, group.by='seurat_clusters', cols=metacolors)
    print(CombinePlots(plots=list(p1, p2), ncol=2))
    dev.off()
    pdf(glue("Fig4A_{label}_label.pdf"), width=12.5, height=4)
    p1 = DimPlot(primary_case, group.by='SCT_snn_res.0.8')
    if (label != 'Tumor24')
        p2 = DimPlot(tumor, group.by='SCT_snn_res.1.8')
    else
        p2 = DimPlot(tumor, group.by='SCT_snn_res.0.8')
    print(CombinePlots(plots=list(p1, p2), ncol=2))
    dev.off()
    pdf(glue("Supp1_{label}_dot.pdf"), width=6.5, height=4)
    p1 <- DotPlot(tumor, features=na.omit(human_ortholog[match(erms.topsign, human_ortholog$Hsortholog), 'Gene']),
                  cols = c("lightgrey", "red"), group.by='seurat_clusters')
    print(p1)
    dev.off()
    pdf(glue("Supp1_{label}_marker.pdf"), width=8.5, height=5)
    p2 <- FeaturePlot(tumor, features=na.omit(human_ortholog[match(c("MYOD1", "MYF5", "MYOG", "MYLPF"), human_ortholog$Hsortholog), 'Gene']))
    print(p2)
    dev.off()
    tumor$seurat_clusters = reorder(tumor$seurat_clusters,
                                    new.order=names(sort(table(tumor$seurat_clusters))))
    gene.modules <- Sys.glob('../final_annotations/gene_modules/*txt')
    gene.list <- lapply(gene.modules, scan, what='')
    names(gene.list) <- basename(gsub('.txt', '', gene.modules))
    list2df <- function(x) {
        ylist <- list()
        for (y in names(x)) {
            ylist[[y]] = data.frame(module=y, gene=x[[y]])
        }
        do.call('rbind', ylist)
    }
    gene.list = list2df(gene.list)
    cores = gene.list[grepl("core", gene.list[,1]), ]
    gene.list <- gene.list[as.vector(gene.list[,1])%in%c(levels(tumor$seurat_clusters), 'ARMS_core', 'ERMS_core'), ]
    ## gene.list <- gene.list[as.vector(gene.list[,1])%in%levels(tumor$seurat_clusters), ]
    gene.list$gene <- human_ortholog[match(gene.list[, 2], human_ortholog[,3]), 2]
    gene.list <- gene.list[complete.cases(gene.list), ]
    geneord = names(sort(table(tumor$seurat_clusters)))
    gene.list[,1] = reorder(droplevels(gene.list[,1]),
                            new.order=c(geneord[geneord!='Ground'],
                                        'ARMS_core', 'ERMS_core'))
    cpm = as.data.frame(apply(as.matrix(tumor$RNA@counts), 2, correct))
    sortgenes = order(gene.list[, 1])
    gene.list = gene.list[sortgenes,]
    if (label != 'Tumor24') {
        sortcells = order(tumor$seurat_clusters, tumor$SCT_snn_res.1.8)
        clusters <- tumor$SCT_snn_res.1.8[sortcells]
    } else {
        sortcells = order(tumor$seurat_clusters, tumor$SCT_snn_res.0.8)
        clusters <- tumor$SCT_snn_res.0.8[sortcells]
    }
    heatdata  <- cpm[as.character(gene.list[,2]), sortcells]
    annrow <- gene.list[, 1, drop=F]
    rownames(annrow) <- paste0(gene.list[, 1], '.', gene.list[, 2])
    ## rownames(annrow) <- gene.list[, 2]
    anncol <- data.frame(cluster=tumor$seurat_clusters[sortcells])
    ## mat2 = t(apply(heatdata, 1, function(x) {
    ##     q10 <- quantile(x, 0.1)
    ##     q90 <- quantile(x, 0.9)
    ##     x[x < q10] <- q10
    ##     x[x > q90] <- q90
    ##     scale(x, scale=T)
    ## }))
    ## mat2 = mat2[selection, ]
    mat2= t(scale(t(heatdata)))
    selection = complete.cases(mat2)
    annrow = annrow[selection, ,drop=F]
    mat2=mat2[selection, ]
    (a<-levels(clusters))
    if (label != 'Tumor24') {
        test=cbind(as.vector(tumor$SCT_snn_res.1.8)[sortcells], as.vector(tumor$seurat_clusters)[sortcells])
    } else {
        test=cbind(as.vector(tumor$SCT_snn_res.0.8)[sortcells], as.vector(tumor$seurat_clusters)[sortcells])
    }
    ## test=table(test[,2], test[,1])
    test=table(test[,1], test[,2])
    test=t(test[,names(sort(table(tumor$seurat_clusters)))])
    cluster_cols=c()
    nn=c()
    for (i in seq_along(rownames(test))) {
        if (length(cluster_cols) > 0 && cluster_cols[length(cluster_cols)]=='gray')
            cluster_cols = c(cluster_cols, rep(c('white', 'gray'), sum(test[i, ]>0))[seq(1, sum(test[i, ]>0))])
        else
            cluster_cols = c(cluster_cols, rep(c('gray', 'white'), sum(test[i, ]>0))[seq(1, sum(test[i, ]>0))])
        nn = c(nn, colnames(test)[test[i, ]>0])
    }
    names(cluster_cols)= nn
    metacolors = c(metacolors, 'ARMS_core'='coral', 'ERMS_core'='chartreuse')
    tiff(glue('Fig4C_{label}_heatmap.tiff'), units="in", width=12, height=5, res=300)
    topha = HeatmapAnnotation(states=as.vector(anncol[,1]),
                              clusters=as.vector(clusters),
                              col=list(states=metacolors,
                                       clusters=cluster_cols),
                              show_legend=T,
                              gp = gpar(col = NA),
                              border = c(states=F, clusters=T))
    leftha = rowAnnotation(modules=annrow[, 1],
                           col=list(modules=metacolors),
                           gp = gpar(col = NA))
    ha = Heatmap(mat2, name = "Scaled Expression",
                 use_raster = TRUE, raster_quality = 2,
                 show_column_dend = FALSE, cluster_rows=F, cluster_columns=F,
                 show_row_names=F,
                 show_column_names = FALSE, #width = unit(8, "cm"),
                 column_title = qq(glue("{label} relative expression for @{ncol(heatdata)} cells")),
                 col = colorRamp2(c(-1, -0.6, -0.4, 0, 1),
                                  c("blue", rgb(97, 233, 234, maxColorValue = 255),
                                    "white",
                                    ## rgb(97, 233, 234, maxColorValue = 255),
                                    "white", "red")),
                 top_annotation=topha,
                 left_annotation=leftha)
    draw(ha)
    dev.off()
}

cols = unique(unlist(lapply(props, names)))
cell.prop = matrix(0, ncol=length(cols), nrow=length(labels))
colnames(cell.prop) = cols
rownames(cell.prop) = labels
for (i in seq_along(labels)) {
    cell.prop[labels[i], names(props[[labels[i]]])] = props[[labels[i]]]
}
write.table(cell.prop, file='fish_tumor_all_cells.xls', quote=F, sep='\t', col.names=NA)
