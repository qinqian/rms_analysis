library(Seurat)
library(reticulate)
library(patchwork)
library(ggplot2)
library(GetoptLong)
library(RColorBrewer)
library(circlize)
library(patchwork)
library(plot3D)
library(dplyr)
library(glue)
library(ComplexHeatmap)
set.seed(100)

get_internal_geneset = function() {
    erms.topsign = c("ITM2A", "FABP5", "EMILIN1", "RRBP1", "CCND1", "EMP3", "EIF4EBP1", "TSTA3", "ADA", "HMGA2")
    arms.topsign = c("HSPB2", "MYL4", "PIPOX", "TNNT2", "MYOG", "ENO3", "NRN1", "GYPC", "TSPAN3", "TFF3")
    human.features = c("MYOD1", "MYF5", "MYOG", "DES", "MYC", "PAX3", "PAX7", "KRAS",
                       "NRAS", "MYLPF",
                       "CDK4", "CDK6", "WEE1")
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
    print('----------')
    gene.list <- list2df(gene.list)
    colnames(gene.list) = c("cellMarker", "geneSymbol")
    tidyr::as_tibble(gene.list) %>% dplyr::mutate(geneID = strsplit(as.character(geneSymbol), ', '))
}


plotmarkers = function(x, meta, cluster='RNA_snn_res.0.8') {
    state = unlist(meta[gsub("MAST-", "MAST",as.character((x@meta.data)[1,1])), , drop=T])
    state = state[!((state=='') | (is.na(state)))]
    print(state)
    print(unique(x[[cluster]]))
    levels(x@meta.data[[cluster]]) = state
    clusters = x@meta.data[[cluster]]
    levels(clusters)[grepl('EMT', levels(clusters))] = 'EMT'
    cols = as.vector(clusters)
    print(cols)
    cols[!(clusters%in%c("Muscle", "EMT", "Prolif", "Ground"))] = "Other"
    print(cols)
    cols = factor(cols, levels=sort(c(unique(as.character(clusters[clusters%in%c("Muscle", "EMT", "Prolif", "Ground")])), "Other")))
    print(levels(cols))
    muscle_markers = Matrix::colMeans(x[paste0("hg19-", muscle),]$RNA@data)
    emt_markers = Matrix::colMeans(x[paste0("hg19-", emt), ]$RNA@data)
    prolif_markers = Matrix::colMeans(x[paste0("hg19-", prolif), ]$RNA@data)
    print(levels(cols))
    named_cols = c("red", "purple", "blue", "brown", "grey")
    names(named_cols) = c("Muscle", "EMT", "Prolif", "Ground", "Other")
    scatter3D(muscle_markers, emt_markers, prolif_markers,
              xlab="Muscle", ylab="EMT", zlab="Prolif",
              xlim=c(0, 2.0), ylim=c(0, 1.3), zlim=c(0, 1.3),
              bty = "g", alpha=0.7, ticktype = "detailed",
              theta = 135, phi = 40, pch = 16, cex=0.5, main=as.character((x@meta.data)[1,1]), colvar=as.integer(cols) , col = named_cols[levels(cols)],
              colkey = list(at = seq(1, length(levels(cols))), side = 1, 
                            addlines = TRUE, length = 0.6, width = 0.4,
                            labels = levels(cols)))
}

integration <- function(x, y, label) {
    Idents(x) = x$seurat_clusters2 = x$seurat_clusters
    states = unlist(as.vector(annotation[label[1], ]))
    states = states[states!='']
    states = states[!is.na(states)]
    print(states)
    levels(x$seurat_clusters) = states
    print(levels(x$seurat_clusters))
    Idents(y) = y$seurat_clusters2 = y$seurat_clusters
    states = unlist(as.vector(annotation[label[2], ]))
    states = na.omit(states[states!=''])
    print(states)
    print(length(states))
    print(levels(y$seurat_clusters))
    levels(y$seurat_clusters) = states
    seurat.pseudo.list = list(x, y)
    print(head(y$seurat_clusters))
    for (i in 1:length(seurat.pseudo.list)) {
        print(i)
        seurat.pseudo.list[[i]] <- NormalizeData(seurat.pseudo.list[[i]], verbose = FALSE)
        seurat.pseudo.list[[i]] <- FindVariableFeatures(seurat.pseudo.list[[i]],
                                                        selection.method = "vst", 
                                                        nfeatures = 2000, verbose = FALSE)
    }
    names(seurat.pseudo.list) <- label
    ## Update to latest api https://satijalab.org/seurat/articles/integration_introduction.html
    features <- SelectIntegrationFeatures(object.list = seurat.pseudo.list)
    anchors <- FindIntegrationAnchors(object.list = seurat.pseudo.list, anchor.features = features)
    integrated <- IntegrateData(anchorset = anchors)
    DefaultAssay(integrated) <- "integrated"
    integrated <- ScaleData(integrated, verbose = FALSE)
    integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
    integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
    integrated
}

annotation = read.delim('clusters.txt', sep='\t', row.names=1, header=T,
                        check.names=F, stringsAsFactors=F)

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
                rgb(0, 255, 253, maxColorValue = 255))
metalabels <- c("Ground", "Hypoxia", "EMT",
                "Prolif", "UNASSIGNED",
                "Prolif",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'UPR', "Neural")
names(metacolors) <- metalabels

## res <- readRDS('results/seurat_sara/MAST-39-Resistant_seurat-object.rds')
## sen <- readRDS('results/seurat_sara/MAST-39-Sensitive_seurat-object.rds')
## res$seurat_clusters = res$RNA_snn_res.0.8
## sen$seurat_clusters = sen$RNA_snn_res.0.8

res <- readRDS('../PDX/results/seurat_sara/MAST-139-Resistant_seurat-object.rds')
sen <- readRDS('../PDX/results/seurat_sara/MAST-139-Sensitive_seurat-object.rds')
res$seurat_clusters = res$RNA_snn_res.0.8
sen$seurat_clusters = sen$RNA_snn_res.0.8

conds = list(res, sen)
objs = integration(conds[[1]], conds[[2]], c('MAST139-Resistant', 'MAST139-Sensitive'))

pdf('ALL_MAST139_batchcorrected_UMAP.pdf', width=8, height=8)
p2=DimPlot(objs, ncol=1, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
## print(p1/p2)
print(p2)
dev.off()

DefaultAssay(objs) <- "RNA"
pdf('ALL_MAST139_markers.pdf', width=8, height=6)
p2=FeaturePlot(objs, features=paste0('hg19-', c("PIK3CA", "NDRG1", "AKT", "mTORC1", "NRF", "NPRG1", "ABC1", "NDR1")), slot='data', split.by='orig.ident',  cols = c("grey", "red"), by.col=F) + theme(legend.position='right')
print(p2)
dev.off()

muscle = scan("/PHShome/qq06/langenau/projects/01_sc_rms/final_annotations/gene_modules/Muscle.txt", what="")
emt = scan("/PHShome/qq06/langenau/projects/01_sc_rms/final_annotations/gene_modules/EMT.txt", what="")
prolif = scan("/PHShome/qq06/langenau/projects/01_sc_rms/final_annotations/gene_modules/Prolif.txt", what="")

pdf("MAST139_3Dplot.pdf", width=6.8, height=8)
par(mfrow=c(2, 1), cex=0.6)
for (index in seq_along(conds)) {
    print(index)
    plotmarkers(conds[[index]], annotation, cluster=c('RNA_snn_res.0.8', 'RNA_snn_res.0.8')[index])
}
dev.off()

DefaultAssay(objs) <- "integrated"
objs <- FindNeighbors(objs, reduction = "pca", dims = 1:30)
objs <- FindClusters(objs, resolution = 0.8)

pdf('ALL_MAST139_batchcorrected_UMAP_integrative_cluster.pdf', width=8, height=4)
p1=DimPlot(objs, reduction='umap', split.by='orig.ident', group.by='integrated_snn_res.0.8', label=T) + theme(legend.position='right')
print(p1)
dev.off()

DefaultAssay(objs) <- "RNA"
cluster.de = list()
for (cluster in levels(objs$integrated_snn_res.0.8)) {
    cluster.markers <- FindConservedMarkers(objs, ident.1 = cluster, grouping.var = "orig.ident", verbose = FALSE)
    cluster.markers$cluster = cluster
    cluster.de[[cluster]] = cluster.markers
}

cluster.de.all = do.call('rbind', cluster.de)

write.table(cluster.de.all, file="MAST139_resistant_sensitive_integrative_markers.xls", quote=F,
            sep='\t')

DefaultAssay(objs) <- "RNA"
mean.exp <- colMeans(x = objs@assays$RNA@data[paste0("hg19-", muscle), ], na.rm = TRUE)
objs@meta.data$muscle.score <- mean.exp
mean.exp <- colMeans(x = objs@assays$RNA@data[rownames(objs@assays$RNA@data)[(rownames(objs@assays$RNA@data)%in%paste0("hg19-", emt))], ], na.rm = TRUE)
objs@meta.data$emt.score <- mean.exp
mean.exp <- colMeans(x = objs@assays$RNA@data[paste0("hg19-", prolif), ], na.rm = TRUE)
objs@meta.data$prolif.score <- mean.exp

pdf("MAST139_average_signatures.pdf", width=12, height=4)
FeaturePlot(objs, features=c("muscle.score", "prolif.score", "emt.score"),
            slot='data', split.by='orig.ident',
            cols = c("grey", "red"), by.col=F) + theme(legend.position='right')
dev.off()

correct <- function(x) { return (1e6*x/sum(x)) }

gene.list = get_internal_geneset()
plot.genes = list()
for (i in gene.list$cellMarker) {
    plot.genes[[i]] = data.frame(module=i, gene=gene.list[gene.list$cellMarker == i, 'geneID']$geneID[[1]])
}
plot.genes = do.call(rbind, plot.genes)

library(clusterProfiler)
library(DOSE)
allgenes = readRDS("allgenes_signatures.rds")


cluster.de.all$Sensitive_diff_pct = cluster.de.all$`MAST-139-Sensitive_pct.1` - cluster.de.all$`MAST-139-Sensitive_pct.2`

cluster.de.all$Resistant_diff_pct = cluster.de.all$`MAST-139-Resistant_pct.1` - cluster.de.all$`MAST-139-Resistant_pct.2`

clusters = list()
for (i in unique(cluster.de.all$cluster)) {
    genes.sub = subset(cluster.de.all,
                       cluster==i & minimump_p_val<=0.01 & (`MAST-139-Resistant_avg_log2FC` > 0 | `MAST-139-Sensitive_avg_log2FC` > 0) & (Sensitive_diff_pct >= 0.1 | Resistant_diff_pct >= 0.1))
    genes = gsub("\\d*\\.hg19-", "", rownames(genes.sub))
    genes = genes[(!is.na(genes)) & (genes!="")]
    cat(i, "\t", length(genes), "\n")
    if (length(genes) != 0){
        y1 = enricher(genes, TERM2GENE=allgenes[,-2], maxGSSize=1000, minGSSize = 5)
        y1 = y1@result
        rownames(y1) = paste0(i, ".", rownames(y1))
        y1$clusters = i
        clusters[[as.character(i)]] <- y1
    }
}

clusters = do.call('rbind', clusters)
clusters = clusters %>% mutate(FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

write.table(clusters, file="MAST139_resistant_sensitive_integrative_markers_gsea.xls", quote=F,
            sep='\t')

ann_all = read.delim("MAST139/MAST139_integrative_clusters.txt", sep='\t', header=T, row.names=1)

plotmarkers.integrate = function(x, cluster='seurat_clusters') {
    clusters = x@meta.data[[cluster]]
    levels(clusters)[grepl('EMT', levels(clusters))] = 'EMT'
    cols = as.vector(clusters)
    cols[!(clusters%in%c("Muscle", "EMT", "Prolif", "Ground"))] = "Other"
    cols = factor(cols, levels=sort(c(unique(as.character(clusters[clusters%in%c("Muscle", "EMT", "Prolif", "Ground")])), "Other")))
    muscle_markers = Matrix::colMeans(x[paste0("hg19-", muscle),]$RNA@data)
    emt_markers = Matrix::colMeans(x[paste0("hg19-", emt), ]$RNA@data)
    prolif_markers = Matrix::colMeans(x[paste0("hg19-", prolif), ]$RNA@data)
    muscle_markers = muscle_markers/max(muscle_markers)
    emt_markers = emt_markers/max(emt_markers)
    prolif_markers = prolif_markers/max(prolif_markers)
    named_cols = c("red", "purple", "blue", "brown", "grey")
    names(named_cols) = c("Muscle", "EMT", "Prolif", "Ground", "Other")
    scatter3D(muscle_markers, emt_markers, prolif_markers,
              xlab="Muscle", ylab="EMT", zlab="Prolif",
              bty = "g", alpha=0.7, ticktype = "detailed",
              theta = 135, phi = 40, pch = 16, cex=0.5, main=as.character((x@meta.data)[1,1]), colvar=as.integer(cols) , col = named_cols[levels(cols)],
              colkey = list(at = seq(1, length(levels(cols))), side = 1, 
                            addlines = TRUE, length = 0.6, width = 0.4,
                            labels = levels(cols)))
}

## after annotation
pdf('ALL_MAST139_batchcorrected_UMAP_integrative_cluster_allcellstates.pdf', width=8, height=4)
Idents(objs) = objs$seurat_clusters = objs$integrated_snn_res.0.8
levels(objs$seurat_clusters) = unlist(ann_all[1,,drop=T])
p1=DimPlot(objs, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
print(p1)
dev.off()

print(table(objs$orig.ident, objs$seurat_clusters))

pdf('ALL_MAST139_batchcorrected_UMAP_integrative_cluster_dominantstates.pdf', width=8, height=4)
Idents(objs) = objs$seurat_clusters = objs$integrated_snn_res.0.8
levels(Idents(objs)) = levels(objs$seurat_clusters) = unlist(ann_all[2,,drop=T])
p1=DimPlot(objs, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors[unique(unlist(ann_all[2, , drop=T]))]) + theme(legend.position='right')
print(p1)
dev.off()

print(table(objs$orig.ident, objs$seurat_clusters))

DefaultAssay(objs) <- "RNA"

pdf("MAST139_3Dplot_integrate.pdf", width=8, height=5)
par(mfrow=c(1, 2), cex=0.6)
plotmarkers.integrate(objs[, objs@meta.data$orig.ident=='MAST-139-Resistant'])
plotmarkers.integrate(objs[, objs@meta.data$orig.ident=='MAST-139-Sensitive'])
dev.off()

pdf("ALL_MAST139_integrate_markers.pdf")
FeaturePlot(objs, features = paste0('hg19-', c("CD44", "MYOG", "MYOD1", "BUB3")), min.cutoff = "q9")
dev.off()

cpm = as.data.frame(apply(as.matrix(objs$RNA@counts), 2, correct))

cpm.sub = cpm[paste0("hg19-", as.character(plot.genes[,2])), ][, objs@meta.data$orig.ident=='MAST-139-Resistant']
sortcells = order(objs[, objs@meta.data$orig.ident=='MAST-139-Resistant']$integrated_snn_res.0.8)
heatdata  <- cpm.sub[, sortcells]
clusters <- objs[, objs@meta.data$orig.ident=='MAST-139-Resistant']$integrated_snn_res.0.8[sortcells]
annrow <- plot.genes[, 1, drop=F]
rownames(annrow) <- paste0(plot.genes[, 1], '.', plot.genes[, 2])

anncol <- data.frame(cluster=clusters)
mat2 = t(apply(heatdata, 1, function(x) {
    q10 <- quantile(x, 0.1, na.rm=T)
    q90 <- quantile(x, 0.9, na.rm=T)
    x[x < q10] <- q10
    x[x > q90] <- q90
    x = (x - mean(x)) / sd(x) ^ as.logical(sd(x))
    ## ## x = log2(x+1)
    ## scale(x)
}))
selection = complete.cases(mat2)
mat2 = mat2[selection, ]
rownames(mat2) = rownames(annrow)[selection]
colnames(mat2) = rownames(anncol)

pdf(glue("MAST139_resistant_heatmap.pdf"), width=22, height=12)
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
             column_title = qq("MAST139 resistant expression for @{ncol(heatdata)} cells"),
             top_annotation=topha,
             col = colorRamp2(c(min(mat2), -0.6, -0.4, 0, max(mat2)),
                              c("blue",
                                rgb(97, 233, 234, maxColorValue = 255),
                                "white",
                                "white", "red")),
             left_annotation=leftha)
draw(ha)
dev.off()


cpm.sub = cpm[paste0("hg19-", as.character(plot.genes[,2])), ][, objs@meta.data$orig.ident=='MAST-139-Sensitive']
sortcells = order(objs[, objs@meta.data$orig.ident=='MAST-139-Sensitive']$integrated_snn_res.0.8)
heatdata  <- cpm.sub[, sortcells]
clusters <- objs[, objs@meta.data$orig.ident=='MAST-139-Sensitive']$integrated_snn_res.0.8[sortcells]

annrow <- plot.genes[, 1, drop=F]
rownames(annrow) <- paste0(plot.genes[, 1], '.', plot.genes[, 2])
anncol <- data.frame(cluster=clusters)
mat2 = t(apply(heatdata, 1, function(x) {
    q10 <- quantile(x, 0.1, na.rm=T)
    q90 <- quantile(x, 0.9, na.rm=T)
    x[x < q10] <- q10
    x[x > q90] <- q90
    scale(x)
}))
selection = complete.cases(mat2)
mat2 = mat2[selection, ]
rownames(mat2) = rownames(annrow)[selection]
colnames(mat2) = rownames(anncol)

pdf(glue("MAST139_sensitive_heatmap.pdf"), width=22, height=12)
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
             column_title = qq("MAST139 resistant expression for @{ncol(heatdata)} cells"),
             top_annotation=topha,
             col = colorRamp2(c(min(mat2), -0.6, -0.4, 0, max(mat2)),
                              c("blue",
                                rgb(97, 233, 234, maxColorValue = 255),
                                "white",
                                "white", "red")),
             left_annotation=leftha)
draw(ha)
dev.off()


objs$celltype.stim <- paste(objs$seurat_clusters, objs$orig.ident, sep = "_")
Idents(objs) <- "celltype.stim"

responses = list()
for (state in levels(objs$seurat_clusters)) {
    response <- FindMarkers(objs,
                            ident.1 = paste0(state, "_MAST-139-Resistant"),
                            ident.2 = paste0(state, "_MAST-139-Sensitive"),
                            verbose = FALSE)
    responses[[state]] = response
}

responses.filter = list()
for (state in levels(objs$seurat_clusters)) {
    responses[[state]]$diff.pct = responses[[state]]$pct.1 - responses[[state]]$pct.2
    responses.filter[[state]] = subset(responses[[state]],
                                       p_val_adj<0.05 & abs(avg_log2FC) > 0.5 & abs(diff.pct) > 0.05)
    write.table(responses.filter[[state]], file=paste0("MAST139_resistant_sensitive_integrative_diffgenes", state, ".xls"), quote=F,
            sep='\t')
}

pdf(glue("MAST139_differential_markers_EMT.pdf"), width=9, height=28)
for (chunk in split(rownames(head(responses.filter[['EMT']], 50)), ceiling(seq(50)/9))) {
    print(chunk)
    print(FeaturePlot(objs, features = chunk, ncol=3,
                      split.by = "orig.ident", max.cutoff = 3, cols = c("grey", "red")))
}
dev.off()


pdf(glue("MAST139_differential_markers_Muscle.pdf"), width=9, height=28)
for (chunk in split(rownames(head(responses.filter[['Muscle']], 50)), ceiling(seq(50)/9))) {
    print(chunk)
    print(FeaturePlot(objs, features = chunk, ncol=3,
                      split.by = "orig.ident", max.cutoff = 3, cols = c("grey", "red")))
}
dev.off()

pdf(glue("MAST139_differential_markers_Prolif.pdf"), width=9, height=28)
for (chunk in split(rownames(head(responses.filter[['Prolif']], 50)), ceiling(seq(50)/9))) {
    print(chunk)
    print(FeaturePlot(objs, features = chunk, ncol=3,
                      split.by = "orig.ident", max.cutoff = 3, cols = c("grey", "red")))
}
dev.off()

pdf(glue("MAST139_differential_markers_Ground.pdf"), width=9, height=28)
for (chunk in split(rownames(head(responses.filter[['Ground']], 50)), ceiling(seq(50)/9))) {
    print(chunk)
    print(FeaturePlot(objs, features = chunk, ncol=3,
                      split.by = "orig.ident", max.cutoff = 3, cols = c("grey", "red")))
}
dev.off()

cat(intersect(intersect(intersect(rownames(head(responses.filter[['Ground']], 50)),
                                  rownames(head(responses.filter[['EMT']], 50))),
          rownames(head(responses.filter[['Muscle']], 50))),
          rownames(head(responses.filter[['Prolif']], 50))), "\n")
