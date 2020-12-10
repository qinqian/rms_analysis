library(Seurat)
library(patchwork)
library(glue)
library(ComplexHeatmap)
library(GetoptLong)
library(circlize)

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

rd = readRDS('../results/seurat_sara/RD_seurat-object.rds')

annotation <- read.delim('../final_annotations/cellline_annotation.txt', row.names=1,
                         stringsAsFactors = F)

label = 'RD'
states = as.vector(unlist(as.vector(annotation[label, ])))
states = as.vector(states[!is.na(states)])

rd$seurat_clusters = rd$RNA_snn_res.0.8

levels(rd$RNA_snn_res.0.8) = states

write.csv(table(rd$RNA_snn_res.0.8), 'RD_cellprop.csv')

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
gene.list <- gene.list[as.vector(gene.list[,1])%in%c(levels(rd$RNA_snn_res.0.8), 'ARMS_core', 'ERMS_core'), ]

geneord = names(sort(table(rd$RNA_snn_res.0.8)))

library(ggplot2)
library(circlize)
library(gdata)

gene.list[,1] = reorder(droplevels(gene.list[,1]),
                        new.order=c("EMT", "Muscle", "Prolif", "Hypoxia",
                                    'ARMS_core', 'ERMS_core'))

source('DEGs_seurat3_sara.R')

cpm = as.data.frame(apply(as.matrix(rd$RNA@counts), 2, correct))

sortgenes = order(gene.list[, 1])
gene.list = gene.list[sortgenes,]
sortcells = order(rd$RNA_snn_res.0.8, as.vector(rd$seurat_clusters))
heatdata  <- cpm[as.character(gene.list[,2]), sortcells]
clusters <- rd$seurat_clusters[sortcells]
annrow <- gene.list[, 1, drop=F]
rownames(annrow) <- paste0(gene.list[, 1], '.', gene.list[, 2])
anncol <- data.frame(cluster=rd$RNA_snn_res.0.8[sortcells])

## mat2= t(scale(t(heatdata)))
mat2 = t(apply(heatdata, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))

(a<-levels(clusters))
test=cbind(as.vector(rd$seurat_clusters)[sortcells], as.vector(rd$RNA_snn_res.0.8)[sortcells])
test=table(test[,1], test[,2])
test=t(test[,levels(rd$RNA_snn_res.0.8[sortcells])])
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

tiff(glue('Fig4C_{label}_heatmap.tiff'), units="in", width=18, height=4.5, res=320)
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
                                "white", "red")),
             top_annotation=topha,
             left_annotation=leftha)
draw(ha)
dev.off()

pdf(glue('Fig4C_{label}_umap.pdf'), width=5, height=4.5)
print(DimPlot(rd, group.by='RNA_snn_res.0.8', cols=metacolors))
dev.off()
