library(Seurat)
library(GetoptLong)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(Seurat)
library(glue)
library(ComplexHeatmap)
library(GetoptLong)
library(ColorBrewer)
library(ggplot2)
library(circlize)

pdxs = Sys.glob('../data/seurat_obj/*rds')[1:10]

annotation = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1, header=T,
                        check.names=F, stringsAsFactors=F)

pdxs = c(pdxs, '../results/seurat_sara/20191031_MSK74711_seurat-object.rds',
         '../results/seurat_sara/MAST118_seurat-object.rds',
         '../results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds',
         '../results/seurat_sara/MAST139_1cells_seurat-object.rds')
pdxs.objs = lapply(pdxs, readRDS)

labels = unlist(lapply(pdxs.objs, function(x) {
    levels(x$orig.ident[1])
}))

labels[9] = 'MAST85-1'
labels[11] = 'MSK74711'
labels[13] = 'MSK72117'
labels[14] = 'MAST139-1'
labels[10] = 'RH74-10'

allpdx.meta = do.call('rbind', list(pdxs.objs[[1]]@meta.data,
                                    pdxs.objs[[2]]@meta.data,
                                    pdxs.objs[[3]]@meta.data,
                                    pdxs.objs[[4]]@meta.data,
                                    pdxs.objs[[5]]@meta.data,
                                    pdxs.objs[[6]]@meta.data,
                                    pdxs.objs[[7]]@meta.data
                                    ))


allpdx.meta2 = do.call('rbind', list(pdxs.objs[[11]]@meta.data,
                                     pdxs.objs[[12]]@meta.data,
                                     pdxs.objs[[13]]@meta.data
                                     ))

nrow(allpdx.meta) + nrow(allpdx.meta2)

summary (c (allpdx.meta [, c ('nFeature_RNA')], allpdx.meta2[, c ('nFeature_RNA')]))
sd (c (allpdx.meta [, c ('nFeature_RNA')], allpdx.meta2[, c ('nFeature_RNA')]))

results = list()
for (i in seq_along(labels)) {
    states = unlist(as.vector(annotation[labels[i], ]))
    states = states[states!='']
    ## print(states)
    pdxs.objs[[i]]$seurat_clusters = pdxs.objs[[i]]$RNA_snn_res.0.8
    levels(pdxs.objs[[i]]$RNA_snn_res.0.8) = states
    print(levels(pdxs.objs[[i]]$RNA_snn_res.0.8))
    results[[labels[i]]] = table(pdxs.objs[[i]]$RNA_snn_res.0.8)
}

cols = unique(unlist(lapply(results, names)))
cell.prop = matrix(0, ncol=length(cols), nrow=length(pdxs.objs))
colnames(cell.prop) = cols
rownames(cell.prop) = labels

for (i in seq_along(labels)) {
    cell.prop[labels[i], names(results[[labels[i]]])] = results[[labels[i]]]
}

write.table(cell.prop, file='PDX_all_cells.xls', quote=F, sep='\t', col.names=NA)

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
                'gray',
                'yellow')
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "INTERFERON", "Prolif",
                "Histone", "Apoptosis",
                'ARMS core', 'ERMS core')

core.sig = read.table('tables_storage/RMS_core_t_test_pval0.05_fold1.5.xls')
test = rowMeans(core.sig[,1:4]) - rowMeans(core.sig[,5:11])
erms.topsign = rownames(core.sig)[test<0]
arms.topsign = rownames(core.sig)[test>0]

names(metacolors) <- metalabels
cols = metacolors[levels(pdxs.objs[[1]]$RNA_snn_res.0.8)]

pdf("Fig1_MAST111.pdf", width=9.2, height=3.5)
p1=DimPlot(pdxs.objs[[1]], group.by='seurat_clusters', label=F)
p2=DimPlot(pdxs.objs[[1]], group.by='RNA_snn_res.0.8', label=F,
           cols=cols)
print(CombinePlots(plots=list(p1, p2), ncol=2))
dev.off()

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

gene.list <- list2df(gene.list)
gene.list <- gene.list[gene.list[,1]%in%levels(pdxs.objs[[1]]$RNA_snn_res.0.8), ]

core.sig = read.table('tables_storage/RMS_core_t_test_pval0.05_fold1.5.xls')
test = rowMeans(core.sig[,1:4]) - rowMeans(core.sig[,5:11])
gene.list = rbind(gene.list, data.frame(module='ERMS core', gene=rownames(core.sig)[test<0]))
gene.list = rbind(gene.list, data.frame(module='ARMS core', gene=rownames(core.sig)[test>0]))

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(gdata)

tumor = pdxs.objs[[1]]
df = read.table(paste0('../results/seurat_sara/20190624_seurat-object_MAST111_SCT_res0.8.xls'))
tumor$RNA_snn_res.0.8 = reorder(tumor$RNA_snn_res.0.8,
                                new.order=c('Prolif', 'Muscle', 'Hypoxia', 'EMT', 'Apoptosis', 'Ground'))

gene.list[,1] = reorder(droplevels(gene.list[,1]),
                        new.order=c('Prolif', 'Muscle', 'Hypoxia', 'EMT', 'Apoptosis', 'ARMS core', 'ERMS core'))

## gene.list[,1] = reorder(droplevels(gene.list[,1]),
##                         new.order=c('Prolif', 'Muscle', 'Hypoxia', 'EMT', 'Apoptosis'))


source('DEGs_seurat3_sara.R')

cpm = as.data.frame(apply(as.matrix(tumor$RNA@counts), 2, correct))

sortgenes = order(gene.list[, 1])
gene.list = gene.list[sortgenes,]

sortcells <- order(tumor$RNA_snn_res.0.8, tumor$seurat_clusters)
## sortcells = order(tumor$seurat_clusters, tumor$RNA_snn_res.0.8)
heatdata  <- cpm[as.character(gene.list[,2]), sortcells]
clusters <- tumor$seurat_clusters[sortcells]

test=(cbind(as.vector(tumor$RNA_snn_res.0.8), as.vector(tumor$seurat_clusters)))
test=table(test[,1], test[,2])

cluster_cols=rep('gray', 9)
names(cluster_cols)= 0:8
for (i in seq_along(rownames(test))) {
    if (sum(test[i, ]>0)>1) {
        if (i%%2)
            cluster_cols[test[i, ]>0] = c('gray', 'white')
        else
            cluster_cols[test[i, ]>0] = c('white', 'gray')
    } else {
        if (i%%2)
            cluster_cols[test[i, ]>0] = c('white')
        else
            cluster_cols[test[i, ]>0] = c('gray')
    }
}

annrow <- gene.list[, 1, drop=F]
## rownames(annrow) <- paste0(gene.list[, 1], '.', gene.list[, 2])
rownames(annrow) <- paste0(gene.list[,1], gene.list[, 2])

anncol <- data.frame(cluster=tumor$RNA_snn_res.0.8[sortcells])

## ha = structure(brewer.pal(length(unique(anncol[, 1])), "Set3"),
##                names=levels(anncol[, 1]))
## colsrow = cols[as.vector(annrow[,1])]
## mat2 = t(apply(log2(heatdata+1), 1, function(x) {
mat2 = t(apply(heatdata, 1, function(x) {
                                         q10 <- quantile(x, 0.1)
                                         q90 <- quantile(x, 0.9)
                                         x[x < q10] <- q10
                                         x[x > q90] <- q90
                                         scale(x, scale=T)
                                         ##    (x-mean(x))/ sd(x) ^ as.logical(sd(x))
                                         }))
selection = complete.cases(mat2)
mat2 = mat2[selection, ]
annrow = annrow[selection, ,drop=F]
## anncol = anncol[selection, ,drop=F]

names (metacolors) = metalabels
## pdf('Fig1_MAST111heatmap.pdf', width=10, height=8)
tiff('Fig1_MAST111heatmap.tif', units="in", width=18, height=6, res=320)
topha = HeatmapAnnotation(states=as.vector(anncol[,1]),
                          clusters=clusters,
                          col=list(states=cols,
                                   clusters=cluster_cols),
                          show_legend=T,
                          gp = gpar(col = NA),
                          border = c(states=F, clusters=T))
leftha = rowAnnotation(modules=as.vector(annrow[,1]),
                       ## col=list(modules=cols),
                       col = list (modules=metacolors),
                       gp = gpar(col = NA))
ha = Heatmap(mat2, name = "Scaled Expression",
             use_raster = TRUE, raster_quality = 2,
             show_column_dend = FALSE, cluster_rows=F, cluster_columns=F,
             show_row_names=F,
             show_column_names = FALSE, width = unit(10, "cm"),
             column_title = qq("MAST111 relative expression for @{ncol(heatdata)} cells"),
             col = colorRamp2(c(-1, -0.6, -0.4, 0, 1),
                              c("blue", rgb(97, 233, 234, maxColorValue = 255),
                                "white",
                                ## rgb(97, 233, 234, maxColorValue = 255),
                                "white", "red")),
             top_annotation=topha,
             left_annotation=leftha)
draw(ha)
dev.off()

