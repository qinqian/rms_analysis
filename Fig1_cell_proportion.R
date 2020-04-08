library(Seurat)
library(GetoptLong)
library(RColorBrewer)
library(ComplexHeatmap)
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

results = list()
for (i in seq_along(labels)) {
    states = unlist(as.vector(annotation[labels[i], ]))
    states = states[states!='']
    print(states)
    levels(pdxs.objs[[i]]$RNA_snn_res.0.8) = states
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
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "INTERFERON", "Prolif",
                "Histone", "TNFA")

names(metacolors) <- metalabels

cols = metacolors[levels(pdxs.objs[[1]]$RNA_snn_res.0.8)]

pdf("Fig1_MAST111.pdf", width=4.6, height=3.5)
DimPlot(pdxs.objs[[1]], group.by='RNA_snn_res.0.8', label=F,
        cols=cols, legend=F)
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

library(ComplexHeatmap)
tumor = pdxs.objs[[1]]
df = read.table(paste0('../results/seurat_sara/20190624_seurat-object_MAST111_SCT_res0.8.xls'))
tumor$RNA_snn_res.0.8 = reorder(tumor$RNA_snn_res.0.8,
                                new.order=c('Prolif', 'Muscle', 'Hypoxia', 'EMT', 'TNFA', 'Ground'))
gene.list[,1] = reorder(droplevels(gene.list[,1]),
                        new.order=c('Prolif', 'Muscle', 'Hypoxia', 'EMT', 'TNFA'))


sortcells <- order(tumor$RNA_snn_res.0.8)
source('DEGs_seurat3_sara.R')

cpm = as.data.frame(apply(as.matrix(tumor$RNA@counts), 2, correct))

sortgenes = order(gene.list[, 1])
gene.list = gene.list[sortgenes,]

heatdata  <- cpm[as.character(gene.list[,2]), sortcells]
annrow <- gene.list[, 1, drop=F]
## rownames(annrow) <- paste0(gene.list[, 1], '.', gene.list[, 2])
rownames(annrow) <- gene.list[, 2]
anncol <- data.frame(cluster=tumor$RNA_snn_res.0.8[sortcells])
## ha = structure(brewer.pal(length(unique(anncol[, 1])), "Set3"),
##                names=levels(anncol[, 1]))
colsrow = cols[as.vector(annrow[,1])]
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

pdf('Fig1_MAST111heatmap.pdf', width=10)
## mat2= t(scale(t(heatdata)))
## mat2 = t(scale(t(log2(heatdata+1))))
topha = HeatmapAnnotation(clusters=as.vector(anncol[,1]), 
                          col=list(clusters=cols[anncol[,1]]), show_legend=T,
                          gp = gpar(col = NA))
leftha = rowAnnotation(modules=annrow[, 1],
                       col=list(modules=colsrow),
                       gp = gpar(col = NA))
## Heatmap(as.vector(annrow[,1]), col=colsrow, width = unit(0.5, "cm"), name='clusters')+
Heatmap(mat2, name = "Scaled Expression", 
        show_column_dend = FALSE, cluster_rows=F, cluster_columns=F,
        show_row_names=F,
        show_column_names = FALSE, width = unit(15, "cm"), 
        column_title = qq("MAST111 relative expression for @{ncol(heatdata)} cells"),
        col = colorRamp2(c(-1.5, -0.6, -0.5, 0, 1.5),
                         c("blue", rgb(97, 233, 234, maxColorValue = 255),
                           "white", "white", "red")),
        top_annotation=topha,
        left_annotation=leftha)
dev.off()
