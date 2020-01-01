library(argparse)
library(velocyto.R)
library(Seurat)
library(gdata)
library(tidyverse)
library(Matrix)
#library(clustree)
#library(pheatmap)
#library(foreach)
library(SeuratWrappers)
#library(UpSetR)
set.seed(100)
library(foreach)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--seuratobj1', dest='seurat1', metavar='N', type="character", nargs='+')
    parser$add_argument('--label', dest='label', default='')

    parser$add_argument('--species', dest='species', type='character', default='human')
    args = parser$parse_args()
    args
}

args = get_args()

if (length(args$seurat1) == 0 || args$label == '') {
    cat('empty argument, exit..')
    q()
}

## args$seurat1 <- "/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds"
## args$label <- "MAST111"

#args$seurat1 <- "../results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds"
#args$label <- "20191031_MSK72117tencell"

source("DEGs_seurat3_sara.R")
seurat1 <- readRDS(args$seurat1) ## sara's version
print(dim(seurat1))

Idents(seurat1) <- seurat1@meta.data$seurat_clusters <- seurat1@meta.data$RNA_snn_res.0.8
allcluster <- names((seurat1@meta.data$seurat_clusters) %>% table())

clusterde <- list()
for (i in allcluster) {
    print(i)
    print(allcluster[-(as.integer(i)+1)])
    de.up <- get_upregulated_genes_in_clusters(seurat1, i, allcluster[-(as.integer(i)+1)])
    de.up$cluster <- i
    clusterde[[i]] <- de.up
}
seurat1.de <- do.call('rbind', clusterde) ## still too many genes
## again, filter by adjusted p value and fraction of cells expressing the genes
seurat1.de <- subset(seurat1.de, p.adjusted <= 0.01 & (fg_fraction >= 0.1 | bg_fraction >= 0.1))

#seurat1.de = FindAllMarkers(seurat1, only.pos=T, min.pct=0.1,
#                            test.use='MAST', min.diff.pct=0.1,
#                            random.seed=100, logfc.threshold=0.2)
#seurat1.de$enrichment = seurat1.de$pct.1 - seurat1.de$pct.2
#seurat1.de = subset(seurat1.de, p_val_adj <= 0.01)

system('mkdir -p ../results/seurat_sara/')
write.table(seurat1.de, file=paste0('../results/seurat_sara/', args$label, '_SCT_res0.8.xls'), sep='\t', quote=F, col.names=NA)
## write.table(seurat1.de, file=paste0('../results/seurat_sara/', args$label, '_sara_v4_res0.8.xls'), sep='\t', quote=F, row.names=F)

pdf(paste0('../results/seurat_sara/', args$label, '_v4_comparison.pdf'), width=14, height=6)
p1 = DimPlot(seurat1, reduction='umap', label=T, group.by='RNA_snn_res.0.8') ## Sara umap with sara cluster
p2 = DimPlot(seurat1, reduction='umap', label=T, group.by='seurat_clusters') ## Alvin umap with sara cluster
p = CombinePlots(plots=list(p1, p2), ncol=2)
print(p)
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
library(pheatmap)
gene.list[, 2] <- as.vector(gene.list[, 2])
gene.list <- gene.list[gene.list[, 2] %in% rownames(seurat1$RNA@scale.data), ]
sortcells <- sort(seurat1$seurat_clusters)
heatdata  <- as.matrix(seurat1$RNA@scale.data)[gene.list[, 2], order(seurat1$seurat_clusters)]

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

pheatmap(heatdata, annotation_row=annrow, annotation_col=anncol,
         scale='none', annotation_colors=ann_colors,
         show_rownames = F, show_colnames = F, filename=paste0('../results/seurat_sara/', args$label, '_v4_heatmap_withmodule.pdf'), width=9, height=9,
         breaks=seq(-2, 2, length.out=31),
         cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("blue", "white", "darkred"))(30))

## pdf(paste0('../results/seurat_sara/', args$label, '_v4_heatmap.pdf'), width=22, height=30)
## par(mar=c(3,5,8,3))
#heat = DoHeatmap(subset(seurat1, downsample = 100), features = seurat1.de$gene, size = 3) + theme(
#do not downsample cells
## heat = DoHeatmap(seurat1, features = seurat1.de$gene, size = 5) + theme(
##   panel.background = element_rect(fill = "white"),
##   plot.margin = margin(2, 3, 5, 1, "cm"),
##   plot.background = element_rect(
##     fill = "white",
##     colour = "white",
##     size = 1
##   )
## )
## print(heat)
## dev.off()
