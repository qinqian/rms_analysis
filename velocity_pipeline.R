library(argparse)
library(org.Hs.eg.db)
library(velocyto.R)
library(Seurat)
library(gdata)
library(tidyverse)
library(Matrix)
library(clustree)
library(pheatmap)
library(foreach)
library(SeuratWrappers)
library(UpSetR)
set.seed(100)

## human_ortholog = read.table('~/alvin_singlecell/01_rms_projects//01_fish/data/ortholog_mapping/Beagle_fish_human_all_genes.txt', header=T, sep='\t', stringsAsFactors=F)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--seuratobj', dest='seurat', metavar='N', type="character", nargs='+')
    parser$add_argument('--velobj', dest='vel', metavar='N', type="character", nargs='+')
    parser$add_argument('--label', dest='label', default='')
    parser$add_argument('--clusterlabel', dest='clusterlabel', nargs='+')

    parser$add_argument('--finalres', dest='res', default=0.8, type='double')
    parser$add_argument('--species', dest='species', type='character', default='human')
    args = parser$parse_args()
    args
}

args = get_args()
if (length(args$seurat) == 0 || args$vel == '' || args$label == '') {
    cat('empty argument, exit..')
    q()
}


#args$seurat = '../results/seurat_intersect_velocity/Tumor24_seu.rds'
#args$vel = '../results/seurat_intersect_velocity/Tumor24_vel.rds'
#args$label = 'Tumor24'
#args$clusterlabel = 'Tumor24'
#args$species = 'fish'
## args$label = '20190624_seurat-object_MAST111_seu.rds'

obj = readRDS(args$seurat)
vel = readRDS(args$vel)
name = args$label

if (length(names(vel@tools))!=0) {
    print('velocity exist...')
} else {
    vel <- RunVelocity(object=vel, deltaT = 1, kCells = 25, fit.quantile = 0.02)
    saveRDS(vel, file=args$vel)
}

velocity = Tool(object=vel, slot = "RunVelocity")

if (args$species == 'human') {
    metacolors = c(rgb(166, 166, 166, maxColorValue = 255),
                   rgb(241, 149, 69, maxColorValue = 255),
                   rgb(103, 35,  102, maxColorValue = 255),
                   rgb(52, 101, 252, maxColorValue = 255),
                   rgb(242, 242, 242, maxColorValue = 255),
                   rgb(52, 101, 252, maxColorValue = 255),
                   rgb(233, 63,  51, maxColorValue = 255),
                   rgb(65, 129,  7, maxColorValue = 255),
                   rgb(52, 101, 252, maxColorValue = 255),
                   rgb(253, 247, 49, maxColorValue = 255))
    metalabels = c("GROUND", "Hypoxia", "EMT",  "G1S", "UNASSIGNED", "G2M",  "MUSCLE", "INTERFERON", "PROLIF", "Histones")
    colortab = read.table('color_table.xls', sep='\t', header=T)
    cluster = as.character(colortab[,args$clusterlabel])
    labels = na.omit(metalabels[match(cluster, metalabels)])
    colors = na.omit(metacolors[match(cluster, metalabels)])
    ident.colors = colors
    names(ident.colors) <- levels(obj@meta.data$RNA_snn_res.0.8)
    cell.colors         <- ident.colors[obj@meta.data$RNA_snn_res.0.8]
} else {

    metacolors = c(rgb(103, 35,  102, maxColorValue = 255),
                   rgb(52, 101, 252, maxColorValue = 255),
                   rgb(241, 149, 69, maxColorValue = 255),
		   rgb(166, 166, 166, maxColorValue = 255),
                   rgb(233, 63,  51, maxColorValue = 255),
                   rgb(52, 101, 252, maxColorValue = 255),
                   rgb(65, 129,  7, maxColorValue = 255),
                   rgb(242, 242, 242, maxColorValue = 255),
                   rgb(103, 35,  102, maxColorValue = 255),
                   rgb(52, 101, 252, maxColorValue = 255), 
                   rgb(103, 35,  102, maxColorValue = 255),
                   rgb(241, 149, 69, maxColorValue = 255),
                   rgb(233, 63,  51, maxColorValue = 255), 
                   rgb(52, 101, 252, maxColorValue = 255)
    )

    print(metacolors)
    colortab = read.delim('fish_color_table.txt', sep='\t', header=T, stringsAsFactors=F)
    metalabels = c("EMT", "MYC-N", "Hypoxia", "Ground", "MUSCLE", "Cell_cycle", "TNF", 'Unassigned', "EMT_ECM", "RhoA_Cell_cycle", "ECM_Invasion", "Hypoxia_TNF_SC", "Muscle_Prolif", "Prolif")

    print(metalabels)
    metalabels = metalabels[metalabels!='']
    cluster = na.omit(as.character(colortab[, args$clusterlabel]))
    cluster = cluster[cluster!='']
    labels = metalabels[match(cluster, metalabels)]
    colors = metacolors[match(cluster, metalabels)]
    ident.colors = colors
    names(ident.colors) <- levels(obj@meta.data$seurat_clusters)
    cell.colors         <- ident.colors[obj@meta.data$seurat_clusters]
}

if (args$species == 'human') {
    ### for human use the seurat umap embedding from cDNA read count by log transformation by sara !!!
    embedvel = Embeddings(object=vel, reduction = "umap")
    embed = Embeddings(object=obj, reduction = "umap")
    names(cell.colors) <- rownames(embedvel)
    rownames(embed) = rownames(embedvel)
} else {
    ### for zebrafish use the seurat umap from splicing read count by SCT transformation by alvin !!!
    embedvel = Embeddings(object=vel, reduction = "umap")
    embed = Embeddings(object=vel, reduction = "umap")
    names(cell.colors) <- rownames(embedvel)
    rownames(embed) = rownames(embedvel)
}

pdf(paste0('../results/', name, '_velocity_tumoronly.pdf'), width=16, height=9)
## par(mfrow=c(1, 2), font=2, cex=1.3)
par(font=2, cex=1.3, xpd=T, mar=c(4, 4, 3, 15))
show.velocity.on.embedding.cor(emb = embed,
    			       vel = velocity,
    			       n = 200, scale = "sqrt", 
    			       cell.colors = ac(cell.colors, alpha = 0.5),
    			       cex=1, n.cores=6,
    			       arrow.scale = 1.2, 
    			       show.grid.flow = TRUE, 
    			       min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, do.par = F, cell.border.alpha = 0.3)
legend(max(embed[,1])*1.15, max(embed[,2])*1.0, legend=paste0(names(ident.colors), '_', labels), pch=19, col=ident.colors, bty = "n")
title(args$label)
dev.off()

