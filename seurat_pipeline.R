library(argparse)
library(Seurat)
library(SeuratWrappers)
## library(gdata)
## library(tidyverse)
## library(Matrix)
library(clustree)
library(pheatmap)
library(foreach)
## library(garnett)
library(fgsea)
## library(cellassign)
library(UpSetR)
library(clusterProfiler)
library(fgsea)
source("functions.R")
set.seed(100)
library(tidyverse)

## human_ortholog = read.table('~/alvin_singlecell/01_rms_projects//01_fish/data/ortholog_mapping/Beagle_fish_human_all_genes.txt', header=T, sep='\t', stringsAsFactors=F)
get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--seuratobj', dest='seurat', metavar='N', type="character", nargs='+')
    parser$add_argument('--label', dest='label', default='')
    parser$add_argument('--assaytype', dest='assaytype', default='RNA')
    parser$add_argument('--finalres', dest='res', default=0.05, type='double')
    parser$add_argument('--tumor', dest='tumor', metavar='T', type='integer', nargs='+')
    parser$add_argument('--species', dest='species', type='character', default='fish')
    parser$add_argument('--transform', dest='trans', type='character', default='SCT')
    args = parser$parse_args()
    args
}

args = get_args()
if (length(args$seurat) == 0 || args$label == '') {
    cat('empty argument, exit..')
    q()
}

## args$seurat = '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color/outs/filtered_feature_bc_matrix'
## args$seurat = '/PHShome/qq06/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/20190801_MAST85-1cell_5Kcells_hg19/outs/filtered_feature_bc_matrix'
## args$label = '20190801_MAST85-1cell_5Kcells_hg19'
## args$species = 'human'
## args$res = 0.8
## args$trans = 'SCT'
## args$trans = 'Log'
## args$res = 0.05
## args$assaytype = 'spliced'

all_markers_final = readRDS("../results/final_manual_markers.RDS")

if (length(args$seurat) > 1) {
    libraries = paste0("Library", seq(2))  # c("sort", 'bulk')
    if (args$assaytype == 'RNA') {
        seurat.obj1 <- Read10X(args$seurat[1])
        seurat.obj1 <- CreateSeuratObject(counts=seurat.obj1,
                                          project=libraries[1], min.cells=3, min.features=10)
        print('load 2')
        seurat.obj2 <- Read10X(args$seurat[2])
        seurat.obj2 <- CreateSeuratObject(counts=seurat.obj2,
                                          project=libraries[2], min.cells=3, min.features=10)
    } else {
        print('load 1')
        seurat.obj1 <- ReadVelocity(args$seurat[1])
        seurat.obj1 <- as.Seurat(seurat.obj1)
        print('load 2')
        seurat.obj2 <- ReadVelocity(args$seurat[2])
        seurat.obj2 <- as.Seurat(seurat.obj2)
    }
    seurat.obj  <- merge(seurat.obj1, seurat.obj2, add.cell.ids = libraries, project = args$label)
} else {
    if (args$assaytype == 'RNA') {
        seurat.obj <- Read10X(args$seurat[1])
        seurat.obj <- CreateSeuratObject(counts=seurat.obj,
                                         project=args$label, min.cells=3, min.features=10)
    } else {
        print(args$seurat)
        seurat.obj <- ReadVelocity(args$seurat[1])
        seurat.obj <- as.Seurat(seurat.obj)
    }
}

print(head(seurat.obj@meta.data))
print(tail(seurat.obj@meta.data))

if (args$species == 'fish') {
    seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern='^mt-')
} else {
    seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern='^MT-')
}

plot.qc = VlnPlot(seurat.obj, features=c(paste0("nFeature_", args$assaytype), paste0("nCount_", args$assaytype), "percent.mt"))
ggsave(paste0("../results/Cell_QC_violin_", args$label, ".pdf"), width=8, height=18)

pdf(paste0("../results/Cell_QC_histgram_", args$label, ".pdf"), width=18, height=12)
label = paste0('QC of ', args$label)
par(mfrow=c(2, 1), font=2, cex=1.5, mar=c(5, 5, 2, 0))
hist(seurat.obj@meta.data[[paste0('nFeature_', args$assaytype)]], xlab='Number of detected gene number', ylab="Cell barcode number",
     main=label, n=50, col="blue")
hist(seurat.obj@meta.data$percent.mt,  xlab="Mitochondria ratio", ylab="Cell barcode number",
     main=label, n=50, col="blue")
dev.off()

print(dim(seurat.obj))
if (args$assaytype=='RNA') {
    if (args$species == 'fish') {
        seurat.obj <- subset(seurat.obj, 
                             subset=nFeature_RNA>1000 & nFeature_RNA<4000 & percent.mt<10)
    } else {
        if (args$trans == 'SCT') {
            seurat.obj <- subset(seurat.obj, 
                                 subset=nFeature_RNA>1000 & nFeature_RNA<7000 & percent.mt<20)
        } else {
            seurat.obj <- subset(seurat.obj, 
                                 subset=nFeature_RNA>1000 & nFeature_RNA<7000 & percent.mt<10)
        }
    }
} else {
    if (args$species == 'fish') {
        seurat.obj <- subset(seurat.obj, 
                             subset=nFeature_spliced>1000 & nFeature_spliced<4000 & percent.mt<10)
    } else {
        cat('skip filtering cells for splicing matrix of human dataset', '\n')
        #seurat.obj <- subset(seurat.obj, 
        #                     subset=nFeature_spliced>1000 & nFeature_spliced<8000 & percent.mt<20)
    }
}

print(dim(seurat.obj))

print('test a................')
if (args$trans == 'SCT') {
    seurat.obj = process_standard(seurat.obj, norm=F, assaytype=args$assaytype,
                                  output=paste0('../results/', args$label, '_pc_jackstraw.pdf'))
} else {
    seurat.obj = NormalizeData(seurat.obj)
    seurat.obj = FindVariableFeatures(seurat.obj, selection.method='vst', nfeatures=2000)
    seurat.obj <- ScaleData(seurat.obj, vars.to.regress = "percent.mt")
}

print('test b................')

pdf(paste0("../results/", args$label, "_decide_resolution.pdf"), width=6, height=10)
clustree(seurat.obj, prefix='SCT_snn_res.')
dev.off()

seurat.obj <- FindClusters(seurat.obj, resolution=args$res) ## set resolution

pdf(paste0("../results/", args$label, "_", args$species, "_umap_vs_tsne.pdf"), width=12, height=5)
p1 = DimPlot(object=seurat.obj, reduction='umap', group.by=c("orig.ident", "seurat_clusters"), ncol=1)
p2 = DimPlot(object=seurat.obj, reduction='tsne', group.by=c("orig.ident", "seurat_clusters"), ncol=1)
plots = CombinePlots(plots=list(p1, p2), ncol=2)
print(plots)
dev.off()

if (args$species == 'fish') {
    library(ggpubr)
    pdf(paste0('../results/', args$label, '_markers_', args$species, '.pdf'), width=30, height=10)
    for (m in gsub('_', '-', all_markers_final$tumor_related$tumor_markers)) {
        if (sum(rownames(seurat.obj) %in% m) == 1) {
            p7 = FeaturePlot(seurat.obj,	  reduction='umap', features=m)
            p8 = FeaturePlot(seurat.obj,	  reduction='tsne', features=m)
            plot.markers <- CombinePlots(list(p7, p8), ncol=2)
            plots.all <- ggarrange(plots, plot.markers, heights = c(2, 1), ncol=1)
            print(plots.all)
        }
    }
    dev.off()

    library(gridExtra)
    library(grid)
    pdf(paste0("../results/", args$label, "nontumor_markers_", args$species, ".pdf"), width=8, height=12)
    top3markers <- unlist(all_markers_final$top3_markers)
    top3markers <- top3markers[top3markers %in% rownames(seurat.obj)]
    for (cell.type in 1:length(top3markers)) {
        j = top3markers[cell.type]
        cell.type = names(top3markers)[cell.type]
        p7 = FeaturePlot(seurat.obj,        reduction='umap', features=j)
        p8 = FeaturePlot(seurat.obj,        reduction='tsne', features=j)
        plot.markers <- CombinePlots(list(p7, p8), ncol=2)
        grid.arrange(plots, plot.markers,
                     ncol=1,
                     top=paste0(cell.type, " ", j),
                     heights = c(2, 1),
                     bottom = textGrob(paste0(cell.type, " ", j), gp = gpar(fontface = 2, fontsize = 14),
                                       hjust = 1,
                                       x = 1))
    }
    dev.off()
}

## After reviewing markers plot, run second time to save the tumor cells seurat object
if (length(args$tumor) >= 1) {
    if (args$tumor[1] == -1) {
        print(dim(seurat.obj))
        saveRDS(seurat.obj, paste0('../results/seurat/', args$label, '_seurat_obj_tumors.rds'))
    } else {
        seurat.obj <- subset(seurat.obj, ident=args$tumor)
        print(dim(seurat.obj))
        saveRDS(seurat.obj, paste0('../results/seurat/', args$label, '_seurat_obj_tumors.rds'))
    }
}

