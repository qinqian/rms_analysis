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
                '#00FFFD')
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'ER stress', 'Neural')
names(metacolors) <- metalabels

## library(readxl)
library(Seurat)
#library(pheatmap)
library(gplots)
library(RColorBrewer)
## library(tidyverse)
## library(argparse)
library(patchwork)
library(gridtext)
library(ggtext)
library(ggplot2)
library(glue)

## get_args = function(x) {
##     require(argparse)
##     parser = ArgumentParser(description='seurat normalization')
##     parser$add_argument('--seuratobj', dest='seurat', metavar='N', type="character", nargs='+')
##     parser$add_argument('--genes', dest='label', metavar='N', type="character", nargs="+")
##     args = parser$parse_args()
##     args
## }

suppfig2 = readRDS('~/langenau/projects/01_sc_rms/results/seurat_sara/20696_seurat-object.rds')
erms.topsign = c("ITM2A", "FABP5", "EMILIN1", "RRBP1", "CCND1", "EMP3", "EIF4EBP1", "TSTA3", "ADA", "HMGA2", "MYOD1", "DES", "MYC")
## arms.topsign = c("HSPB2", "MYL4", "PIPOX", "TNNT2", "MYOG", "ENO3", "NRN1", "GYPC", "TSPAN3", "TFF3")

suppfig2@meta.data$FN.core = colMeans(as.matrix(suppfig2[erms.topsign[1:10], ]$RNA@data))

pdf("suppfig2_markers.pdf", width=24, height=11.5)
## FeaturePlot(suppfig2, arms.topsign, ncol=5)
FeaturePlot(suppfig2, erms.topsign, ncol=5)
dev.off()

pdf("suppfig2_markers_FNcore.pdf", width=5, height=4)
## FeaturePlot(suppfig2, arms.topsign, ncol=5)
FeaturePlot(suppfig2, features='FN.core')
dev.off()


annotation = read.delim('~/langenau/projects/01_sc_rms/final_annotations/Final_clusters_NG.txt', sep='\t', row.names=1, header=T, check.names=F, stringsAsFactors=F)

pdxs = Sys.glob('~/langenau/projects/01_sc_rms/data/seurat_obj/*rds')[1:10]

pdxs = c(pdxs, '~/langenau/projects/01_sc_rms/results/seurat_sara/20191031_MSK74711_seurat-object.rds',
         '~/langenau/projects/01_sc_rms/results/seurat_sara/MAST118_seurat-object.rds',
         '~/langenau/projects/01_sc_rms/results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds',
         '~/langenau/projects/01_sc_rms/results/seurat_sara/MAST139_1cells_seurat-object.rds')

mast95 = readRDS(pdxs[6])
states = unlist(as.vector(annotation["MAST95", ]))
states = states[states!='']
print(states)

mast95$seurat_clusters = mast95$RNA_snn_res.0.8
levels(mast95$RNA_snn_res.0.8) = states
Idents(mast95) = mast95$seurat_clusters

pdf("SuupFig5B_MAST95.pdf", width=5, height=4)
DimPlot(mast95, cols=metacolors, group.by='RNA_snn_res.0.8')
dev.off()

mast111 = readRDS(pdxs[1])
states = unlist(as.vector(annotation["MAST111", ]))
states = states[states!='']
print(states)

mast111$seurat_clusters = mast111$RNA_snn_res.0.8
levels(mast111$RNA_snn_res.0.8) = states
Idents(mast111) = mast111$seurat_clusters

pdf("Fig1_MAST111.pdf", width=5, height=4)
DimPlot(mast111, cols=metacolors, group.by='RNA_snn_res.0.8')
dev.off()

pdxs = pdxs[c(grep("MAST139", pdxs)[1], grep("MSK74711", pdxs)[1])]

## pdxs = pdxs[c(1, 2, 11, 12)]
pdxs.objs = lapply(pdxs, readRDS)

labels = unlist(lapply(pdxs.objs, function(x) {
    levels(x$orig.ident[1])
}))

## labels[9] = 'MAST85-1'
## labels[11] = 'MSK74711'
## labels[13] = 'MSK72117'
## labels[14] = 'MAST139-1'
## labels[10] = 'RH74-10'

labels[2] = 'MSK74711'

results = list()
for (i in seq_along(labels)) {
    states = unlist(as.vector(annotation[labels[i], ]))
    states = states[states!='']
    print(states)
    pdxs.objs[[i]]$seurat_clusters = pdxs.objs[[i]]$RNA_snn_res.0.8
    levels(pdxs.objs[[i]]$RNA_snn_res.0.8) = states
    Idents(pdxs.objs[[i]]) = pdxs.objs[[i]]$seurat_clusters
    results[[labels[i]]] = table(pdxs.objs[[i]]$RNA_snn_res.0.8)
}

## ALDH1 not exist (use ALDH1B1, ALDH1L1), IBSP not exist
## see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5549701/
## markers1 = c("ALDH1B1", "ALDH1L1", "NANOG", "GLI1", "PROM1", "HHIP", "SOX2",
##              "CXCR4", "POU5F1", "MYC", "PAX3", "YAP1", "NOTCH1", "NOTCH3", "VANGL2",
##              "HEY1", "NES", "ABCG2", "BGLAP", "IBSP")
## markers2 = c("OGN", "DCN", "COL3A1", "VIM", "COL1A1", "SPARC")

## Selected From Dave and Claudia
## markers1 = c("PROM1", "NANOG", "POU5F1", "SOX2", "GLI1", "HHIP", "NOTCH1", "CXCR4", "VANGL2")
markers1 = c("THY1", "CD44", "LRRN1", "TSPAN33", "NDRG1", "EGFR", "TNNT3", "MX1", "MKI67", "CHODL")

#pdf("sciencereviewer_comments.pdf", width=16.5, height=5.5)
pdf("reviewer_comments2.pdf", width=19.5, height=5.5)
for (i in seq_along(labels)) {
    plots = list()
    plots[['clusters']]=DimPlot(pdxs.objs[[i]], group.by = 'RNA_snn_res.0.8', cols=metacolors)+ggtitle(labels[i])
    for (d in markers1) {
        if (d %in% rownames(pdxs.objs[[i]])) {
            plots[[d]]=FeaturePlot(pdxs.objs[[i]], d)+labs(title = glue("*{d}*"))+theme(plot.title = element_markdown())
        } else {
            plots[[d]]=ggplot() + theme_void()
        }
    }
    ## for (w in markers2[markers2 %in% rownames(pdxs.objs[[i]])]) {
    ##     plots[[w]]=FeaturePlot(pdxs.objs[[i]], w)
    ## }
    print(wrap_plots(plots, ncol=6))
}
dev.off()

## bone_markers = c("OGN", "MGP", "SPARC")
## pdf("bone_markers.pdf", width=16, height=3)
## for (i in seq_along(labels)) {
##     plots = list()
##     plots[['clusters']]=DimPlot(pdxs.objs[[i]], group.by = 'RNA_snn_res.0.8', cols=metacolors)+ggtitle(labels[i])
##     for (d in bone_markers[bone_markers %in% rownames(pdxs.objs[[i]])]) {
##         plots[[d]]=FeaturePlot(pdxs.objs[[i]], d, min.cutoff=0.1, repel=T)
##     }
##     print(wrap_plots(plots, ncol=4))
## }
## dev.off()
