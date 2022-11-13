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
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'UPR', "Neural")
names(metacolors) <- metalabels

library(Seurat)
library(reticulate)
library(patchwork)
library(ggplot2)
library(patchwork)

res <- readRDS('results/seurat_sara/MAST-39-Resistant_seurat-object.rds')
sen <- readRDS('results/seurat_sara/MAST-39-Sensitive_seurat-object.rds')
res$seurat_clusters = res$RNA_snn_res.0.8
sen$seurat_clusters = sen$RNA_snn_res.0.8

cell4 = list(res, sen)

## Integration by CCA analysis
integration <- function(x, y, label, method='seurat') {
    Idents(x) = x$seurat_clusters2 = x$seurat_clusters
    states = unlist(as.vector(annotation[label[1], ]))
    states = states[states!='']
    states = states[!is.na(states)]
    levels(x$seurat_clusters) = states
    print('-------')
    ## levels(Idents(x)) = states
    Idents(y) = y$seurat_clusters2 = y$seurat_clusters
    states = unlist(as.vector(annotation[label[2], ]))
    states = na.omit(states[states!=''])
    print(length(states))
    print(levels(y$seurat_clusters))
    levels(y$seurat_clusters) = states
    ## levels(Idents(y)) = states
    seurat.pseudo.list = list(x, y)
    for (i in 1:length(seurat.pseudo.list)) {
        print(i)
        seurat.pseudo.list[[i]] <- NormalizeData(seurat.pseudo.list[[i]], verbose = FALSE)
        seurat.pseudo.list[[i]] <- FindVariableFeatures(seurat.pseudo.list[[i]],
                                                        selection.method = "vst", 
                                                        nfeatures = 2000, verbose = FALSE)
    }
    names(seurat.pseudo.list) <- label
    if (method == 'seurat') {
        anchors <- FindIntegrationAnchors(object.list = seurat.pseudo.list, dims = 1:30)
        integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
        DefaultAssay(integrated) <- "integrated"
        integrated <- ScaleData(integrated, verbose = FALSE)
        integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
    } else if (method == 'conos') {
        seurat.pseudo.list <- lapply(seurat.pseudo.list, function(x) {
            ScaleData(x) %>% RunPCA(verbose=F)
        })
        pseudo.con <- Conos$new(seurat.pseudo.list)
        pseudo.con$buildGraph(k=15, k.self=5, space="PCA", ncomps=30, n.odgenes=2000, matching.method='mNN',
                              metric = 'angular', score.component.variance=T, verbose=T)
        pseudo.con$findCommunities()
        pseudo.con$embedGraph()
        integrated = as.Seurat(pseudo.con)
    }
    integrated
}

set.seed(100)
annotation = read.delim('clusters.txt', sep='\t', row.names=1, header=T,
                        check.names=F, stringsAsFactors=F)

integrated.seurat4 = integration(cell4[[1]], cell4[[2]], c('MAST39-Resistant', 'MAST39-Sensitive'))

pdf('ALL_PDX_MAST39_batchcorrected_UMAP.pdf', width=16, height=10)
p4d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='seurat_clusters2', label=T) + theme(legend.position='right')
p5d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
print(p4d / p5d)
dev.off()

res <- readRDS('results/seurat_sara/MAST-139-Resistant_seurat-object.rds')
sen <- readRDS('results/seurat_sara/MAST-139-Sensitive_seurat-object.rds')
res$seurat_clusters = res$RNA_snn_res.0.8
sen$seurat_clusters = sen$RNA_snn_res.0.8

cell4 = list(res, sen)
integrated.seurat5 = integration(cell4[[1]], cell4[[2]], c('MAST139-Resistant', 'MAST139-Sensitive'))

pdf('ALL_PDX_MAST139_batchcorrected_UMAP.pdf', width=16, height=10)
p4d=DimPlot(integrated.seurat5, reduction='umap', split.by='orig.ident', group.by='seurat_clusters2', label=T) + theme(legend.position='right')
p5d=DimPlot(integrated.seurat5, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
print(p4d / p5d)
dev.off()
