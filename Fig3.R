library(Seurat)
library(reticulate)
library(tidyverse)
library(patchwork)

cell1 = list(readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190815_seurat-object_MAST85-1cell.rds'),
             readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST85.rds'))

cell2 = list(readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190815_seurat-object_RH74-10cells.rds'),
             readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_RH74.rds'))

cell3 = list(readRDS('../results/seurat_sara/MAST139_1cells_seurat-object.rds'),
             readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds'))

cell4 = list(readRDS('../results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds'),
             readRDS('../results/seurat_sara/C12SC2_seurat-object.rds'))

integration <- function(x, y, label, method='seurat') {
    x$seurat_clusters = x$RNA_snn_res.0.8
    y$seurat_clusters = y$RNA_snn_res.0.8
    states = unlist(as.vector(annotation[label[1], ]))
    states = states[states!='']
    print(states)
    print(levels(x$seurat_clusters))
    levels(x$seurat_clusters) = states
    ## levels(Idents(x)) = states
    states = unlist(as.vector(annotation[label[2], ]))
    states = states[states!='']
    print(states)
    levels(y$seurat_clusters) = states
    ## levels(Idents(y)) = states
    seurat.pseudo.list = list(x, y)
    for (i in 1:length(seurat.pseudo.list)) {
        print(1)
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
annotation = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1, header=T,
                        check.names=F, stringsAsFactors=F)

integrated.seurat1 = integration(cell1[[1]], cell1[[2]], c('MAST85-1', 'MAST85'))

integrated.seurat2 = integration(cell2[[1]], cell2[[2]], c('RH74-10', 'RH74'))

integrated.seurat3 = integration(cell3[[1]], cell3[[2]], c('MAST139-1', 'MAST139'))

integrated.seurat4 = integration(cell4[[2]], cell4[[1]], c('MSK72117-1', 'MSK72117'))

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
                rgb(241, 250, 100, maxColorValue = 255))
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'UPR')
names(metacolors) <- metalabels

integrated.seurat3@meta.data$orig.ident[integrated.seurat3@meta.data$orig.ident=='MAST139_1cells'] = 'MAST139-1cell'

integrated.seurat4@meta.data$orig.ident[integrated.seurat4@meta.data$orig.ident=='C12SC2'] = 'MSK72117-1cell'

integrated.seurat4@meta.data$orig.ident[integrated.seurat4@meta.data$orig.ident!='MSK72117-1cell'] = 'MSK72117'

write.csv(table(integrated.seurat4$seurat_clusters, integrated.seurat4@meta.data[,1]), file='MSK72117.csv', quote=F)
## table(integrated.seurat3$seurat_clusters, integrated.seurat3@meta.data[,1])

integrated.seurat4$seurat_clusters[grep("Unique", integrated.seurat4$seurat_clusters)] = 'Unique'

## pdf('Fig3D.pdf', width=5.2, height=8.5)
pdf('Fig3D.pdf', width=8, height=9.5)
p1d=DimPlot(integrated.seurat1, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors)+ theme(legend.position='right')
p3d=DimPlot(integrated.seurat3, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
p4d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
## pd = CombinePlots(plots=list(p1d, p2d, p3d), ncol=1)
pd = CombinePlots(plots=list(p1d, p3d, p4d), ncol=1)
print(pd)
dev.off()

pdf("Supp3D_RH74.pdf", width=8, height=5)
p2d=DimPlot(integrated.seurat2, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
print(p2d)
dev.off()

write.csv(table(integrated.seurat2$seurat_clusters, integrated.seurat2@meta.data[,1]), file='RH74_10cells.csv', quote=F)
