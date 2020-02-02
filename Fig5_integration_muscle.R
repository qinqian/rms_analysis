library(Seurat)
library(reticulate)
library(tidyverse)
library(patchwork)
library(readxl)

cells.pre = list(readRDS('../results/seurat_sara/MAST139_Muscle_plus_seurat-object.rds'),
                 readRDS('../results/seurat_sara/MAST139_Muscle_minus_seurat-object.rds'))

cells = list(CreateSeuratObject(Read10X('MAST139_Muscle_plus/outs/filtered_feature_bc_matrix'), project='Muscle+'),
             CreateSeuratObject(Read10X('MAST139_Muscle_minus/outs/filtered_feature_bc_matrix'), project='Muscle-'))

pdf('mast139_muscle_biomarkers.pdf', width=18, height=10)
p1 = DimPlot(cells.pre[[1]], reduction = "umap",
             group.by = "RNA_snn_res.0.8")
p2 = FeaturePlot(object = cells.pre[[1]], reduction = "umap",
                 cols = c("lightgrey", "red"), order = TRUE,
                 features = c("TSPAN33", "LRRN1"))
p3 = DimPlot(cells.pre[[2]], reduction = "umap",
             group.by = "RNA_snn_res.0.8")
p4 = FeaturePlot(object = cells.pre[[2]], reduction = "umap",
                 cols = c("lightgrey", "red"), order = TRUE,
                 features = c("TSPAN33", "LRRN1"))
(p1 + p2) / (p3 + p4)
dev.off()

integration <- function(x, y, label, method='seurat') {
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
integrated.seurat = integration(cells[[1]], cells[[2]], c('Muscle+', 'Muscle-'))

DefaultAssay(integrated.seurat) <- "RNA"
p1 = FeaturePlot(integrated.seurat, features = c("TSPAN33", "LRRN1"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))
## p2 <- VlnPlot(integrated.seurat, features = c("TSPAN33", "LRRN1"), split.by = "orig.ident", group.by = "orig.ident",
##     pt.size = 0, combine = FALSE)
p1
## ggsave('test.pdf', width=12)
ggsave('test2.pdf', width=12, height=12)

levels(cells.pre[[1]]@meta.data$RNA_snn_res.0.8) = c("EMT", "MUSCLE", "G1S", "MUSCLE", "G2M")
levels(cells.pre[[2]]@meta.data$RNA_snn_res.0.8) = c("Hypoxia", "G1S", "MUSCLE", "EMT", "EMT", "GROUND", "G2M", "MUSCLE", "EMT", "INTERFERON")

pdf('mast139_muscle_biomarkers_labelled.pdf', width=18, height=10)
p1 = DimPlot(cells.pre[[1]], reduction = "umap",
             group.by = "RNA_snn_res.0.8")
p2 = FeaturePlot(object = cells.pre[[1]], reduction = "umap",
                 cols = c("lightgrey", "red"), order = TRUE,
                 features = c("TSPAN33", "LRRN1"))
p3 = DimPlot(cells.pre[[2]], reduction = "umap",
             group.by = "RNA_snn_res.0.8")
p4 = FeaturePlot(object = cells.pre[[2]], reduction = "umap",
                 cols = c("lightgrey", "red"), order = TRUE,
                 features = c("TSPAN33", "LRRN1"))
(p1 + p2) / (p3 + p4)
dev.off()

set.seed(100)
integrated.seurat2 = integration(cells.pre[[1]], cells.pre[[2]], c('Muscle+', 'Muscle-'))

DefaultAssay(integrated.seurat2) <- "RNA"
p1 = DimPlot(integrated.seurat2, group.by='RNA_snn_res.0.8')
p2 = DimPlot(integrated.seurat2, group.by='orig.ident') + NoLegend()
p3 = FeaturePlot(integrated.seurat2, features = c("TSPAN33", "LRRN1"), split.by = "orig.ident", max.cutoff = 3,
            cols = c("grey", "red"))+NoLegend()
p4 <- VlnPlot(integrated.seurat2, features = c("TSPAN33", "LRRN1"), split.by = "orig.ident", group.by = "RNA_snn_res.0.8",
              pt.size = 0, combine = FALSE)
(((p2 + p1 + plot_layout(ncol=2)) / p3) + p4) + plot_layout(heights=c(3, 3), ncol=2)
ggsave('test2.pdf', width=20, height=12)
