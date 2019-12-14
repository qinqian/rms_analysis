library(Seurat)

MAST111 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds')
MAST139 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds')
MAST39 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST39.rds')
RH74 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_RH74.rds')
MAST95 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST95.rds')
MAST85 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST85.rds')
MSK82489 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MSK82489.rds')
MAST35 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST35.rds')

all.obj = Reduce('merge', list(MAST111, MAST139, MAST39, RH74, MAST95, MAST85, MSK82489, MAST35))

all.obj <- NormalizeData(object = all.obj, normalization.method = "LogNormalize", scale.factor = 10000)

all.obj <- FindVariableFeatures(object = all.obj, selection.method = "vst",
                                mean.function = ExpMean, dispersion.function = LogVMR,
                                x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)

all.obj <- ScaleData(object = all.obj, genes.use = rownames(all.obj),
                     vars.to.regress = c("nUMI", "nGene", "percent.mito", "percent.ribo", "fraction.mouse"),
                     model.use = "linear", use.umi = FALSE) 

highvar.genes <- head(VariableFeatures(object = all.obj), 1000)

all.obj <- RunPCA(object = all.obj, features = highvar.genes,
                  npcs = 50, ndims.print = 1:5, nfeatures.print = 1:5,
                  reduction.name = "pca", reduction.key = "PC_", seed.use = 123)

all.obj <- JackStraw(object = all.obj, num.replicate = 100)

all.obj <- ScoreJackStraw(object = all.obj, dims = 1:20)

all.obj <- FindNeighbors(object = all.obj, k.param = 20, dims = 1:20, reduction = "pca")
## all.obj <- FindClusters(object = all.obj, reduction.type = "pca", dims.use = 1:20, 
##                         algorithm = 1, 
##                         resolution = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
##                         random.seed = 123)
                     
all.obj <- RunUMAP(all.obj, reduction.use = "pca", dims = 1:20, reduction.name = "umap")

saveRDS(all.obj, file='all_eight_samples.RDS')

all.obj = readRDS('all_eight_samples.RDS')

pdf('all_eight_QC.pdf', width=60, height=6)
VlnPlot(object = all.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo","fraction.mouse"), 
        pt.size = 0.001, ncol = 6, group.by='orig.ident')
dev.off()

pdf('eight_human_integration.pdf')
DimPlot(all.obj, reduction='umap', group.by = 'orig.ident')
dev.off()

tumor24 = readRDS('../../results/seurat_intersect_velocity/Tumor24_seu.rds')
tumor21 = readRDS('../../results/seurat_intersect_velocity/Tumor21_seu.rds')
tumor22 = readRDS('../../results/seurat_intersect_velocity/Tumor22_seu.rds')

tumor21@meta.data$orig.ident='Tumor21'
tumor22@meta.data$orig.ident='Tumor22'

fish.obj <- Reduce("merge", list(tumor21, tumor22, tumor24))

fish.obj <- NormalizeData(object = fish.obj, normalization.method = "LogNormalize", scale.factor = 10000)

fish.obj <- FindVariableFeatures(object = fish.obj, selection.method = "vst",
                                mean.function = ExpMean, dispersion.function = LogVMR,
                                x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)

fish.obj <- ScaleData(object = fish.obj, genes.use = rownames(fish.obj),
                     vars.to.regress = c("percent.mt"),
                     model.use = "linear", use.umi = FALSE)

highvar.genes <- head(VariableFeatures(object = fish.obj), 1000)

fish.obj <- RunPCA(object = fish.obj, features = highvar.genes,
                  npcs = 50, ndims.print = 1:5, nfeatures.print = 1:5,
                  reduction.name = "pca", reduction.key = "PC_", seed.use = 123)

fish.obj <- JackStraw(object = fish.obj, num.replicate = 100)

fish.obj <- ScoreJackStraw(object = fish.obj, dims = 1:20)

fish.obj <- FindNeighbors(object = fish.obj, k.param = 20, dims = 1:20, reduction = "pca")

## fish.obj <- FindClusters(object = fish.obj, reduction.type = "pca", dims.use = 1:20, 
##                         algorithm = 1, 
##                         resolution = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
##                         random.seed = 123)
fish.obj <- RunUMAP(fish.obj, reduction.use = "pca", dims = 1:20, reduction.name = "umap")

saveRDS(fish.obj, file='all_three_zebrafish.RDS')

pdf('three_zebrafish_integration.pdf')
DimPlot(fish.obj, reduction='umap', group.by = 'orig.ident')
dev.off()
