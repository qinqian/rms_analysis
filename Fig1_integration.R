library(Seurat)

MAST111 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds')
MAST139 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds')
MAST39 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST39.rds')
RH74 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_RH74.rds')
MAST95 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST95.rds')
MAST85 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST85.rds')
MSK82489 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MSK82489.rds')
MAST35 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST35.rds')

## new added samples
MSK72117_10cells = readRDS('/data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds')
MAST118 = readRDS('/data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/MAST118_seurat-object.rds')
MSK74711 = readRDS('/data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/20191031_MSK74711_seurat-object.rds')

MSK72117_10cells@meta.data$orig.ident = gsub('20191031_', '', MSK72117_10cells@meta.data$orig.ident)
MSK74711@meta.data$orig.ident = gsub('20191031_', '', MSK74711@meta.data$orig.ident)

all.obj = Reduce('merge', list(MAST111, MAST139, MAST39, RH74, MAST95, MAST85, MSK82489, MAST35,
                               MSK72117_10cells, MAST118, MSK74711))

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

## all.obj <- JackStraw(object = all.obj, num.replicate = 100)
## all.obj <- ScoreJackStraw(object = all.obj, dims = 1:20)

all.obj <- FindNeighbors(object = all.obj, k.param = 20, dims = 1:20, reduction = "pca")

## all.obj <- FindClusters(object = all.obj, reduction.type = "pca", dims.use = 1:20, 
##                         algorithm = 1, 
##                         resolution = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
##                         random.seed = 123)
all.obj <- RunUMAP(all.obj, reduction.use = "pca", dims = 1:20, reduction.name = "umap")

saveRDS(all.obj, file='../data/seurat_obj/all_eight_samples.RDS')

pdf('../figures/all_eight_QC.pdf', width=60, height=6)
VlnPlot(object = all.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo","fraction.mouse"),
        pt.size = 0.001, ncol = 6, group.by='orig.ident')
dev.off()

pdf('../figures/all_mouse_pdx_integration.pdf', width=10, height=8)
DimPlot(all.obj, reduction='umap', group.by = 'orig.ident')
dev.off()

## tumor24 = readRDS('../results/seurat_intersect_velocity/Tumor24_seu.rds')
## tumor21 = readRDS('../results/seurat_intersect_velocity/Tumor21_seu.rds')
## tumor22 = readRDS('../results/seurat_intersect_velocity/Tumor22_seu.rds')

tumor24 = readRDS('../results/seurat/Tumor24_unfilter_seurat_obj_tumors.rds')
tumor21 = readRDS('../results/seurat/Tumor21_unfilter_seurat_obj_tumors.rds')
tumor22 = readRDS('../results/seurat/Tumor22_unfilter_seurat_obj_tumors.rds')

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
## fish.obj <- JackStraw(object = fish.obj, num.replicate = 100)
## fish.obj <- ScoreJackStraw(object = fish.obj, dims = 1:20)
fish.obj <- FindNeighbors(object = fish.obj, k.param = 20, dims = 1:20, reduction = "pca")

## fish.obj <- FindClusters(object = fish.obj, reduction.type = "pca", dims.use = 1:20, 
##                         algorithm = 1, 
##                         resolution = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
##                         random.seed = 123)
fish.obj <- RunUMAP(fish.obj, reduction.use = "pca", dims = 1:20, reduction.name = "umap")

saveRDS(fish.obj, file='../data/seurat_obj/all_three_zebrafish.rds')

fish.obj@meta.data$orig.ident = gsub('_unfilter', '', fish.obj@meta.data$orig.ident)

pdf('../figures/three_zebrafish_integration.pdf', width=10, height=8)
DimPlot(fish.obj, reduction='umap', group.by = 'orig.ident')
dev.off()


primary.obj <- Reduce("merge", list(
                                readRDS('../results/seurat_sara/20082_hg19_premrna_seurat-object.rds'),
                                readRDS('../results/seurat_sara/21202_hg19_premrna_seurat-object.rds'),
                                readRDS('../results/seurat_sara/20696_seurat-object.rds'),
                                readRDS('../results/seurat_sara/29806_hg19_premrna_seurat-object.rds'))
                      )

primary.obj@meta.data$orig.ident <- gsub('_hg19_premrna', '', primary.obj@meta.data$orig.ident)

primary.obj <- NormalizeData(object = primary.obj, normalization.method = "LogNormalize", scale.factor = 10000)
primary.obj <- FindVariableFeatures(object = primary.obj, selection.method = "vst",
                                    mean.function = ExpMean, dispersion.function = LogVMR,
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)

primary.obj <- ScaleData(object = primary.obj, genes.use = rownames(primary.obj),
                         vars.to.regress = c("percent.mito"),
                         model.use = "linear", use.umi = FALSE)

highvar.genes <- head(VariableFeatures(object = primary.obj), 1000)

primary.obj <- RunPCA(object = primary.obj, features = highvar.genes,
                      npcs = 50, ndims.print = 1:5, nfeatures.print = 1:5,
                      reduction.name = "pca", reduction.key = "PC_", seed.use = 123)

## primary.obj <- JackStraw(object = primary.obj, num.replicate = 100)
## primary.obj <- ScoreJackStraw(object = primary.obj, dims = 1:20)

primary.obj <- FindNeighbors(object = primary.obj, k.param = 20, dims = 1:20, reduction = "pca")

## primary.obj <- FindClusters(object = primary.obj, reduction.type = "pca", dims.use = 1:20, 
##                         algorithm = 1, 
##                         resolution = c(0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
##                         random.seed = 123)

## primary.obj <- RunUMAP(primary.obj, reduction.use = "pca", dims = 1:15, reduction.name = "umap", seed.use=123)
primary.obj <- RunUMAP(primary.obj, reduction.use = "pca", dims = 1:15, reduction.name = "umap", seed.use=42)

saveRDS(primary.obj, file='../data/seurat_obj/all_primary_patient.RDS')

meta.copy = primary.obj@meta.data

library(glue)
for (i in unique(primary.obj@meta.data$orig.ident)) {
    print(i)
    meta.data = meta.copy
    meta.data[primary.obj@meta.data$orig.ident!=i, 'RNA_snn_res.0.8'] = -1
    primary.obj@meta.data = meta.data
    pdf(glue('../figures/primary_patient_integration_label_{i}.pdf'), width=16, height=8)
    p1=DimPlot(primary.obj, reduction='umap', group.by = 'orig.ident')
    p2=DimPlot(primary.obj, reduction='umap', group.by = 'RNA_snn_res.0.8', label=T, legend=F)
    print(CombinePlots(plots=list(p1, p2)))
    dev.off()
}
