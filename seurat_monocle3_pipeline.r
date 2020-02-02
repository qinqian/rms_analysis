library(Seurat)
library(monocle3)
library(SeuratWrappers)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')

    parser$add_argument('--dir10x', dest='data', default='')
    parser$add_argument('--velloom', dest='vel', default='')
    args = parser$parse_args()
    args
}

args = get_args()
if (args$data == '' || args$vel == '') {
    cat('empty argument, exit..')
    q()
}

## intersect seurat object and velocity object
## use velocity cell barcodes for create cell dataset
seurat.obj <- readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds')
seurat.obj = FindClusters(seurat.obj, resolution=0.8)

velocity  <- as.Seurat(ReadVelocity(file = "/data/langenau/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/20190418_MAST139_5Kcells_hg19/velocyto/20190418_MAST139_5Kcells_hg19.loom"))
velocity[["percent.mt"]] <- PercentageFeatureSet(velocity, pattern='^mt-')

velocity <- subset(velocity,
                      subset=nFeature_spliced > 1000 & nFeature_spliced <4000 & percent.mt<10)
source('functions.R')

velocity <- process_standard(velocity, assaytype='spliced', output='../results/MAST139_velocity_jackstraw.pdf', norm=F)

tumor.cells = paste0("20190418_MAST139_5Kcells_hg19:", colnames(seurat.obj), 'x')

velocity = subset(velocity, cells=tumor.cells)
tumor.cells = colnames(seurat.obj)[tumor.cells %in% colnames(velocity)]

seurat.obj = subset(seurat.obj, cells=tumor.cells)

construct.cds = function(vel, seurat) {
    gene_annotation <- as.data.frame(rownames(vel@reductions[["pca"]]@feature.loadings), row.names = rownames(vel@reductions[["pca"]]@feature.loadings))
    colnames(gene_annotation) <- "gene_short_name"

    cell_metadata <- as.data.frame(vel@assays[["SCT"]]@counts@Dimnames[[2]], row.names = vel@assays[["SCT"]]@counts@Dimnames[[2]])
    colnames(cell_metadata) <- "barcode"

    New_matrix <- vel@assays[["SCT"]]@counts
    New_matrix <- New_matrix[rownames(vel@reductions[["pca"]]@feature.loadings), ]
    expression_matrix <- New_matrix

    cds_from_seurat <- new_cell_data_set(expression_matrix,
                                         cell_metadata = cell_metadata,
                                         gene_metadata = gene_annotation)

    recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
    names(recreate.partition) <- cds_from_seurat@colData@rownames
    recreate.partition <- as.factor(recreate.partition)
    cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

    list_cluster <- seurat@meta.data$seurat_clusters
    names(list_cluster) <- colnames(vel@assays[["SCT"]])

    cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
    cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

    cds_from_seurat@reducedDims@listData[["UMAP"]] <- vel@reductions[["umap"]]@cell.embeddings
    cds_from_seurat@preprocess_aux$gene_loadings <- vel@reductions[["pca"]]@feature.loadings

    cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)

    cds_from_seurat
}

source('functions.R')
velocity = recluster.withtree(velocity, name='MAST139')
## velocity = FindClusters(object=velocity, resolution=0.8)

monoclecds = construct.cds(velocity, seurat.obj)
## root_cell_list <- grep("Ground", colData(cds_from_seurat)$celltype)
## root_cell_list <- counts(cds_from_seurat)@Dimnames[[2]][root_cell_list]
root_cell_list = names(monoclecds@clusters@listData[["UMAP"]][["clusters"]])[monoclecds@clusters@listData[["UMAP"]][["clusters"]] == 4]

monoclecds <- order_cells(monoclecds, root_cells=root_cell_list[503])

output.dir = '../results/'
tumorname = 'MAST139'
Dim = '2D'
pdf(sprintf("%s/%s.with.trajectory.%s.pdf", output.dir, tumorname, Dim), width = 9, height = 3.5)
clus <- plot_cells(monoclecds,
                   color_cells_by = 'cluster', reduction_method='UMAP',
                   label_groups_by_cluster=T,
                   label_leaves=T,
                   label_branch_points=T)
ptime <- plot_cells(monoclecds,
                    color_cells_by = 'pseudotime',
                    label_groups_by_cluster=FALSE,
                    label_leaves=FALSE,
                    label_branch_points=FALSE)+ggplot2::xlim(-6, 12)
CombinePlots(plots=list(clus, ptime))
dev.off()
