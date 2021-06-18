library(Seurat)
library(monocle3)
library(SeuratWrappers)
library(patchwork)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')

    parser$add_argument('--seuratobj', dest='data', default='')
    parser$add_argument('--id', dest='id', default='')
    args = parser$parse_args()
    args
}

args = get_args()
if (args$data == '' || args$id == '') {
    cat('empty argument, exit..')
    q()
}

## args$data='/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds'
## args$id = 'MAST139'

args$data = '../results/seurat_sara/20191031_MSK74711_seurat-object.rds'
args$id = 'MSK74711'

seurat.obj <- readRDS(args$data)
obj <- args$id

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
                'gray',
                'gray',
                'gray',
                'gray',
                'gray')
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'UPR',
                "Unique #1",
                "Unique #2",
                "Unique #3",
                "Unique #4",
                "Unique #6")
names(metacolors) <- metalabels

cols = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1, header=T,
                  check.names=F, stringsAsFactors=F)

construct.cds = function(seurat) {
    gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
    colnames(gene_annotation) <- "gene_short_name"
    cell_metadata <- as.data.frame(seurat@assays$RNA@counts@Dimnames[[2]], row.names = seurat@assays$RNA@counts@Dimnames[[2]])
    colnames(cell_metadata) <- "barcode"
    New_matrix <- seurat@assays$RNA@counts
    New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
    expression_matrix <- New_matrix
    cds_from_seurat <- new_cell_data_set(expression_matrix,
                                         cell_metadata = cell_metadata,
                                         gene_metadata = gene_annotation)
    recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
    names(recreate.partition) <- cds_from_seurat@colData@rownames
    recreate.partition <- as.factor(recreate.partition)
    cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
    ## list_cluster <- seurat@meta.data$seurat_clusters
    labels = unlist(cols[obj, ,drop=T])
    list_cluster <- labels[as.vector(seurat.obj@meta.data$RNA_snn_res.0.8)]
    names(list_cluster) <- colnames(seurat)
    cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
    cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
    print(table(list_cluster))
    print('test------')
    cds_from_seurat@reducedDims@listData[["UMAP"]] <- seurat@reductions[["umap"]]@cell.embeddings
    cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings
    print('test------')
    cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)
    cds_from_seurat
}

monoclecds = construct.cds(seurat.obj)

output.dir = '../results/'
tumorname = obj
Dim = '2D'
print('test------')

if (sum(monoclecds@clusters@listData[["UMAP"]][["clusters"]] == 'EMT') > 0) {
    pdf(sprintf("%s/%s.with.monocle.%s.pdf", output.dir, tumorname, Dim), width = 24, height = 5)
} else {
    pdf(sprintf("%s/%s.with.monocle.%s.pdf", output.dir, tumorname, Dim), width = 20, height = 5)
}
clus <- plot_cells(monoclecds,
                   color_cells_by = 'cluster', reduction_method='UMAP',
                   label_groups_by_cluster=F,
                   label_cell_groups=F,
                   group_label_size=0.01,
                   label_leaves=F, show_trajectory_graph=F,
                   label_branch_points=F)+ggplot2::scale_color_manual(values=metacolors)+ggplot2::theme(legend.position = 'right')
root_cell_list = names(monoclecds@clusters@listData[["UMAP"]][["clusters"]])[monoclecds@clusters@listData[["UMAP"]][["clusters"]] == 'Prolif']
monoclecds <- order_cells(monoclecds, root_cells=root_cell_list)
ptime <- plot_cells(monoclecds,
                    color_cells_by = 'pseudotime',
                    label_groups_by_cluster=FALSE,
                    label_leaves=FALSE,
                    label_branch_points=FALSE)+ggplot2::xlim(-6, 12)
pout = pseudotime(monoclecds)
if (sum(monoclecds@clusters@listData[["UMAP"]][["clusters"]] == 'EMT') > 0) {
    root_cell_list = names(monoclecds@clusters@listData[["UMAP"]][["clusters"]])[monoclecds@clusters@listData[["UMAP"]][["clusters"]] == 'EMT']
    monoclecds <- order_cells(monoclecds, root_cells=root_cell_list)
    ptime2 <- plot_cells(monoclecds,
                         color_cells_by = 'pseudotime',
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE, label_roots = F,
                         label_branch_points=FALSE)+ggplot2::xlim(-6, 12)
    pout2 = pseudotime(monoclecds)
}
root_cell_list = names(monoclecds@clusters@listData[["UMAP"]][["clusters"]])[monoclecds@clusters@listData[["UMAP"]][["clusters"]] == 'Ground']
monoclecds <- order_cells(monoclecds, root_cells=root_cell_list)
ptime3 <- plot_cells(monoclecds,
                     color_cells_by = 'pseudotime',
                     label_groups_by_cluster=FALSE,
                     label_leaves=FALSE, label_roots = F,
                     label_branch_points=FALSE)+ggplot2::xlim(-6, 12)
pout3 = pseudotime(monoclecds)
if (sum(monoclecds@clusters@listData[["UMAP"]][["clusters"]] == 'EMT') > 0) {
    print((clus +  ptime + ptime2 + ptime3) + plot_layout(widths = c(1.6, 2.3, 2.3, 2.3)))
    pouts = data.frame(prolif=pout, emt=pout2, ground=pout3, clusters=monoclecds@clusters@listData[["UMAP"]][["clusters"]])
} else {
    print((clus +  ptime + ptime3) + plot_layout(widths = c(1.6, 2.3, 2.3)))
    pouts = data.frame(prolif=pout, ground=pout3, clusters=monoclecds@clusters@listData[["UMAP"]][["clusters"]])
}
dev.off()

print(head(pouts))
write.csv(pouts, file=paste0(args$id, '_pt.csv'))
