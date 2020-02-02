## use ggplot2 color to velocity plot
## https://github.com/velocyto-team/velocyto.R/issues/16
## devtools::install_github('cole-trapnell-lab/monocle3')

## from https://github.com/satijalab/seurat/files/3418898/Final_Seurat_to_Monocle3_2D_and_3D_190719.2.txt
library(Seurat)
library(monocle3)
library(slingshot)
library(htmlwidgets)
library(SeuratWrappers)
library(infercnv)

Dim = "2D"
input.dir <- getwd()
nPC <- 15
cluster.res <- 0.7
dir.create(file.path('../results/monocle_zebrafish'), showWarnings = FALSE)

### Reading in Seurat object
print("Readingin Seurat objects")

tumorname = 'Tumor24'
##tumorname = 'Tumor21'
input.dir <- '../results/final_seurat_normalize_cluster_obj.RDS'
seurat <- readRDS(file = input.dir)
veloc  = readRDS('../results/final_velocity_obj.RDS')

if (tumorname=='Tumor24') {
    allcells  <- seurat[[1]]
    tumor.vel <- veloc[[1]]
    tumor <- subset(allcells, idents=c(0, 3))
    tumor.cell <- paste0('Tumor24_zebrafish_with_orf_color_v2:', colnames(tumor), 'x')
    tumor.vel <- subset(tumor.vel, cells=tumor.cell)
    tumor.cell = colnames(tumor)[tumor.cell %in% colnames(tumor.vel)]
} else if (tumorname == 'Tumor21') {
    allcells  <- seurat[[2]]
    tumor.vel <- veloc[[2]]
    tumor <- subset(allcells, idents=c(0, 1))
    tumor.cell <- rep("", length(colnames(tumor)))
    tumor.cell[grepl('sort_', colnames(tumor))] = paste0('sort_Tumor21_zebrafish_with_orf_color_v2:', gsub('sort_', '', colnames(tumor)[grepl('sort_', colnames(tumor))], 'x'))
    tumor.cell[grepl('bulk_', colnames(tumor))] = paste0('bulk_Tumor21_2ndlibrary_zebrafish_with_orf_color_v2:', gsub('bulk_', '', colnames(tumor)[grepl('bulk_', colnames(tumor))]), 'x')
    tumor.vel <- subset(tumor.vel, cells=tumor.cell)
    tumor.cell = colnames(tumor)[tumor.cell %in% colnames(tumor.vel)]
}

tumor = subset(tumor, cells=tumor.cell)

### Re-dimension reduction for 3D rendering
if (Dim == "3D") {
  print ("Running UMAP 3D")
  seurat <- RunUMAP(object = seurat, reduction = "pca", dims = 1:nPC, n.components = 3)
  print("Clustering 3D")
  seurat <- FindNeighbors(object=seurat, dims=1:nPC)
  seurat <- FindClusters(object=seurat, resolution=cluster.res)
  seurat[[sprintf("ClusterNames_%.1f_%dPC", cluster.res, nPC)]] <- Idents(object = seurat)
}

### Building the necessary parts for a basic cds
# part one, gene annotations
gene_annotation <- as.data.frame(rownames(tumor.vel@reductions[["pca"]]@feature.loadings), row.names = rownames(tumor.vel@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
# part two, cell information
cell_metadata <- as.data.frame(tumor.vel@assays[["SCT"]]@counts@Dimnames[[2]], row.names = tumor.vel@assays[["SCT"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
# part three, counts sparse matrix
New_matrix <- tumor.vel@assays[["SCT"]]@counts
New_matrix <- New_matrix[rownames(tumor.vel@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix
### Construct the basic cds object
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


### Construct and assign the made up partition, this is used to exclude unrelated cell populations..
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

### Assign the cluster info from Seurat object 
source('functions.R')
tumor <- recluster.withtree(tumor, name=tumorname)
tumor <- FindClusters(tumor, resolution=0.8)

sum(gsub('x', '', gsub('Tumor24_zebrafish_with_orf_color_v2:', '', colnames(tumor.vel))) == colnames(tumor))
Idents(tumor.vel) = Idents(tumor)
tumor.vel$seurat_clusters = tumor$seurat_clusters

tumor.vel.loom = as.loom(tumor.vel, filename='../results/tumor24_seurat_vel_monocle_obj.loom', verbose=T)
tumor.vel.loom$close_all()

list_cluster <- tumor@meta.data$seurat_clusters
names(list_cluster) <- colnames(tumor.vel@assays[["SCT"]]@counts)

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
### Could be a space-holder, but essentially fills out louvain parameters
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

### Assign UMAP coordinate from velocity analysis
## cds_from_seurat@reducedDims@listData[["UMAP"]] <- tumor.vel@reductions[["umap"]]@cell.embeddings[,c('UMAP_2', 'UMAP_1')]
cds_from_seurat@reducedDims@listData[["UMAP"]] <- tumor.vel@reductions[["umap"]]@cell.embeddings

### Assign feature loading for downstream module analysis
cds_from_seurat@preprocess_aux$gene_loadings <- tumor.vel@reductions[["pca"]]@feature.loadings

### Learn graph, this step usually takes a significant period of time for larger samples
print("Learning graph, which can take a while depends on the sample")
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)

print("Plotting clusters")
if (Dim == "2D") {
  pdf(sprintf("%s/%s.with.trajectory.%s.pdf", output.dir, tumorname, Dim), width = 6, height = 3.5)
  clus <- plot_cells(cds_from_seurat, 
                color_cells_by = 'cluster',reduction_method='UMAP',
                label_groups_by_cluster=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE) + ggplot2::xlim(-6, 9) #+ ggplot2::ylim(-6, 5)
  print(clus)
  dev.off()
}

if (Dim == "3D") {
  cluster <- plot_cells_3d(cds_from_seurat, 
                           color_cells_by = 'cluster',
                           label_groups_by_cluster=FALSE,
                           label_leaves=FALSE,
                           label_branch_points=FALSE)
  saveWidget(cluster, file = sprintf("%s/clusters.with.trajectory.%s.html", output.dir, Dim))
}

### Plot cell type info with trajectory
print("Plotting cell type info")
print("Reading in cell type information")
if (tumorname == 'Tumor24') {
    ## lineage.table <- read.table(file = Celltypefile, sep = "", header = T)
    ## lineage.table$cell.barcodes.i. <- paste(lineage.table$cell.barcodes.i., "-1", sep = "")
    ## lineage.table <- lineage.table[, 1:2]
    ## lineage.table$lineage1 <- gsub("_\\d+","", lineage.table$lineage1)
    ## lineage.table$lineage1 <- gsub("TCELLA","T-CELL", lineage.table$lineage1)
    ## lineage.table$lineage1 <- gsub("BCELLA","B-CELL", lineage.table$lineage1)
    ## lineage.table$lineage1 <- gsub("NKA","NK", lineage.table$lineage1)
    ## lineage.table$lineage1 <- gsub("DENDRA","DENDRITIC", lineage.table$lineage1)
    ## lineage.table$lineage1 <- gsub("DENDA","DENDRITIC", lineage.table$lineage1)
    ## lineage.table$lineage1 <- gsub("[0-9]+", "", lineage.table$lineage1)
    ## lineage.table$cell.barcodes.i. <- gsub("-1", "", lineage.table$cell.barcodes.i.)
    ## indices <- match(pData(cds_from_seurat)$barcode, lineage.table$cell.barcodes.i.)
    ## lineage.table.refined <- lineage.table[indices,]
    cluster_names = c("EMT", "EMT", "EMT", "MYC-N", "Hypoxia", "GroundState", "EMT", "DiffMuscle", "Hypoxia", "CellCycle")
    lineage.table.final <- cluster_names[as.numeric(list_cluster)]
    colData(cds_from_seurat)$celltype <- lineage.table.final
}

if (Dim == "2D") {
  pdf(sprintf("%s/%s.celltype.with.trajectory.%s.pdf", output.dir, tumorname, Dim), width = 6, height = 3.5)
  ctype <- plot_cells(cds_from_seurat, 
                color_cells_by = 'celltype',
                label_groups_by_cluster=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE)+ggplot2::xlim(-6, 12) #+ ggplot2::ylim(-6, 5)
  print(ctype)
  dev.off()
}

if (Dim == "3D") {
  celltype <- plot_cells_3d(cds_from_seurat, 
                            color_cells_by = 'celltype',
                            label_groups_by_cluster=FALSE,
                            label_leaves=FALSE,
                            label_branch_points=FALSE)
  saveWidget(celltype, file = sprintf("%s/celltype.with.trajectory.%s.html", output.dir, Dim))
}


### You can also plot out different subclone information, which is not shown here


# Before we order the cells in pseudotime, we need to identify which is the pricipal node/cells
# We can either use their default GUI (which unfortunately is not compatible with 3D), or in the case of NBM, select the HSCs
# For the sake of completeness, we will feed in the cells they are labeled as HSCs, but in other samples for example HSCs, it could be a different story
# plus I am not sure if the GUI will show up on the cluster
# This is using the GUI
## cds_from_seurat <- order_cells(cds_from_seurat)

root_cell_list <- grep("Ground", colData(cds_from_seurat)$celltype)
root_cell_list <- counts(cds_from_seurat)@Dimnames[[2]][root_cell_list]

# Sometimes, in this case, the root cells are dispered so single cell has to be identified to input into the order cells function
cds_from_seurat <- order_cells(cds_from_seurat, root_cells = root_cell_list) # counts(cds_from_seurat)@Dimnames[[2]][1])
## cds_from_seurat <- order_cells(cds_from_seurat, root_pr_nodes = "5") # counts(cds_from_seurat)@Dimnames[[2]][1])

### Now we plot pseudotime
print("Plotting pseudotime")
if (Dim == "2D") {
  pdf(sprintf("%s/%s.pseudotime.with.trajectory.%s.pdf", output.dir, tumorname, Dim), width = 6, height = 3.5)
  ptime <- plot_cells(cds_from_seurat, 
                color_cells_by = 'pseudotime',
                label_groups_by_cluster=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE)+ggplot2::xlim(-6, 12) #+ ggplot2::ylim(-6, 5)
  print(ptime)
  dev.off()
}

## if (Dim == "3D") {
##   pseudotime <- plot_cells_3d(cds_from_seurat, 
##                               color_cells_by = 'pseudotime',
##                               label_groups_by_cluster=FALSE,
##                               label_leaves=FALSE,
##                               label_branch_points=FALSE)
##   saveWidget(pseudotime, file = sprintf("%s/pseudotime.with.trajectory.%s.html", output.dir, Dim))
## }

## ### Now we can subset the cells and do branch expression/pseudotime analysis
## ### This first part is to identify genes that are differentially expressed on the different paths through the trajectory

## # This function is disabled for now

## p <- 0
## if (p > 0) {
  
##   print("calculating most significant genes across the whole trajectory, EXRREMELY time-consuming")
  
##   cds_from_seurat_pr_test_res = graph_test(cds_from_seurat, neighbor_graph="principal_graph", cores=4)
##   pr_deg_ids = row.names(subset(cds_from_seurat_pr_test_res, q_value < 0.05))
  
##   gene_list_1 <- pr_deg_ids[1:4]
  
  
##   ### Plotting individual genes on the map
  
##   print("Plotting individual genes across the whole trajectory")
  
##   if (Dim == "2D") {
##     pdf(sprintf("%s/top4.highly.sig.genes.with.trajectory.%s.pdf", output.dir, Dim), width = 20, height = 20)
##     hgenes <- plot_cells(cds_from_seurat, 
##                   genes = gene_list_1,
##                   show_trajectory_graph=FALSE,
##                   label_cell_groups=FALSE,
##                   label_leaves=FALSE)
##     print(hgenes)
##     dev.off()
##   }
  
##   if (Dim == "3D") {
##     indi_genes <- plot_cells_3d(cds_from_seurat, 
##                                 genes = gene_list_1,
##                                 label_groups_by_cluster=FALSE,
##                                 label_leaves=FALSE,
##                                 label_branch_points=FALSE)
##     saveWidget(indi_genes, file = sprintf("%s/top4.highly.sig.genes.with.trajectory.%s.html", output.dir, Dim))
##   }
  
##   ### Now we can collect the trajectory-variable genes into modules for co-reg elements
  
##   gene_module_df = monocle3:::find_gene_modules(cds_from_seurat[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
  
##   cell_group_df = tibble::tibble(cell=row.names(colData(cds_from_seurat)), cell_group=colData(cds_from_seurat)$cell.type)
##   agg_mat = aggregate_gene_expression(cds_from_seurat, gene_module_df, cell_group_df)
##   row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
##   pheatmap::pheatmap(agg_mat,
##                      scale="column", clustering_method="ward.D2")
  
##   ### Now we can plot these modules
  
##   if (Dim == "2D") {
##     pdf(sprintf("%s/sig.modules.with.trajectory.%s.pdf", output.dir, Dim), width = 20, height = 20)
##     smodules <- plot_cells(cds_from_seurat,
##                   genes=gene_module_df %>% filter(module %in% c(29, 20, 11, 22)),
##                   label_cell_groups=FALSE,
##                   show_trajectory_graph=FALSE)
##     print(smodules)
##     dev.off()
##   }
  
##   if (Dim == "3D") {
##     indi_modules <- plot_cells_3d(cds_from_seurat, 
##                                   genes=gene_module_df %>% filter(module %in% c(29, 20, 11, 22)),
##                                   label_groups_by_cluster=FALSE,
##                                   label_leaves=FALSE,
##                                   label_branch_points=FALSE)
##     saveWidget(indi_modules, file = sprintf("%s/sig.modules.with.trajectory.%s.html", output.dir, Dim))
##   }
  
##   ### Or we can plot them against pseudotime
##   # Say you want to study the expression of certain gene (MPO) in certain lineages (HSC - MEP - CEP - Monocyte)
  
##   Mono_genes = c("MPO")
##   Mono_lineage_cds_from_seurat = cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% Mono_genes,
##                          colData(cds_from_seurat)$celltype %in% c("MONO", "MEP", "CEP", "HSC")]
  
##   print("Plotting genes against pseudotime")
  
##   pdf(sprintf("%s/Genes.against.pseudotime.%s.pdf", output.dir, Dim), width = 20, height = 20)
##   plot_genes_in_pseudotime(Mono_lineage_cds_from_seurat,
##                            color_cells_by="celltype",
##                            min_expr=0.5)
##   dev.off()
  
## }


## ### Now we can subset the cells and do branch expression/pseudotime analysis
## ### This second part is to identify genes that are differentially expressed only within this selected region but not necessarily with the trajectories

## # This function is disabled for now as well

## if (p > 1) {
  
##   # Subsetting the cells
##   # Yet this step required GUI, which I am not sure how it will pan out on the cluster
  
##   cds_from_seurat_subset = choose_cells(cds_from_seurat)
  
  
##   # Identify the significant genes
  
##   print("Subset analysis")
  
##   subset_pr_test_res = graph_test(cds_from_seurat_subset, cores=4)
##   pr_deg_ids = row.names(subset(subset_pr_test_res, q_value < 0.05))
  
  
##   # Group them into modules
##   # Here the resolution parameters and k might vary so feel free to adjust them as well
##   # Sometimes when there are too few genes we might want to 
 
##   gene_module_df = monocle3:::find_gene_modules(cds_from_seurat_subset[pr_deg_ids,], resolution=1e-5, verbose = T, k = 10)
  
  
##   # Plot the scores for the modules
##   cell_group_df = tibble::tibble(cell=row.names(colData(cds_from_seurat_subset)), cell_group=clusters(cds_from_seurat_subset)[colnames(cds_from_seurat_subset)])
##   agg_mat = aggregate_gene_expression(cds_from_seurat_subset, gene_module_df, cell_group_df)
##   module_dendro = hclust(dist(agg_mat))
##   gene_module_df$module <- factor(gene_module_df$module, levels = row.names(agg_mat)[module_dendro$order])
  
  
##   # Plot the modules along the trajectory
##   if (Dim == "2D") {
##     pdf(sprintf("%s/subset.sig.modules.with.trajectory.%s.pdf", output.dir, Dim), width = 20, height = 20)
##     submodules <- plot_cells(cds_from_seurat,
##                   genes=gene_module_df,
##                   label_cell_groups=FALSE,
##                   show_trajectory_graph=FALSE)
##     print(submodules)
##     dev.off()
##   }
  
##   if (Dim == "3D") {
##     subset_indi_modules <- plot_cells_3d(cds_from_seurat, 
##                                          genes=gene_module_df,
##                                          label_groups_by_cluster=FALSE,
##                                          label_leaves=FALSE,
##                                          label_branch_points=FALSE)
##     saveWidget(subset_indi_modules, file = sprintf("%s/subset.sig.modules.with.trajectory.%s.html", output.dir, Dim))
##   }
## }

## ### End of the script
