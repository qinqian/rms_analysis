translate.dim.code=function(reduction.use) {
  return.code="PC"
  if (reduction.use=="ica") return.code="IC"
  if (reduction.use=="tsne") return.code="tSNE_"
  if (reduction.use=="mds") return.code="MDS"
  return(return.code)
}

GeoFeaturePlot <-  function(object, features, reduction.use = "tsne", dim.1 = 1,
                            dim.2 = 2, pt.size = 1, cells.use = NULL, use.imputed = FALSE,
                            cols.use = c("grey","blue"), pch.use = 16) {
  dim.code = translate.dim.code(reduction.use)
  dim.codes = paste(dim.code, c(dim.1, dim.2), sep = "")
  data.plot = FetchData(object, dim.codes)
  x1 = paste(dim.code, dim.1, sep = "")
  x2 = paste(dim.code, dim.2, sep = "")
  data.plot$x = data.plot[, x1]
  data.plot$y = data.plot[, x2]
  data.plot$pt.size = pt.size
  data.use = data.frame(t(FetchData(object, features, cells = cells.use))) 
  geomeans = colMeans(data.use)
  if (length(cols.use) == 1) {
    brewer.gran <- brewer.pal.info[cols.use, ]$maxcolors
  } else {
    brewer.gran <- length(cols.use)
  }
  data.plot$geomean = geomeans
  data.cut = as.numeric(cut(as.numeric(geomeans), breaks = 2))
  data.plot$cut = data.cut
  data.plot$col = as.factor(data.cut)
  data.plot$ident = as.factor(Idents(seurat.obj)) 
  data_plot_cell <- data.plot %>% tibble::rownames_to_column(var = "cell")
  metadata_cell <- object@meta.data %>% tibble::rownames_to_column(var="cell")
  data.plot <- data_plot_cell %>% left_join(metadata_cell, by = "cell")
  return(data.plot)
}

correct <- function(x) { return (1e6*x/sum(x))}

get_upregulated_genes_in_clusters <- function(seurat_object, forgeround_clusters, control_clusters, logfc_threshold=1, fg_expression_threshold=10)
{
  cluster_cells_fg <- WhichCells(seurat_object, idents = forgeround_clusters) 
  control_all <- WhichCells(seurat_object, idents = control_clusters) 
  return (get_upregulated_genes_in_cells(seurat_object, cluster_cells_fg, control_all, logfc_threshold, fg_expression_threshold))
}

get_upregulated_genes_in_cells <- function(seurat_object, foreground_cells, control_cells, logfc_threshold=1, fg_expression_threshold=10)
{
  raw_data_cluster_cells_fg <- seurat_object[["RNA"]]@counts[rownames(seurat_object[["RNA"]]@data), foreground_cells] 
  raw_data_control_all <- seurat_object[["RNA"]]@counts[rownames(seurat_object[["RNA"]]@data), control_cells] 
  raw_data_cluster_fg_gene_counts <- rowSums(as.matrix(raw_data_cluster_cells_fg)) 
  raw_data_all_control_gene_counts <- rowSums(as.matrix(raw_data_control_all)) 
  all_cells <- data.frame("cluster_cells_fg"=raw_data_cluster_fg_gene_counts, "control_all"=raw_data_all_control_gene_counts)
  all_cells_corrected <- as.data.frame(apply(all_cells, 2, correct))
  all_cells_corrected$logfc_all_cluster_fg_vs_control <- log2(all_cells_corrected$cluster_cells_fg+0.1) - log2(all_cells_corrected$control_all+0.1)
  all_cells_corrected$fg_fraction <- apply(seurat_object[["RNA"]]@data[,foreground_cells],1,function(x)return(length(x[x>1])/length(x)))
  all_cells_corrected$bg_fraction <- apply(seurat_object[["RNA"]]@data[,control_cells],1,function(x)return(length(x[x>1])/length(x)))
  all_cells_corrected$diff_in_enrichment <- all_cells_corrected$fg_fraction - all_cells_corrected$bg_fraction
  all_cells_corrected$fg_cells_expressing <- apply(seurat_object[["RNA"]]@data[,foreground_cells],1,function(x)return(length(x[x>1])))
  all_cells_corrected$bg_cells_expressing <- apply(seurat_object[["RNA"]]@data[,control_cells],1,function(x)return(length(x[x>1])))
  all_cells_corrected$fg_cells_notexpressing <- length(foreground_cells) - all_cells_corrected$fg_cells_expressing
  all_cells_corrected$bg_cells_notexpressing <- length(control_cells) - all_cells_corrected$bg_cells_expressing
  all_cells_corrected$fisher_exact_p <- foreach(i = 1:nrow(all_cells_corrected), .combine = c) %do%
    {
      m <- matrix(c(all_cells_corrected$fg_cells_expressing[i],
                    all_cells_corrected$bg_cells_expressing[i],
                    all_cells_corrected$fg_cells_notexpressing[i],
                    all_cells_corrected$bg_cells_notexpressing[i]), nrow=2, byrow = TRUE)
      fisher.test(m)$p.value
    }
  all_cells_corrected$p.adjusted <- p.adjust(all_cells_corrected$fisher_exact_p, method = "BH")
  print(head(rownames(all_cells_corrected)[order(all_cells_corrected$logfc_all_cluster_fg_vs_control, decreasing = T)], 20))
  marker_gene_rows <- all_cells_corrected[all_cells_corrected$logfc_all_cluster_fg_vs_control > logfc_threshold & all_cells_corrected$cluster_cells_fg > fg_expression_threshold,]
  marker_genes <- marker_gene_rows[order(marker_gene_rows$diff_in_enrichment, decreasing = TRUE),] %>%
    tibble::rownames_to_column(var="genename") 
  return (marker_genes)
}

