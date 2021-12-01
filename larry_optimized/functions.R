library(ggplot2)
#human_ortholog = read.table('~/alvin_singlecell/01_rms_projects//01_fish/data/ortholog_mapping/Beagle_fish_human_all_genes.txt', header=T, sep='\t', stringsAsFactors=F)

process_standard <- function(obj, output, assaytype='RNA', regress_mt=T, norm=T) {
    if (!norm) {
    #if (F) { 
        ## input should always be normalized
        ## SCTransform is incompatible with JackStraw after Seurat 3.1.5
        obj <- SCTransform(obj, assay=assaytype, vars.to.regress = "percent.mt")
        ## Use the Seurat old normalization
        ## obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000)
        ## obj <- FindVariableFeatures(object = obj, selection.method = "vst",
        ##                             mean.function = ExpMean, dispersion.function = LogVMR,
        ##                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
        ## obj <- ScaleData(object = obj, genes.use = rownames(obj), #vars.to.regress = c("nUMI", "nGene", "percent.mito"),
        ##                  model.use = "linear", use.umi = FALSE) 
    }
    obj <- RunPCA(object=obj)
    # this cannot be runned on the SCT... after 3.1.5 Seurat
    #obj <- JackStraw(obj, num.replicate = 100)
    #obj <- ScoreJackStraw(obj, dims=1:20)
    #plot1 <- JackStrawPlot(obj, dims=1:20) + ylab("100 resampling")
    #plot2 = ElbowPlot(obj)
    #p <- CombinePlots(plots=list(plot1, plot2))
    #pdf(output, width=18, height=6)
    #print(p)
    #dev.off()
    
    obj <- FindNeighbors(object=obj) ## use dims later...
    for (i in seq(0, 1, 0.05)) {
        obj <- FindClusters(object=obj, resolution=i) ## should add algorithm=2 later for better louvain
    }

    #pc.pval  <- obj@reductions$pca@jackstraw@overall.p.values
    #pc.num  <- max(pc.pval[pc.pval[,2] <= 1e-3, 1])
    #obj <- RunTSNE(object=obj, dims=1:pc.num)
    #obj <- RunUMAP(object=obj, dims=1:pc.num)
    ## obj <- RunTSNE(object=obj, dims=1:20)
    obj <- RunUMAP(object=obj, dims=1:20)
    obj
}

scDE.output <- function(obj, x, logfcThresh=0.2, expression_threshold=0.5, fg_expression_threshold=4, fdr=0.05) { 
    ##should be used with logfc differential expression results x
    ##append more information for Seurat v3 output
    results <- foreach(i=levels(x$cluster)) %do% {
        x.i <- x %>% filter(cluster==i)
        x.i <- merge(x.i, human_ortholog, by.x='gene', by.y='Gene') ## would lose some clusters
	#x.i <- cbind(x.i, human_ortholog[match(x.i$gene, human_ortholog$Gene), ])
	print(head(x.i))
	fg.cells <- WhichCells(obj, ident=i)
	bg.cells <- WhichCells(obj, ident=unique(x$cluster)[-(as.numeric(i)+1)])
        #fg.normexp = rowSums(obj$SCT@data[x.i$gene, fg.cells])
        #bg.normexp = rowSums(obj$SCT@data[x.i$gene, bg.cells])
        fg.normexp = rowSums(obj$SCT@counts[x.i$gene, fg.cells])
        bg.normexp = rowSums(obj$SCT@counts[x.i$gene, bg.cells])
	diff.enrich = x.i$pct.1 - x.i$pct.2
	fg.expressed.cells <- apply(obj$SCT@data[x.i$gene, fg.cells], 1, function(x) {length(x[x>expression_threshold])})
	fg.nonexpressed.cells <- length(fg.cells) - fg.expressed.cells
	bg.expressed.cells <- apply(obj$SCT@data[x.i$gene, bg.cells], 1, function(x) {length(x[x>expression_threshold])})
	bg.nonexpressed.cells <- length(bg.cells) - bg.expressed.cells
        fisher.input <- cbind(fg.expressed.cells, fg.nonexpressed.cells, bg.expressed.cells, bg.nonexpressed.cells)
        fisher.p <- apply(fisher.input, 1, function(x) {
          fisher.test(matrix(x, nrow=2, byrow=T))$p.value
        })
        fisher.fdr <- p.adjust(fisher.p, method='fdr')
	y <- cbind(fg.normexp, bg.normexp, x.i, diff.enrich, fisher.p, fisher.fdr)
        df <- data.frame(gene=as.vector(y$gene), cluster_cells_fg=y$fg.normexp, 
                         control_all=y$bg.normexp,
                         logfc_all_cluster_fg_vs_control=y$avg_logFC,
                         fg_fraction=y$pct.1,
                         bg_fraction=y$pct.2,
                         diff_enrichment=y$diff.enrich,
                         fisher_exact_p=y$fisher.p,
                         fisher_fdr=y$fisher.fdr,
                         MAST_diff_p=y$p_val,
                         MAST_diff_fdr=y$p_val_adj,
                         EnsembleID=y$EnsembleID,
                         HSortholog=y$Hsortholog, 
                         HSclosest=y$Hs_closest)
      subset(df, cluster_cells_fg>=fg_expression_threshold & logfc_all_cluster_fg_vs_control>logfcThresh & MAST_diff_fdr<fdr)
    }
    results
}


recluster.withtree <- function(obj, name) {
  obj <- RunPCA(object=obj)
  obj <- JackStraw(obj, num.replicate = 100)
  obj <- ScoreJackStraw(obj, dims=1:20)
  pc.pval  <- obj@reductions$pca@jackstraw@overall.p.values
  pc.num  <- max(pc.pval[pc.pval[,2] <= 1e-3, 1])
  print(pc.num)
  obj <- FindNeighbors(object=obj) ## use the default 10 to keep consistent with previous seurat results, use dims=1:pc.num later 
  for (i in seq(0, 2.0, 0.05)) {
      obj <- FindClusters(object=obj, resolution=i)
  }
  #obj <- RunTSNE(object=obj, dims=1:pc.num)
  obj <- RunUMAP(object=obj, dims=1:pc.num)
  obj
}

