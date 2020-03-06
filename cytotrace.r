library(CytoTRACE)
library(uwot)

## all PDX samples
cols <- read.table('color_table.xls', sep='\t', header=T, check.names=F)

## raw.erms.results <- conos:::papply(colnames(cols)[2:9], function(sam) {
raw.erms.results <- conos:::papply(c("MAST35", "MAST85", "MAST139"), function(sam) {
    ## rds = Sys.glob(paste0('../results/seurat_intersect_velocity/*', sam, '_seu.rds'))
    rds = Sys.glob(paste0('/data/langenau/human_rms_pdxs/seurat_objects/*', sam, '.rds'))
    primary1 <- readRDS(rds)
    outputDir <- gsub('.rds', '', basename(rds))
    outputDir2 <- gsub('.rds', '_2', basename(rds))
    print(outputDir)
    input.mat <- as.matrix(primary1@assays$RNA@counts)
    colnames(input.mat) = paste0(colnames(input.mat), '_', sam)
    cyto <- CytoTRACE(input.mat)
    levels(primary1@meta.data$RNA_snn_res.0.8) = as.vector(cols[, sam][cols[, sam]!=""])
    pheno <- as.vector(primary1@meta.data$RNA_snn_res.0.8)
    names(pheno) <- colnames(input.mat)
    sum(names(cyto$CytoTRACE) == names(pheno))
    emb = primary1@reductions$umap@cell.embeddings
    rownames(emb) = names(pheno)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOD1", outputDir=outputDir,
                  emb=emb)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "VIM", outputDir=outputDir2,
                  emb=emb)
    df <- data.frame(CytoTRACE=cyto$CytoTRACE,
                     pheno=pheno)
    library(ggplot2)
    ggplot(df, aes(y=CytoTRACE, x=pheno)) + geom_boxplot()
    ggsave(paste0(sam, '_box.pdf'))
    list(input.mat, pheno)
})

erms.results <- iCytoTRACE(lapply(raw.erms.results, function(x) x[[1]]))

erms.pheno <- unlist(lapply(raw.erms.results, function(x) {x[2]}))

erms.results.bak <- erms.results

umap <- umap(erms.results.bak$coord, n_neighbors=50)

rownames(umap) <- colnames(erms.results.bak$exprMatrix)[!colnames(erms.results.bak$exprMatrix)%in%erms.results.bak$filteredCells]

erms.results.bak$CytoTRACE <- erms.results.bak$CytoTRACE[rownames(umap)]
erms.results.bak$exprMatrix <- erms.results.bak$exprMatrix[, rownames(umap)]

plotCytoTRACE(erms.results.bak, phenotype=erms.pheno[rownames(umap)], gene = "VIM", outputDir='erms_results',
              emb=umap)

erms.df <- data.frame(CytoTRACE=erms.results$CytoTRACE,
                      pheno=as.factor(erms.pheno))

library(ggplot2)
library(tidyverse)

erms.df %>% drop_na() %>%
    mutate(pheno = reorder(pheno, CytoTRACE, .fun=median, .desc =F)) %>%
    ggplot(aes(x=fct_reorder(pheno, CytoTRACE, .fun = median, .desc =TRUE), y=CytoTRACE)) + geom_boxplot(aes(fill=pheno))+    geom_jitter(aes(alpha=0.1), position=position_jitter(0.25)) + xlab("Phenotype") + ggpubr:::theme_pubr()
ggsave('erms_box.pdf', width=12)

raw.erms.results <- conos:::papply(colnames(cols)[2:9], function(sam) {
    ## rds = Sys.glob(paste0('../results/seurat_intersect_velocity/*', sam, '_seu.rds'))
    rds = Sys.glob(paste0('/data/langenau/human_rms_pdxs/seurat_objects/*', sam, '.rds'))
    primary1 <- readRDS(rds)
    outputDir <- gsub('.rds', '', basename(rds))
    outputDir2 <- gsub('.rds', '_2', basename(rds))
    print(outputDir)
    input.mat <- as.matrix(primary1@assays$RNA@counts)
    colnames(input.mat) = paste0(colnames(input.mat), '_', sam)
    cyto <- CytoTRACE(input.mat)
    levels(primary1@meta.data$RNA_snn_res.0.8) = as.vector(cols[, sam][cols[, sam]!=""])
    pheno <- as.vector(primary1@meta.data$RNA_snn_res.0.8)
    names(pheno) <- colnames(input.mat)
    sum(names(cyto$CytoTRACE) == names(pheno))
    emb = primary1@reductions$umap@cell.embeddings
    rownames(emb) = names(pheno)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOD1", outputDir=outputDir,
                  emb=emb)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "VIM", outputDir=outputDir2,
                  emb=emb)
    df <- data.frame(CytoTRACE=cyto$CytoTRACE,
                     pheno=pheno)
    library(ggplot2)
    ggplot(df, aes(y=CytoTRACE, x=pheno)) + geom_boxplot()
    ggsave(paste0(sam, '_box.pdf'))
    list(input.mat, pheno)
})

raw.arms.results <- conos:::papply(colnames(cols)[10:ncol(cols)], function(sam) {
    rds = Sys.glob(paste0('/data/langenau/human_rms_pdxs/seurat_objects/*', sam, '.rds'))
    primary1 <- readRDS(rds)
    outputDir <- gsub('.rds', '', basename(rds))
    outputDir2 <- gsub('.rds', '_2', basename(rds))
    print(outputDir)
    input.mat <- as.matrix(primary1@assays$RNA@counts)
    colnames(input.mat) = paste0(colnames(input.mat), '_', sam)
    cyto <- CytoTRACE(input.mat)
    levels(primary1@meta.data$RNA_snn_res.0.8) = as.vector(cols[, sam][cols[, sam]!=""])
    pheno <- as.vector(primary1@meta.data$RNA_snn_res.0.8)
    names(pheno) <- colnames(input.mat)
    sum(names(cyto$CytoTRACE) == names(pheno))
    emb = primary1@reductions$umap@cell.embeddings
    rownames(emb) = names(pheno)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOD1", outputDir=outputDir,
                  emb=emb)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "VIM", outputDir=outputDir2,
                  emb=emb)
    df <- data.frame(CytoTRACE=cyto$CytoTRACE,
                     pheno=pheno)
    library(ggplot2)
    ggplot(df, aes(y=CytoTRACE, x=pheno)) + geom_boxplot()
    ggsave(paste0(sam, '_box.pdf'))
    list(input.mat, pheno)
})

arms.results <- iCytoTRACE(lapply(raw.arms.results, function(x) x[[1]]))

arms.pheno <- lapply(1:length(raw.arms.results),function(i)
    raw.arms.results[[i]][[2]][na.omit(match(colnames(arms.results$exprMatrix), colnames(raw.arms.results[[i]][[1]])))])
arms.pheno <- do.call(c, arms.pheno)

umap = umap(arms.results$coord, n_neighbors=50)

rownames(umap) <- colnames(arms.results$exprMatrix)

plotCytoTRACE(arms.results, phenotype = arms.pheno, gene = "VIM", outputDir='arms_results',
              emb=umap)

arms.df <- data.frame(CytoTRACE=arms.results$CytoTRACE,
                      pheno=arms.pheno)

library(ggplot2)
arms.df %>% drop_na() %>%
    mutate(pheno = reorder(pheno, CytoTRACE, .fun=median, .desc =F)) %>%
    ggplot(aes(x=fct_reorder(pheno, CytoTRACE, .fun = median, .desc =TRUE), y=CytoTRACE)) + geom_boxplot(aes(fill=pheno))+    geom_jitter(aes(alpha=0.1), position=position_jitter(0.25)) + xlab("Phenotype") + ggpubr:::theme_pubr()
ggsave('arms_box.pdf', width=12)


fish.cols <- read.delim('fish_color_table.txt', sep='\t', header=T, check.names=F)

raw.erms.results <- conos:::papply(colnames(fish.cols)[2:ncol(fish.cols)], function(sam) {
    rds = Sys.glob(paste0('../results/seurat_intersect_velocity/*', sam, '_seu.rds'))
    print(rds)
    ## rds = Sys.glob(paste0('/data/langenau/human_rms_pdxs/seurat_objects/*', sam, '.rds'))
    primary1 <- readRDS(rds)
    outputDir <- gsub('.rds', '', basename(rds))
    outputDir2 <- gsub('.rds', '_2', basename(rds))
    print(outputDir)
    input.mat <- as.matrix(primary1@assays$RNA@counts)
    colnames(input.mat) = paste0(colnames(input.mat), '_', sam)
    cyto <- CytoTRACE(input.mat)
    levels(primary1@meta.data$seurat_clusters) = as.vector(fish.cols[, sam][fish.cols[, sam]!=""])
    pheno <- as.vector(primary1@meta.data$seurat_clusters)
    names(pheno) <- colnames(input.mat)
    sum(names(cyto$CytoTRACE) == names(pheno))
    emb = primary1@reductions$umap@cell.embeddings
    rownames(emb) = names(pheno)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOD1", outputDir=outputDir,
                  emb=emb)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "VIM", outputDir=outputDir2,
                  emb=emb)
    df <- data.frame(CytoTRACE=cyto$CytoTRACE,
                     pheno=pheno)
    library(ggplot2)
    ggplot(df, aes(y=CytoTRACE, x=pheno)) + geom_boxplot()
    ggsave(paste0(sam, '_box.pdf'))
    list(input.mat, pheno)
})
