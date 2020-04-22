library(CytoTRACE)
library(foreach)
library(doMC)
library(uwot)
library(glue)
library(patchwork)
registerDoMC(3)

metacolors <- c(rgb(166, 166, 166, maxColorValue = 255),
                rgb(241, 149, 69, maxColorValue = 255),
                rgb(103, 35,  102, maxColorValue = 255),
                rgb(52, 101, 252, maxColorValue = 255),
                rgb(242, 242, 242, maxColorValue = 255),
                rgb(52, 101, 252, maxColorValue = 255),
                rgb(233, 63,  51, maxColorValue = 255),
                rgb(65, 129,  7, maxColorValue = 255),
                rgb(52, 101, 252, maxColorValue = 255),
                rgb(253, 247, 49, maxColorValue = 255),
                'purple', 'gray', 'gray', 'gray', 'gray', 'gray')
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "INTERFERON", "Prolif",
                "Histone", "TNFA", 'Unique #5', 'Unique #4', 'Unique #3', 'Unique #2', 'Unique #6')
names(metacolors) <- metalabels

cols = read.delim('../final_annotations/fish_clusters.txt', sep='\t', row.names=1, header=T,
                  check.names=F, stringsAsFactors=F)
primary <- lapply(c('../results/seurat/Tumor21_unfilter_seurat_obj_tumors.rds',
                    '../results/seurat/Tumor22_unfilter_seurat_obj_tumors.rds',
                    '../results/seurat/Tumor24_unfilter_seurat_obj_tumors.rds'), readRDS)
primary.tumors <- lapply(c('../results/seurat_v6/Tumor21_recluster1.8.rds',
                           '../results/seurat_v6/Tumor22_recluster1.8.rds',
                           '../results/seurat_intersect_velocity/Tumor24_seu.rds'), readRDS)
labels <- c("Tumor21", "Tumor22", "Tumor24")

raw.erms.results <- foreach(i=1:length(labels)) %dopar% {
    primary1 = primary.tumors[[i]]
    input.mat <- as.matrix(primary1@assays$RNA@counts)
    colnames(input.mat) = paste0(colnames(input.mat), '_', labels[i])
    cyto <- CytoTRACE(input.mat)
    states = unlist(as.vector(cols[labels[i], ]))
    states = states[states!='']
    levels(primary1$seurat_clusters) = states
    pheno <- as.vector(primary1$seurat_clusters)
    names(pheno) <- colnames(input.mat)
    emb <- primary1@reductions$umap@cell.embeddings
    rownames(emb) <- names(pheno)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOD1", outputDir=glue('{labels[i]}_cytotrace'),
                  emb=emb)
    ## plotCytoTRACE(cyto, phenotype = pheno, gene = "MYF5", outputDir=glue('{labels[i]}_cytotrace2'),
    ##               emb=emb, colors=metacolors)
    ## plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOG", outputDir=glue('{labels[i]}_cytotrace3'),
    ##               emb=emb, colors=metacolors)
    ## plotCytoTRACE(cyto, phenotype = pheno, gene = "MYLPF", outputDir=glue('{labels[i]}_cytotrace4'),
    ##               emb=emb, colors=metacolors)
    df <- data.frame(CytoTRACE=cyto$CytoTRACE,
                     pheno=pheno)
    library(ggplot2)
    ## ggplot(df, aes(y=CytoTRACE, x=pheno)) + geom_boxplot()
    df %>% mutate(pheno = reorder(pheno, CytoTRACE, .fun=median, .desc =F)) %>%
    ggplot(aes(x=fct_reorder(pheno, CytoTRACE, .fun = median, .desc =TRUE), y=CytoTRACE)) + geom_boxplot(aes(fill=pheno)) + scale_fill_manual(values=metacolors) + xlab("Phenotype") + ggpubr:::theme_pubr() #+ geom_jitter(aes(alpha=0.1), position=pos
    ggsave(glue('{labels[i]}_cytotrace_box.pdf'))
    list(input.mat, pheno, df, emb, labels[i])
}

lapply(raw.erms.results, function(x) {
    emb = x[[4]]
    res = cbind(emb, x[[3]])
    print(head(res))
    p1=ggplot(res, aes(UMAP_1, UMAP_2, colour=pheno)) + geom_point() +
        scale_color_manual(values=metacolors) + ggpubr::theme_pubr()
    p2=ggplot(res, aes(UMAP_1, UMAP_2, colour=CytoTRACE)) + geom_point() +
        scale_colour_gradient2(low="white", high="red") + ggpubr::theme_pubr()
    p1 + p2
    ggsave(glue('{x[[5]]}_umap_cytotrace.pdf'), width=12.5, height=5)
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
    ggplot(aes(x=fct_reorder(pheno, CytoTRACE, .fun = median, .desc =TRUE), y=CytoTRACE)) + geom_boxplot(aes(fill=pheno)) + scale_fill_manual(values=metacolors) + xlab("Phenotype") + ggpubr:::theme_pubr() #+ geom_jitter(aes(alpha=0.1), position=position_jitter(0.25)) 
ggsave('Fig5_fish_erms_box.pdf', width=6.5, height=5)


pdxs = Sys.glob('../data/seurat_obj/*rds')[1:10]
cols = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1,
                  header=T,
                  check.names=F, stringsAsFactors=F)

pdxs = c(pdxs, '../results/seurat_sara/20191031_MSK74711_seurat-object.rds',
         '../results/seurat_sara/MAST118_seurat-object.rds',
         '../results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds',
         '../results/seurat_sara/MAST139_1cells_seurat-object.rds')

pdxs.objs = lapply(pdxs, readRDS)

labels = unlist(lapply(pdxs.objs, function(x) {
    levels(x$orig.ident[1])
}))
labels[9] = 'MAST85-1'
labels[11] = 'MSK74711'
labels[13] = 'MSK72117'
labels[14] = 'MAST139-1'
labels[10] = 'RH74-10'

labels.show = labels[-grep("-", labels)]
pdxs.objs = pdxs.objs[-grep("-", labels)]

subtype = c("ERMS",
            "ERMS",
            "ERMS",
            "ERMS",
            "ERMS",
            "ARMS",
            "ARMS",
            "ERMS",
            "ERMS",
            "ARMS", "ARMS")

## raw.erms.results <- foreach(i=which(subtype=='ERMS')) %dopar% {
## raw.erms.results <- foreach(i=which(subtype=='ERMS')) %dopar% {
raw.erms.results <- foreach(i=which(labels.show %in% c("MAST111", "MAST139", "MSK74711", 'MAST39'))) %dopar% {
    primary1 = pdxs.objs[[i]]
    input.mat <- as.matrix(primary1@assays$RNA@counts)
    colnames(input.mat) = paste0(colnames(input.mat), '_', labels.show[i])
    cyto <- CytoTRACE(input.mat)
    states = unlist(as.vector(cols[labels.show[i], ]))
    states = states[states!='']
    primary1$seurat_clusters = primary1$RNA_snn_res.0.8
    levels(primary1$seurat_clusters) = states
    pheno <- as.vector(primary1$seurat_clusters)
    names(pheno) <- colnames(input.mat)
    emb <- primary1@reductions$umap@cell.embeddings
    rownames(emb) <- names(pheno)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOD1", outputDir=glue('{labels.show[i]}_cytotrace'),
                  emb=emb)
    ## plotCytoTRACE(cyto, phenotype = pheno, gene = "MYF5", outputDir=glue('{labels[i]}_cytotrace2'),
    ##               emb=emb, colors=metacolors)
    ## plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOG", outputDir=glue('{labels[i]}_cytotrace3'),
    ##               emb=emb, colors=metacolors)
    ## plotCytoTRACE(cyto, phenotype = pheno, gene = "MYLPF", outputDir=glue('{labels[i]}_cytotrace4'),
    ##               emb=emb, colors=metacolors)
    df <- data.frame(CytoTRACE=cyto$CytoTRACE,
                     pheno=pheno)
    library(ggplot2)
    ggplot(df, aes(y=CytoTRACE, x=pheno)) + geom_boxplot()
    ggsave(glue('{labels.show[i]}_cytotrace_box.pdf'))
    list(input.mat, pheno, df, emb, labels.show[i])
}

lapply(raw.erms.results, function(x) {
    emb = x[[4]]
    res = cbind(emb, x[[3]])
    print(head(res))
    p1=ggplot(res, aes(UMAP_1, UMAP_2, colour=pheno)) + geom_point() +
        scale_color_manual(values=metacolors) + ggpubr::theme_pubr()
    p2=ggplot(res, aes(UMAP_1, UMAP_2, colour=CytoTRACE)) + geom_point() +
        scale_colour_gradient2(low="white", high="red") + ggpubr::theme_pubr()
    p1 + p2
    ggsave(glue('{x[[5]]}_umap_cytotrace.pdf'), width=12.5, height=5)
})

erms.results <- iCytoTRACE(lapply(raw.erms.results, function(x) x[[1]]))

erms.pheno <- unlist(lapply(raw.erms.results, function(x) {x[2]}))
erms.results.bak <- erms.results
umap <- umap(erms.results.bak$coord, n_neighbors=50)
rownames(umap) <- colnames(erms.results.bak$exprMatrix)[!colnames(erms.results.bak$exprMatrix)%in%erms.results.bak$filteredCells]
erms.results.bak$CytoTRACE <- erms.results.bak$CytoTRACE[rownames(umap)]
erms.results.bak$exprMatrix <- erms.results.bak$exprMatrix[, rownames(umap)]
plotCytoTRACE(erms.results.bak, phenotype=erms.pheno[rownames(umap)], gene = "VIM", outputDir='pdxs_erms_results',
              emb=umap)
erms.df <- data.frame(CytoTRACE=erms.results$CytoTRACE,
                      pheno=as.factor(erms.pheno))
library(ggplot2)
library(tidyverse)
erms.df %>% drop_na() %>%
    mutate(pheno = reorder(pheno, CytoTRACE, .fun=median, .desc =F)) %>%
    ggplot(aes(x=fct_reorder(pheno, CytoTRACE, .fun = median, .desc =TRUE), y=CytoTRACE)) + geom_boxplot(aes(fill=pheno)) + scale_fill_manual(values=metacolors) + xlab("Phenotype") + ggpubr:::theme_pubr() #+ geom_jitter(aes(alpha=0.1), position=position_jitter(0.25)) 
ggsave('Fig5_pdxs_erms_box.pdf', width=6.5, height=5)

raw.arms.results <- foreach(i=which(subtype=='ARMS')) %dopar% {
    primary1 = pdxs.objs[[i]]
    input.mat <- as.matrix(primary1@assays$RNA@counts)
    colnames(input.mat) = paste0(colnames(input.mat), '_', labels.show[i])
    cyto <- CytoTRACE(input.mat)
    states = unlist(as.vector(cols[labels.show[i], ]))
    states = states[states!='']
    primary1$seurat_clusters = primary1$RNA_snn_res.0.8
    levels(primary1$seurat_clusters) = states
    pheno <- as.vector(primary1$seurat_clusters)
    names(pheno) <- colnames(input.mat)
    emb <- primary1@reductions$umap@cell.embeddings
    rownames(emb) <- names(pheno)
    plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOD1", outputDir=glue('{labels.show[i]}_cytotrace'),
                  emb=emb)
    ## plotCytoTRACE(cyto, phenotype = pheno, gene = "MYF5", outputDir=glue('{labels[i]}_cytotrace2'),
    ##               emb=emb, colors=metacolors)
    ## plotCytoTRACE(cyto, phenotype = pheno, gene = "MYOG", outputDir=glue('{labels[i]}_cytotrace3'),
    ##               emb=emb, colors=metacolors)
    ## plotCytoTRACE(cyto, phenotype = pheno, gene = "MYLPF", outputDir=glue('{labels[i]}_cytotrace4'),
    ##               emb=emb, colors=metacolors)
    df <- data.frame(CytoTRACE=cyto$CytoTRACE,
                     pheno=pheno)
    library(ggplot2)
    ggplot(df, aes(y=CytoTRACE, x=pheno)) + geom_boxplot()
    ggsave(glue('{labels.show[i]}_cytotrace_box.pdf'))
    list(input.mat, pheno, df, emb, labels.show[i])
}

lapply(raw.arms.results, function(x) {
    emb = x[[4]]
    res = cbind(emb, x[[3]])
    print(head(res))
    p1=ggplot(res, aes(UMAP_1, UMAP_2, colour=pheno)) + geom_point() +
        scale_color_manual(values=metacolors) + ggpubr::theme_pubr()
    p2=ggplot(res, aes(UMAP_1, UMAP_2, colour=CytoTRACE)) + geom_point() +
        scale_colour_gradient2(low="white", high="red") + ggpubr::theme_pubr()
    p1 + p2
    ggsave(glue('{x[[5]]}_ARMS_umap_cytotrace.pdf'), width=12.5, height=5)
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
    ggplot(aes(x=fct_reorder(pheno, CytoTRACE, .fun = median, .desc =TRUE), y=CytoTRACE)) + geom_boxplot(aes(fill=pheno))+scale_fill_manual(values=metacolors) + xlab("Phenotype") +  ggpubr:::theme_pubr()
ggsave('Fig5_arms_box.pdf', width=12)
