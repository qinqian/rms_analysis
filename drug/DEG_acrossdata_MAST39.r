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
                rgb(0, 255, 253, maxColorValue = 255))
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'UPR', "Neural")
names(metacolors) <- metalabels

library(Seurat)
library(reticulate)
library(patchwork)
library(ggplot2)


res <- readRDS('results/seurat_sara/MAST-39-Resistant_seurat-object.rds')
sen <- readRDS('results/seurat_sara/MAST-39-Sensitive_seurat-object.rds')

res$seurat_clusters = res$RNA_snn_res.0.8
sen$seurat_clusters = sen$RNA_snn_res.0.8
cell4 = list(res, sen)

integration <- function(x, y, label, method='seurat') {
    Idents(x) = x$seurat_clusters2 = x$seurat_clusters
    states = unlist(as.vector(annotation[label[1], ]))
    print(states)
    states = states[states!='']
    states = states[!is.na(states)]
    print('----------')
    print(states)
    print(levels(x$seurat_clusters))
    levels(x$seurat_clusters) = states
    Idents(x) = x$seurat_clusters
    print('-------')
    ## levels(Idents(x)) = states
    Idents(y) = y$seurat_clusters2 = y$seurat_clusters
    states = unlist(as.vector(annotation[label[2], ]))
    states = na.omit(states[states!=''])
    print(length(states))
    print(levels(y$seurat_clusters))
    levels(y$seurat_clusters) = states
    Idents(y) = y$seurat_clusters
    ## levels(Idents(y)) = states
    seurat.pseudo.list = list(x, y)
    for (i in 1:length(seurat.pseudo.list)) {
        print(i)
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
annotation = read.delim('clusters.txt', sep='\t', row.names=1, header=T,
                        check.names=F, stringsAsFactors=F)

integrated.seurat4 = integration(cell4[[1]], cell4[[2]], c('MAST39-Resistant', 'MAST39-Sensitive'))

DefaultAssay(integrated.seurat4) <- "RNA"

integrated.seurat4$cellstate.stim = paste(Idents(integrated.seurat4), integrated.seurat4$orig.ident, sep = "_")
integrated.seurat4$cellstate = Idents(integrated.seurat4)
Idents(integrated.seurat4) <- "cellstate.stim"

plots <- list()

for (cell in c("Ground", "Hypoxia", "Prolif")) {
    EMT1 = subset(integrated.seurat4, ident=paste0(cell, '_MAST-39-Resistant'))
    EMT2 = subset(integrated.seurat4, ident=paste0(cell, '_MAST-39-Sensitive'))
    EMT = merge(EMT1, EMT2)
    avg.emt <- log1p(AverageExpression(EMT, verbose = FALSE)$RNA)
    targets = c("PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "AKT1", "AKT2", "AKT3", "MTOR", "ABCB1", "ABCC1", "ABCG2", "FUT1", "DYDC2", "ROR2")
    targets = paste0('hg19-', targets)
    avg.emt <- data.frame(avg.emt)
    print(head(avg.emt))
    colnames(avg.emt) = c("Resistant", "Sensitive")
    print(head(avg.emt))
    print(avg.emt[targets, ])
    p1 <- ggplot(avg.emt,
                 aes(Sensitive, Resistant)) + geom_point() + ggtitle(cell)
    p1 <- LabelPoints(plot = p1, points = targets, repel = TRUE)
    plots[[cell]] = p1
}

wrap_plots(plots = plots, ncol = 4)
ggsave('MAST39_qiqi_target_genes_scatter.pdf', height=5, width=16)

plots <- VlnPlot(integrated.seurat4, features = targets, split.by = "orig.ident", group.by = "cellstate", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
ggsave('MAST39_qiqi_target_genes.pdf', height=40, width=8)


pdf('ALL_MAST39_batchcorrected_targets_umap.pdf', width=12, height=5.5)
p4d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
print(p4d)
## for (tar in targets) {
##     head(integrated.seurat4@meta.data)
##     if (sum(integrated.seurat4[tar, integrated.seurat4$orig.ident=='RH41-R']$RNA@data)!=0) {
##         p5a=FeaturePlot(integrated.seurat4[, integrated.seurat4$orig.ident=='RH41-R'],
##                         tar, reduction='umap', label=F, ncol=1) + theme(legend.position='right')
##     } else {
##         p5a=FeaturePlot(integrated.seurat4[, integrated.seurat4$orig.ident=='RH41-R'],
##                         tar, reduction='umap', label=F, ncol=1, col=c("lightgrey", "lightgrey")) + theme(legend.position='right')
##     }
##     if (sum(integrated.seurat4[tar, integrated.seurat4$orig.ident=='RH41-S']$RNA@data)!=0) {
##         p5b=FeaturePlot(integrated.seurat4[, integrated.seurat4$orig.ident=='RH41-S'],
##                         tar, reduction='umap', label=F, ncol=1) + theme(legend.position='right')
##     } else {
##         p5b=FeaturePlot(integrated.seurat4[, integrated.seurat4$orig.ident=='RH41-S'],
##                         tar, reduction='umap', label=F, ncol=1, col=c("lightgrey", "lightgrey")) + theme(legend.position='right')
##     }
##     print(p5a|p5b)
## }
dev.off()


pdf('ALL_MAST39_batchcorrected_targets_umap_drugmarkers.pdf', width=10, height=10)
## p5d=FeaturePlot(integrated.seurat4, targets[targets!='PIK3CG'], reduction='umap', split.by='orig.ident')
p5d=FeaturePlot(integrated.seurat4, c("hg19-FUT1", "hg19-DYDC2", "hg19-ROR2"), reduction='umap', split.by='orig.ident')
print(p5d)
dev.off()


DEGs = list()
for (cell in c("Ground", "Hypoxia", "Prolif")) {
    response <- FindMarkers(integrated.seurat4, ident.1 = paste0(cell, "_MAST-39-Resistant"), ident.2 = paste0(cell, "_MAST-39-Sensitive"), verbose = FALSE)
    response$cellstate = cell
    response$enrich = response$pct.1 - response$pct.2
    response$genes = rownames(response)
    DEGs[[cell]] = response
}

DEGs.df = dplyr::bind_rows(DEGs)

get_msigdb = function(category='c2') {
    c2 = read.gmt(Sys.glob(glue("/PHShome/qq06/langenau/01_rms_projects/02_human/data/msigdb/{category}.all.v7.0.symbols.gmt")))
    c2
}

library(clusterProfiler)
library(DOSE)
library(ComplexHeatmap)
library(GetoptLong)
library(ggplot2)
library(circlize)
library(glue)

c2 = get_msigdb()
h = get_msigdb('h')
c5 = get_msigdb('c5')

allgenes = dplyr::bind_rows(c2, h, c5)

GSEAsig = list()
for (cell in c("Ground", "Hypoxia", "Prolif")) { 
    y = DEGs.df%>%filter(p_val_adj<0.01)%>%filter(cellstate==cell)
    y.up = y %>% filter(abs(enrich)>=0.1)%>%filter(avg_log2FC>0)
    y.dn = y %>% filter(abs(enrich)>=0.1)%>%filter(avg_log2FC<0)
    print(head(y.dn$genes))
    y.dn$genes = gsub('hg19-', '', y.dn$genes)
    y.up$genes = gsub('hg19-', '', y.up$genes)
    y.up = enricher((y.up%>%select(genes))$genes,
                    TERM2GENE=allgenes, maxGSSize=1000, minGSSize = 5)@result
    y.dn = enricher((y.dn%>%select(genes))$genes,
                    TERM2GENE=allgenes, maxGSSize=1000, minGSSize = 5)@result
    y.up$cell = paste0(cell, '_UP')
    y.dn$cell = paste0(cell, '_DN')
    GSEAsig[[cell]] = rbind(y.up, y.dn)
}

GSEAsig.df = dplyr::bind_rows(GSEAsig)
## GSEAsig.df = GSEAsig.df[GSEAsig.df$p.adjust <= 0.01, ]

GSEAsig.df = GSEAsig.df[order(GSEAsig.df$cell, GSEAsig.df$p.adjust), ]

readr::write_csv(DEGs.df, "MAST39_qiqi_DEG_across_cond.csv")
readr::write_csv(GSEAsig.df, "MAST39_qiqi_GSEAsig_across_cond.csv")
