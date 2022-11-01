library(glue)
## library("ggVennDiagram")
library(rlang)
library("ggvenn")
library(patchwork)
library(Seurat)
library(clusterProfiler)
library(fgsea)

allemt = list()
for (wk in c("wk6-7", "wk9", "wk12-14", "wk17-18")) {
    skm = Sys.glob(glue('~/langenau/projects/02_atlas/normal_muscle_cellbrowser/skeletal-muscle/*/*-{wk}-myogenic/SkM*.tsv.gz'))
    if (length(skm) > 0) {
        skm = read.table(skm, sep='\t', header=T, stringsAsFactors=F)
        allemt[[paste0(as_name(wk), 'SkM')]]=skm$symbol
    }
    skm = Sys.glob(glue('~/langenau/projects/02_atlas/normal_muscle_cellbrowser/skeletal-muscle/*/*-{wk}-myogenic/*MC*.tsv.gz'))[1]
    if (length(skm) > 0) {
        skm = read.table(skm, sep='\t', header=T, stringsAsFactors=F)
        allemt[[paste0(as_name(wk), "MBMC")]]=skm$symbol
    }
    skm = Sys.glob(glue('~/langenau/projects/02_atlas/normal_muscle_cellbrowser/skeletal-muscle/*/*-{wk}-myogenic/*MP*.tsv.gz'))[1]
    if (length(skm) > 0) {
        skm = read.table(skm, sep='\t', header=T, stringsAsFactors=F)
        allemt[[paste0(as_name(wk), "MP")]]=skm$symbol
    }
}


pdf("MSK74711_GSEAplot_top10.pdf", width=8, height=12.5)
msk74.test.list = foreach(state = c("EMT", "Muscle", "Prolif")) %do% {
    markers = read.csv(glue('MSK74711_{state}_DE_genes.csv'))
    markers = sort(setNames(markers$avg_log2FC, markers$X), decreasing = T)
    ## markers = sort(setNames(rank(markers$avg_log2FC), markers$X), decreasing = T)
    markers = markers[!duplicated(markers)]
    fgseaRes = fgsea(allemt, markers, nperm=1e5)
    fgseaRes.pos = fgseaRes[NES>0][order(pval), pathway]
    fgseaRes.neg = fgseaRes[NES<0][order(pval), pathway]
    fgseaRes.pathway = c(fgseaRes.pos, rev(fgseaRes.neg))
    wrap_elements(plotGseaTable(allemt[fgseaRes.pathway], markers, fgseaRes, colwidths=c(1, 3, 0.6, 0.6, 0.6), gseaParam = 0.5, render=F)) + ggtitle(state)
    ## fgseaRes
}
msk74.test.list[[1]]/msk74.test.list[[2]]/msk74.test.list[[3]]
dev.off()

pdf("MAST139_GSEAplot_top10.pdf", width=8, height=12.5)
mast.test.list = foreach(state = c("EMT", "Muscle", "Prolif")) %do% {
    markers = read.csv(glue('MAST139_{state}_DE_genes.csv'))
    markers = sort(setNames(markers$avg_log2FC, markers$X), decreasing = T)
    ## markers = sort(setNames(rank(markers$avg_log2FC), markers$X), decreasing = T)
    markers = markers[!duplicated(markers)]
    fgseaRes = fgsea(allemt, markers, nperm=1e5)
    fgseaRes.pos = fgseaRes[NES>0][order(pval), pathway]
    fgseaRes.neg = fgseaRes[NES<0][order(pval), pathway]
    fgseaRes.pathway = c(fgseaRes.pos, rev(fgseaRes.neg))
    wrap_elements(plotGseaTable(allemt[fgseaRes.pathway], markers, fgseaRes, colwidths=c(1, 3, 0.6, 0.6, 0.6), gseaParam = 0.5, render=F)) + ggtitle(state)
    ## fgseaRes
}
mast.test.list[[1]]/mast.test.list[[2]]/mast.test.list[[3]]
dev.off()


## topPathwaysUp <- fgseaRes[NES > 0][head(order(pval), n=10), pathway]
