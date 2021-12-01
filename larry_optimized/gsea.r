library(Seurat)
library(tidyverse)
library(vroom)
library(clusterProfiler)
library(DOSE)
library(glue)
library(ComplexHeatmap)
library(GetoptLong)
library(ggplot2)
library(circlize)

allmarkers = read.table("integrative_clusters_markers.txt")

allgenes = readRDS("allgenes_signatures.rds")

clusters = list()
for (i in unique(allmarkers$cluster)) {
   genes = as.character(allmarkers[allmarkers$cluster==i, ]$genes)
   genes = genes[(!is.na(genes)) & (genes!="")]
   cat(i, "\t", length(genes), "\n")
   ### cat(genes, sep='\n', file=glue('lisa/{args$label}_cluster{i}.symbols'))
   y1 = enricher(genes, TERM2GENE=allgenes[,-2], maxGSSize=1000, minGSSize = 5)
   y1 = y1@result
   rownames(y1) = paste0(i, ".", rownames(y1))
   y1$clusters = i
   clusters[[as.character(i)]] <- y1
}

clusters = do.call('rbind', clusters)

clusters = clusters %>% mutate(FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))

readr::write_csv(clusters, path=glue('integrated_res1_all_annotation.csv'))

clusters = subset(clusters, p.adjust <= 0.01 & qvalue <= 0.01 & Count >= 3 & FoldEnrichment >= 1.2)
readr::write_csv(clusters, path=glue('integrated_res1_filter_annotation.csv'))
