library(Seurat)
library(reticulate)
library(tidyverse)
library(patchwork)

cells = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds')

annotation = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1, header=T,
                        check.names=F, stringsAsFactors=F)

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
                rgb(241, 250, 100, maxColorValue = 255))
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'UPR')
names(metacolors) <- metalabels

states = unlist(as.vector(annotation["MAST139", ]))
states = states[states!='']

levels(cells$RNA_snn_res.0.8) = states

Idents(cells) = cells$RNA_snn_res.0.8

de = FindMarkers(cells, only.pos=F, min.pct=0,
                 ident.1='EMT',
                 test.use='MAST', min.diff.pct=0,
                 random.seed=100, logfc.threshold = 0)

write.csv(de, file='MAST139_DE_genes.csv', row.names=T, quote=F)

library(msigdbr)
library(DOSE)
library(org.Hs.eg.db)

m_t2g <- msigdbr(species = "Homo sapiens", category = c("H")) %>% dplyr::select(gs_name, entrez_gene)

geneList = sort(setNames(de$avg_logFC, rownames(de)), decreasing=T)

symbol2entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                       key=names(geneList),
                                       columns="ENTREZID",
                                       keytype="SYMBOL")

names(geneList) = symbol2entrez[match(names(geneList), symbol2entrez[,1]), 2]
geneList = geneList[!is.na(names(geneList))]
geneList = geneList[!duplicated(names(geneList))]
geneList = geneList[!duplicated(geneList)]

geneList.out = geneList
names(geneList.out) = symbol2entrez[match(names(geneList), symbol2entrez[, 2]), 1]

write.csv(geneList.out, file='MAST139_DE_genes_preranked.csv', row.names=T, quote=F)

library(clusterProfiler)
library(fgsea)
## https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
em2 <- GSEA(geneList, TERM2GENE = m_t2g,  minGSSize = 5, seed=9, nPerm=100000, by='fgsea')
result = em2@result
result.sortes = result[order(result$NES, decreasing=T), ]

p1 <- gseaplot(em2, geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', by = "all", title = result.sortes$Description[1])
ggsave("MAST139_GSEAplot.pdf", width=12, height=12)

pathway = list()
for (i in 1:nrow(m_t2g)) {
    pathway[[m_t2g$gs_name[i]]] = c(pathway[[m_t2g$gs_name[i]]], m_t2g$entrez_gene[i])
}
fgseaRes = fgsea(pathways=pathway,
      stats=geneList, minSize=5, maxSize=500, nperm = 1e5)
topPathwaysUp <- fgseaRes[NES > 0][head(order(pval), n=10), pathway]

pdf("MAST139_GSEAplot_top10.pdf", width=12, height=5)
plotGseaTable(pathway[topPathways], geneList, fgseaRes, colwidths=c(6, 3, 0.6, 1, 1),
              gseaParam = 0.5)
dev.off()
