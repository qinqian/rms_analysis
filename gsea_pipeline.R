## https://stephenturner.github.io/deseq-to-fgsea/
## https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
## https://yulab-smu.github.io/clusterProfiler-book/chapter5.html
library(fgsea)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(reactome.db)
library(tidyverse)
library(dplyr)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--seuratde', dest='seurat', metavar='N', type="character", nargs='+')
    parser$add_argument('--label', dest='label', default='')
    parser$add_argument('--species', dest='species', type='character', default='human')
    args = parser$parse_args()
    args
}

args = get_args()
if (length(args$seurat) == 0 || args$label == '') {
    cat('empty argument, exit..')
    q()
}

markers = read.table(args$seurat, sep='\t', stringsAsFactors=F)
markers$cluster = as.character(markers$cluster)
print(head(markers))

pathways.hallmark <- gmtPathways("../data/msigdb/c2.all.v7.0.symbols.gmt")

if (args$species == 'human') {
   symbol2entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                          key=markers$gene,
                                          columns="ENTREZID",
                                          keytype="SYMBOL")
} else {
   symbol2entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                          key=markers$Hsortholog,
                                          columns="ENTREZID",
                                          keytype="SYMBOL")
}

markers.final = merge(markers, symbol2entrez, by.x='gene', by.y='SYMBOL')

system('mkdir -p ../results/gsea')
system('mkdir -p ../results/gsea/lisa')

gsea.results = list()
print(gsea.results)


for (cluster in sort(unique(markers.final$cluster))) {
    print(cluster)
    gsea.results[[cluster]] = list()
    a=as.numeric(cluster)

    if (args$species == 'fish') {
        cluster_stat = markers.final%>%filter(cluster==a) %>%
            distinct() %>%
            group_by(Hsortholog) %>%
            summarize(stat=mean(avg_logFC)) %>%
            arrange(desc(stat)) %>% filter(stat>0.2)
    } else {
        cluster_stat = markers.final%>%filter(cluster==a) %>%
            distinct() %>%
            group_by(gene) %>%
            summarize(stat=mean(avg_logFC)) %>%
            arrange(desc(stat)) %>% filter(stat>0.2)
    }
    cluster_stat2 = markers.final%>%filter(cluster==a) %>%
        distinct() %>%
        group_by(ENTREZID) %>%
        summarize(stat=mean(avg_logFC)) %>%
        arrange(desc(stat)) %>% filter(stat>0.2)
    cluster_stat = deframe(cluster_stat)


    cluster_stat2 = deframe(cluster_stat2)
    test = fgsea(pathways=pathways.hallmark, stats=cluster_stat, nperm=1000)
    ## my_pathways = reactomePathways(names(cluster_stat))
    ## test = fgsea(pathways=my_pathways, stats=cluster_stat, minSize=5, maxSize=500,
    ##              nperm=1000)
    gsea.results[[cluster]][['fgsea']] = test[order(test$padj), ]
    gsea.results[[cluster]][['fgsea']] = head(gsea.results[[cluster]][['fgsea']], 20)

    x <- enrichPathway(gene=names(cluster_stat2), pvalueCutoff=0.05, readable=T)
    gsea.results[[cluster]][['reactomepa']] = x
    ## dotplot(x, showCategory=15)
    y <- gsePathway(cluster_stat2, nPerm=10000,
                    pvalueCutoff=0.2,
                    pAdjustMethod="BH", verbose=FALSE)
    gsea.results[[cluster]][['reactomepa2']] = y
    ## viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=cluster_stat2)
    ## emapplot(y, color="pvalue")
    ## gseaplot(y, geneSetID = "R-HSA-69242")
    ## print(head(y))
    #z1 <- enrichGO(gene=names(cluster_stat2),
    #              ## universe=
    #              ont="CC",
    #              OrgDb=org.Hs.eg.db,
    #              pAdjustMethod="BH",
    #              pvalueCutoff=0.1,
    #              qvalueCutoff=0.1,
    #              readable=TRUE)
    #gsea.results[[cluster]][['go_cc']] = z1
    #z2 <- enrichGO(gene=names(cluster_stat2),
    #              ## universe=
    #              ont="MF",
    #              OrgDb=org.Hs.eg.db,
    #              pAdjustMethod="BH",
    #              pvalueCutoff=0.1,
    #              qvalueCutoff=0.1,
    #              readable=TRUE)
    #gsea.results[[cluster]][['go_mf']] = z2
    #z3 <- enrichGO(gene=names(cluster_stat2),
    #              ## universe=
    #              ont="BP",
    #              OrgDb=org.Hs.eg.db,
    #              pAdjustMethod="BH",
    #              pvalueCutoff=0.1,
    #              qvalueCutoff=0.1,
    #              readable=TRUE)
    #gsea.results[[cluster]][['go_bp']] = z3

    pdf(paste0('../results/gsea/', args$label, '_gsea_cluster', cluster, '.pdf'), width=18, height=8)
    par(mar=c(3, 12, 2, 3), font=2)
    topPathways = head(gsea.results[[cluster]][['fgsea']][,pathway], 8)
    plotGseaTable(pathways.hallmark[topPathways], cluster_stat, gsea.results[[cluster]][['fgsea']], 
                  gseaParam = 0.5)
    barplot(x, showCategory=20)
    ## emapplot(x)
    ## cnetplot(x, categorySize="pvalue", foldChange=geneList)
    dev.off()
    ## output gene set for each cluster in a lisa directory
    cat(names(cluster_stat), sep='\n', file=paste0('../results/gsea/lisa/', args$label, '_', cluster, 'cluster_genes.symbols'))
}

saveRDS(gsea.results, file=paste0('../results/gsea/', args$label, '_gsea_go_enrichment.rds'))
