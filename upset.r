devtools::install_github("hms-dbmi/UpSetR")
devtools::install_github("mw201608/SuperExactTest")
library(devtools)
install_github("jokergoo/ComplexHeatmap")
BiocManager::install("clusterProfiler")
install_github("YuLab-SMU/clusterProfiler")
devtools::install_github("thomasp85/patchwork")


## library(UpSetR)
library(ComplexHeatmap)
library(clusterProfiler)
library(SuperExactTest)
library(patchwork)

gmt = read.gmt('~/langenau/01_rms_projects/02_human/data/msigdb/c5.all.v7.0.symbols.gmt')
neural = unique(subset(gmt, ont=='GO_NEUROGENESIS')[,2])

## neural = gmt[gmt[,1]%in%c('GO_NEUROGENESIS', 'GO_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION'), ][,2]

gsf = Sys.glob('../final_annotations/gene_modules/*')
gs = lapply(gsf, scan, what='')
gs = lapply(gs, unique)

names(gs) = gsub('.txt', '', basename(gsf))
gsm = list_to_matrix(gs, universal_set = unique(unlist(gs)))

## see mode: https://jokergoo.github.io/ComplexHeatmap-reference/book/08-upset_files/figure-html/unnamed-chunk-7-1.png
set.seed(123)
## allm = make_comb_mat(gsm, mode='distinct') # default
## a bug for intersect and union mode https://github.com/jokergoo/ComplexHeatmap/issues/382

allm = make_comb_mat(gsm, mode='intersect')

## y1 <- enricher(gs[['UNIQUE_1']],
##                ## TERM2GENE = gmt[gmt[,1]%in%c('GO_NEUROGENESIS', 'GO_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION'), ],
##                TERM2GENE = gmt[gmt[,1]%in%c('GO_NEUROGENESIS'), ],
##                minGSSize=1, pvalueCutoff=1, universe=unique(gmt[,2]),
##                qvalueCutoff=1)
## y2 <- enricher(gs[['UNIQUE_1']],
##                ## TERM2GENE = data.frame(module='UNIQUE_5', gene=neural),
##                TERM2GENE = gmt,
##                minGSSize=1, pvalueCutoff=1,
##                qvalueCutoff=1)

## df = data.frame(gene.in.intersest=c(length(intersect(gs[['UNIQUE_1']], neural)), length(gs[['UNIQUE_1']])-length(intersect(gs[['UNIQUE_1']], neural))),
##                 gene.not.interest=c(length(neural)-length(intersect(gs[['UNIQUE_1']], neural)),
##                                     length(unique(gmt[,2]))-length(gs[['UNIQUE_1']])-(length(neural)-length(intersect(gs[['UNIQUE_1']], neural)))))
## rownames(df) = c('GO:NEURO', 'NON-NEURO')
## fisher.test(df, alternative='greater')

## df = data.frame(gene.in.intersest=c(length(intersect(gs[['UNIQUE_3']], neural)), length(gs[['UNIQUE_3']])-length(intersect(gs[['UNIQUE_3']], neural))),
##                 gene.not.interest=c(length(neural)-length(intersect(gs[['UNIQUE_3']], neural)),
##                                     length(unique(gmt[,2]))-length(gs[['UNIQUE_3']])-(length(neural)-length(intersect(gs[['UNIQUE_3']], neural)))))
## rownames(df) = c('GO:NEURO', 'NON-NEURO')
## fisher.test(df, alternative='greater')

## df = data.frame(gene.in.intersest=c(length(intersect(gs[['UNIQUE_6']], neural)), length(gs[['UNIQUE_6']])-length(intersect(gs[['UNIQUE_6']], neural))),
##                 gene.not.interest=c(length(neural)-length(intersect(gs[['UNIQUE_6']], neural)),
##                                     length(unique(gmt[,2]))-length(gs[['UNIQUE_6']])-(length(neural)-length(intersect(gs[['UNIQUE_6']], neural)))))
## rownames(df) = c('GO:NEURO', 'NON-NEURO')
## fisher.test(df, alternative='greater')

pdf("Supp_all_modules_upset.pdf")
print(UpSet(allm[(comb_size(allm)>1 & comb_degree(allm)<=2)]))
## print(upset(as.data.frame(gsm), order.by = "freq", ## weird results, ARMS.core 6 genes
##             empty.intersections = "on", sets.bar.color = "#56B4E9"))
dev.off()

result = supertest(gs, n=length(unique(unlist(gs))))

pdf("Supp_all_modules_upset_test.pdf", width=13, height=13)
plot(result, Layout="landscape", sort.by="size",
     keep=FALSE,
     bar.split=c(150, 1555),
     show.elements=TRUE,
     elements.cex=0.5,
     ## bar.area.range=0.1,
     min.intersection.size=1,
     elements.list=subset(summary(result)$Table, Observed.Overlap <= 15),
     ## show.expected.overlap=TRUE,
     ## expected.overlap.style="hatchedBox",
     minMinusLog10PValue=1.3)
     ## color.expected.overlap='red')
dev.off()

## gs[['GO_NEUROGENESIS']] = neural

de = c('../results/seurat_sara/20190624_seurat-object_MAST111_SCT_res0.8.xls', '../results/seurat_sara/20190624_seurat-object_MAST139_SCT_res0.8.xls', '../results/seurat_sara/20190624_seurat-object_MAST35_SCT_res0.8.xls', '../results/seurat_sara/20190624_seurat-object_MAST39_SCT_res0.8.xls', '../results/seurat_sara/20190624_seurat-object_MAST85_SCT_res0.8.xls', '../results/seurat_sara/20190624_seurat-object_MAST95_SCT_res0.8.xls', '../results/seurat_sara/20190624_seurat-object_MSK82489_SCT_res0.8.xls', '../results/seurat_sara/20190624_seurat-object_RH74_SCT_res0.8.xls', '../results/seurat_sara/20191031_MSK72117tencell_sara_v4_res0.8.xls', '../results/seurat_sara/20191031_MSK74711_seurat-object_sara_v4_res0.8.xls', '../results/seurat_sara/MAST118_seurat-object_sara_v4_res0.8.xls')

cols = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1,
                  header=T,
                  check.names=F, stringsAsFactors=F)
labels = c("MAST111", "MAST139", "MAST35", "MAST39",
           "MAST85", "MAST95", "MSK82489", "RH74",
           "MSK72117", "MSK74711", "MAST118")
deal = list()
for (i in seq_along(labels)) {
    dea = read.table(de[i], header=T, stringsAsFactors = F)
    clusters = unique(dea$cluster)
    dea = lapply(clusters, function(x) {
        unique(dea[dea$cluster == x, ]$genename)
    })
    names(dea) = paste0(labels[i], '_cluster', clusters)
    deal = append(deal, dea)
}

deal[['GO_NEUROGENESIS']] = neural

dfs = as.data.frame(list_to_matrix(deal))

targetm = deal[#c('GO_NEUROGENESIS',
    c(
                 'MSK74711_cluster8',
                 'MSK72117_cluster4',
                 'MAST118_cluster12',
        'MAST95_cluster8')]
targetm[['MSK74711_cluster3/5']] = unique(c(deal[['MSK74711_cluster3']], deal[['MSK74711_cluster5']]))

## supertest is too strict for gene ontology analysis
## since universe gene set is restricted
resultm = supertest(targetm,  n=length(unique(unlist(targetm))))
pdf("Supp_all_neural_upset_test.pdf", height=6, width=8)
plot(resultm, Layout="landscape", sort.by="size",
     keep=FALSE,
     ## bar.split=c(630, 1550),
     show.elements=TRUE,
     elements.cex=0.5,
     ## bar.area.range=0.1,
     min.intersection.size=1,
     margin=c(3, 9, 3, 5),
     elements.list=subset(summary(resultm)$Table, Observed.Overlap <= 15),
     ## show.expected.overlap=TRUE,
     ## expected.overlap.style="hatchedBox",
     minMinusLog10PValue=1.3)
dev.off()

## targetm = deal[#c('GO_NEUROGENESIS',
##     c(
##                  'MSK74711_cluster8',
##                  'MSK72117_cluster4',
##                  'MAST118_cluster12',
##                  'MAST95_cluster8')]

test1 = lapply(targetm, function(x) {
    result = enricher(x, TERM2GENE=gmt, minGSSize = 5, maxGSSize = 3000)
    head(result@result[order(result@result$qvalue), ], 500)
})

test1neural = lapply(test1, function(x)
    ## x[grep("NEURO", rownames(x))[1], ]
    x[grep("GO_NEUROGENESIS", rownames(x)), ]$p.adjust
)

## use enricher for gene ontology analysis
## geneset = c('GO_NEUROGENESIS',
##             'MSK74711_cluster8',
##             'MSK72117_cluster4',
##             'MAST118_cluster12',
##             'MAST95_cluster8'
##             )
## targetm = deal[geneset]

geneset = c('GO_NEUROGENESIS',
            'MSK74711_cluster8',
            'MSK72117_cluster4',
            'MAST118_cluster12',
            'MAST95_cluster8',
            "MSK74711_cluster3/5"
            )
targetm[['GO_NEUROGENESIS']] = neural

allf = make_comb_mat(list_to_matrix(targetm)[, geneset], mode='intersect')

## pdf("Supp_all_modules_upset.pdf")
## print(UpSet(allf[(comb_size(allf)>1 & comb_degree(allf)==2)]))
## ## print(upset(as.data.frame(gsm), order.by = "freq", ## weird results, ARMS.core 6 genes
## ##             empty.intersections = "on", sets.bar.color = "#56B4E9"))
## dev.off()

library(RColorBrewer)
cols = colorRampPalette(c('yellow', 'red'))(4)
combs = comb_name(allf)[stringr::str_count(comb_name(allf), '1') == 2 & (stringr::str_locate(comb_name(allf), '1')[,1]==1)]

targetf = allf[geneset, combs]

name = set_name(allf)[2:6][stringr::str_locate(substr(combs, 2, 6), '1')[,1]][order(comb_size(targetf))]

library(circlize)
col_fun = colorRamp2(c(0, 3, 5), c("white", "white", "red"))

pdf("Supp_neural_rms_clusters.pdf")
## %v% use the first heatmap order to align the second one
## that is name[order(comb_name(targetf, readable=T))]
pvals = -log10(unlist(test1neural))
UpSet(targetf, pt_size=unit(5, "mm"), lwd=3)%v%
    Heatmap(rbind(pvals), name = "-log10(p-value)",
            height = unit(6, "mm"), col=col_fun)
dev.off()

## print(upset(dfs,
##             sets=c('GO_NEUROGENESIS',
##                    'MSK74711_cluster8',
##                    'MSK72117_cluster4',
##                    'MAST118_cluster12',
##                    'MAST95_cluster8'),
            ## grep('MSK74711', colnames(dfs), value=T),
            ## grep('MSK72117', colnames(dfs), value=T)),
            ## grep('MAST118', colnames(dfs), value=T),
            ## grep('MAST95', colnames(dfs), value=T)),
            ## sets.bar.color = "#56B4E9", order.by = "freq", nsets=50)) # empty.intersections = "on"
            ## grep(, colnames(dfs), value=T)),
            ## grep('MAST118', colnames(dfs), value=T),
            ## grep('MAST95', colnames(dfs), value=T))

## for (i in c("MSK74711", 'MSK72117', "MAST95", "MAST118")) {
## print(upset(dfs,
##             sets=c('Neural',
##                    grep(i, colnames(dfs), value=T)),
##             sets.bar.color = "#56B4E9", order.by = "freq", nsets=50)) # empty.intersections = "on"
## }
