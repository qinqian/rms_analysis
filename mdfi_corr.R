library(Seurat)
objs = Sys.glob("/data/langenau/human_rms_pdxs/seurat_objects/*")

results = data.frame(tumor=as.character(),
                     cell =as.character(),
                     gene=as.character(),
                     cor=as.double(),
                     pval=as.double(),
                     cell=as.double(),
                     cor0=as.double(),
                     pval0=as.double(),
                     cell0=as.double())

metalabels = c("GROUND", "Hypoxia", "EMT",  "G1S", "UNASSIGNED", "G2M",  "MUSCLE", "INTERFERON", "PROLIF", "Histones")
colortab = read.table('color_table.xls', sep='\t', header=T, check.names=F)

for (o in objs) {
    obj = readRDS(o)
    clusterlabel = sub('.rds', '', unlist(strsplit(basename(o), '_'))[3])
    cluster = as.character(colortab[, clusterlabel])
    labels = na.omit(metalabels[match(cluster, metalabels)])
    levels(obj@meta.data$RNA_snn_res.0.8) = labels
    for (gene in c('MYOD1', 'MYOG', 'MYF5')) {
        for (cell in levels(obj@meta.data$RNA_snn_res.0.8)) {
            if (sum(rownames(obj$RNA@data) %in% gene)>=1) {
                mdfi = obj$RNA@data['MDFI', obj@meta.data$RNA_snn_res.0.8==cell]
                other = obj$RNA@data[gene, obj@meta.data$RNA_snn_res.0.8==cell]
                vals = cor.test(mdfi,
                                other, method='spearman')
                allgene = cbind(mdfi, other)
                print(sum(rowSums(allgene)!=0))
                allgene = allgene[rowSums(allgene)!=0, ]
                mdfi0 = allgene[, 1]
                other0 = allgene[, 2]
                vals0 = cor.test(mdfi0,
                                 other0, method='spearman')
                newRow <- data.frame(tumor = basename(o),
                                     cell = cell,
                                     gene = gene,
                                     cor=vals$estimate,
                                     pval=vals$p.value,
                                     cell=sum(obj@meta.data$RNA_snn_res.0.8==cell),
                                     cor0=vals0$estimate,
                                     pval0=vals0$p.value,
                                     cell0=nrow(allgene))
            }else{
                newRow <- data.frame(tumor = basename(o),
                                     cell = cell,
                                     gene = gene,
                                     cor=NA,
                                     pval=NA,
                                     cell=NA,
                                     cor0=NA,
                                     pval0=NA,
                                     cell0=NA)
            }
            results = rbind(results, newRow)
        }
    }
}

## write.table(results, sep='\t', quote=F, file='mdfi_cor.xls', row.names=F)
write.table(results, sep='\t', quote=F, file='mdfi_cor_cells.xls', row.names=F)
