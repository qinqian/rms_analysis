
library(Seurat)
library(foreach)
library(doMC)

obj.list = list()
for (rds in Sys.glob('/data/langenau/human_rms_pdxs/seurat_objects/*rds')) {
    obj.list[[basename(rds)]] = readRDS(rds)
    Idents(obj.list[[basename(rds)]]) = obj.list[[basename(rds)]]@meta.data$seurat_clusters = obj.list[[basename(rds)]]$RNA_snn_res.0.8
}

registerDoMC(10)
obj.de = foreach(i=names(obj.list)) %dopar% {
       FindAllMarkers(obj.list[[i]], logfc.threshold = 0.05,
                      test.use = 'MAST', min.pct = 0.05, only.pos=T)
}

names(obj.de) <- c("MAST111", "MAST139", "MAST35", "MAST39", "MAST85",
                   "MAST95",  "MSK82489", "RH74",  "MAST85.1cell", "RH74.10cells")

saveRDS(obj.de, 'human_10samples_allde.rds')

meta = read.delim('color_table.xls', stringsAsFactors=F)

## system('mkdir -p lisa_mergedcluster_gene_set')
system('mkdir -p lisa_mergedcluster_gene_set2')
for (n in names(obj.de)) {
    x = obj.de[[n]]
    x$state = meta[match(as.vector(x$cluster), meta[, 1]), n]
    for (s in unique(x$state)) {
        sx = subset(x, state==s)
        sx$enrich = sx$pct.1 - sx$pct.2
        ## sx = sx[order(sx$enrich, sx$avg_logFC, decreasing=T),][1:300,]
        sx = sx[order(sx$avg_logFC, sx$enrich, decreasing=T),] # [1:500,]
        sx = subset(sx, avg_logFC>=0.2 & sx$enrich>=0.1 & p_val_adj <= 0.05)
        print(dim(sx))
        out = sx$gene
        ## cat(unique(out), sep='\n', file=paste0('lisa_mergedcluster_gene_set/', n, '_', s))
        cat(unique(out), sep='\n', file=paste0('lisa_mergedcluster_gene_set2/', n, '_', s))
    }
}
