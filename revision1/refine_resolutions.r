library(Seurat)
library(clustree)
library(tidyverse)

ds = readRDS('/data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds')

ds = FindClusters(ds, resolution=seq(0.1, 2, 0.1))

pdf("MSK72117_resolutions.pdf", height=13.5, width=9)
clustree(ds, prefix='RNA_snn_res.')
dev.off()

## look at above figure
## Claudia suggestion 0.9
ds = FindClusters(ds, resolution=0.9)
source("DEGs_seurat3_sara.R")
Idents(ds) <- ds@meta.data$seurat_clusters <- ds@meta.data$RNA_snn_res.0.9
allcluster <- names((ds@meta.data$seurat_clusters) %>% table())

clusterde <- list()
for (i in allcluster) {
    print(i)
    print(allcluster[-(as.integer(i)+1)])
    de.up <- get_upregulated_genes_in_clusters(ds, i, allcluster[-(as.integer(i)+1)])
    if (nrow(de.up) > 0) {
        de.up$cluster <- i
    }
    clusterde[[i]] <- de.up
}

ds.de <- do.call('rbind', clusterde)
ds.de <- subset(ds.de, p.adjusted <= 0.01 & (fg_fraction >= 0.1 | bg_fraction >= 0.1))
write.csv(ds.de, file='MSK72117_res0.9.csv', quote=F)

#ds.de = FindAllMarkers(ds, only.pos=T, min.pct=0.1,
#                       test.use='MAST', min.diff.pct=0.1,
#                       random.seed=100, logfc.threshold = 0.2)
#ds.de$enrichment = ds.de$pct.1 - ds.de$pct.2
#write.csv(ds.de, file='MSK72117_res0.9.csv', quote=F)

## 1.4 might be a good resolution
#ds = FindClusters(ds, resolution=1.4)
#ds.de = FindAllMarkers(ds, only.pos=T, min.pct=0.1,
#                       test.use='MAST', min.diff.pct=0.1,
#                       random.seed=100, logfc.threshold = 0.2)
#ds.de$enrichment = ds.de$pct.1 - ds.de$pct.2
#write.csv(ds.de, file='MSK72117_res1.4.csv', quote=F)

#ds = FindClusters(ds, resolution=2.0)
#ds.de = FindAllMarkers(ds, only.pos=T, min.pct=0.1,
#                       test.use='MAST', min.diff.pct=0.1,
#                       random.seed=100, logfc.threshold = 0.2)
#ds.de$enrichment = ds.de$pct.1 - ds.de$pct.2
#write.csv(ds.de, file='MSK72117_res2.0.csv', quote=F)
