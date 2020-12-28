library(readr)

primary.tumors <- lapply(c('20082_recluster2_tumor_only.rds', # '../figures/20082_hg19_premrna_tumoronly_res0.8_umap.rds',
                           '../figures/20696_hg19_tumoronly_res0.8_umap.rds',
                           '../figures/21202_hg19_premrna_tumoronly_res0.8_umap.rds',
                           '../figures/29806_hg19_premrna_tumoronly_res0.8_umap.rds'), readRDS)

allpdx.meta2 = do.call('rbind', list(primary.tumors[[1]]@meta.data,
                                     primary.tumors[[2]]@meta.data,
                                     primary.tumors[[3]]@meta.data,
                                     primary.tumors[[4]]@meta.data
                                     ))

allpdx.meta2$orig.ident = gsub('_hg19_premrna', '', allpdx.meta2$orig.ident)

dfs = list()
for (i in Sys.glob('../results/doublets/2*_doublet_doublet.csv')) {
    dfs[[i]] = read.csv(i)
    dfs[[i]]$cells = gsub('\\-\\d', '', dfs[[i]]$cells)
    dfs[[i]]$orig.ident = unlist(strsplit(basename(i), '_'))[1]
    ## print(i)
    ## print(table(df[,4])/nrow(df))
}

dfs = do.call('rbind', dfs)
length(intersect(paste0(dfs$cells, dfs$orig.ident), paste0(rownames(allpdx.meta2), allpdx.meta2$orig.ident)))
print(table(dfs[match(paste0(rownames(allpdx.meta2), allpdx.meta2$orig.ident), paste0(dfs$cells, dfs$orig.ident)), 'prediction']))

MAST111 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds')
MAST139 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds')
MAST39 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST39.rds')
MAST95 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST95.rds')
MAST85 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST85.rds')
MSK82489 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MSK82489.rds')
MAST35 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST35.rds')
## new added samples
MSK72117_10cells = readRDS('/data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds')
MAST118 = readRDS('/data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/MAST118_seurat-object.rds')
MSK74711 = readRDS('/data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/20191031_MSK74711_seurat-object.rds')

all.obj = Reduce('merge', list(MAST111, MAST139, MAST39, MAST95, MAST85, MSK82489, MAST35,
                               MSK72117_10cells, MAST118, MSK74711))

allpdx.metadata = (all.obj@meta.data)
allpdx.metadata$orig.ident[allpdx.metadata$orig.ident=="20191031_MSK72117tencell"] = 'MSK72117'
allpdx.metadata$orig.ident[allpdx.metadata$orig.ident=="20191031_MSK74711"] = 'MSK74711'

dfs = list()
labels <- c("MAST111", "MAST139", "MAST35", "MAST71", "MAST85", "MSK72117", "MSK82489", "RH74-10", "RH74", "MAST39", "MAST95", "MAST85-1", "MSK72117", "MSK74711", "C12SC1", "C12SC2")
n = 1
for (i in Sys.glob('../results/*_doublet_doublet.csv')) {
    print(i)
    dfs[[i]] = read.csv(i)
    dfs[[i]]$cells = gsub('\\-\\d*', '', dfs[[i]]$cells)
    dfs[[i]]$orig.ident = labels[n]
    print(i)
    print(table(dfs[[i]][,4])/nrow(dfs[[i]]))
    n <- n+1
}
dfs = do.call('rbind', dfs)

rownames(allpdx.metadata) = paste0(gsub("_\\d", "", gsub("_\\d_1_1_1_1", "", rownames(allpdx.metadata))), allpdx.metadata$orig.ident)

length(intersect(paste0(dfs$cells, dfs$orig.ident), rownames(allpdx.metadata)))

print(table(dfs[match(paste0(rownames(allpdx.metadata)), paste0(dfs$cells, dfs$orig.ident)), 'prediction']))
