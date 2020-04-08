library(Seurat)
tumorl = Sys.glob(paste0('../figures/*_tumoronly_res0.8_umap.rds'))
tumorl = tumorl[-1]
tumorl = c('20082_recluster2_tumor_only.rds', tumorl)

tumor = lapply(tumorl, readRDS)

library(ggplot2)
pdf('primary_tumor_mylpf.pdf', width=12, height=12)
print(CombinePlots(plots=lapply(tumor, function(x) {
    FeaturePlot(x, features=c("MYLPF", "MYF5", "MYOG"), ncol=3)
}), ncol=1))
dev.off()
