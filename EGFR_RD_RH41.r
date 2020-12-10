library(Seurat)
library(patchwork)

rd = readRDS('../results/seurat_sara/RD_seurat-object.rds')
rh = readRDS('../results/seurat_sara/RH41_seurat-object.rds')

pdf('RD_RH41_EGFR.pdf', width=12, height=8)
p1.0=DimPlot(rd, group.by = 'RNA_snn_res.0.8')
p1=FeaturePlot(rd, 'EGFR')
p2.0=DimPlot(rh, group.by = 'RNA_snn_res.0.8')
p2=FeaturePlot(rh, 'EGFR')
print(p1.0+p1+p2.0+p2+plot_layout(ncol=2))
dev.off()

pdf('RD_RH41_FGFR4.pdf', width=12, height=8)
p1.0=DimPlot(rd, group.by = 'RNA_snn_res.0.8')
p1=FeaturePlot(rd, 'FGFR4')
p2.0=DimPlot(rh, group.by = 'RNA_snn_res.0.8')
p2=FeaturePlot(rh, 'FGFR4')
print(p1.0+p1+p2.0+p2+plot_layout(ncol=2))
dev.off()


final = c(sum(as.vector(rd['EGFR']$RNA@data)>0)/dim(rd)[2],
          sum(as.vector(rh['EGFR']$RNA@data)>0)/dim(rh)[2])

a = Sys.glob('/data/langenau/human_rms_pdxs/data_april/*Robj')
for (i in a) {
    load(i)
    print(i)
}

rm(a)
result = list()
rm(i)
for (i in  ls()) {
    if (grepl("AP", i)) {
        print(i)
        seu = UpdateSeuratObject(get(i))
        ## seu = as.loom(seu, filename=paste0(i, '.loom'), verbose=T)
        ## seu$close_all()
        result[[i]] = seu
        }
}


for (i in names(result)) {
    final <- c(final, sum(result[[i]]['EGFR']$RNA@data > 0)/dim(result[[i]])[2])
}

names(final) <- c("RD", 'RH41', names(result))

pdf('cellproportion_RMS_normal.pdf', height=6, width=12)
par(mar=c(3, 18, 2, 2))
barplot(final, beside=T, horiz=T, las=2)
title('RMS and normal cell proportion with EGFR expression')
dev.off()
