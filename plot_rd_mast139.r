library(Seurat)
library(tidyr)
library(patchwork)

rd = readRDS('../results/seurat_sara/RD_seurat-object.rds')
annotation <- read.delim('../final_annotations/cellline_annotation.txt', row.names=1,
                         stringsAsFactors = F)
label = 'RD'
states = as.vector(unlist(as.vector(annotation[label, ])))
states = as.vector(states[!is.na(states)])
states = states[states!='']
levels(rd$RNA_snn_res.0.8) = states

MAST139 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds')
label = 'MAST139'
annotation <- read.delim('../final_annotations/Final_clusters.txt', row.names=1,
                         stringsAsFactors = F)
states = as.vector(unlist(as.vector(annotation[label, ])))
states = as.vector(states[!is.na(states)])
states = states[states!='']
levels(MAST139$RNA_snn_res.0.8) = states

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

pdf("RD_MAST139_markers.pdf", width=12, height=10)
## THY1: CD90
## Ki-67: MKI67
DimPlot(rd, cols=metacolors, group='RNA_snn_res.0.8')
FeaturePlot(rd, c('THY1', 'CD44', 'LRRN1', 'TSPAN33', 'NDRG1', 'EGFR', 'TNNT3', 'MX1', 'MKI67'), ncol=3)
DimPlot(MAST139, cols=metacolors, group='RNA_snn_res.0.8')
FeaturePlot(MAST139, c('THY1', 'CD44', 'LRRN1', 'TSPAN33', 'NDRG1', 'EGFR', 'TNNT3', 'MX1', 'MKI67'), ncol=3)
dev.off()

write.csv(table(MAST139$RNA_snn_res.0.8), 'MAST139_cellprop.csv')
write.csv(table(rd$RNA_snn_res.0.8), 'RD_cellprop.csv')
