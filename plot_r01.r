library(Seurat)
library(tidyr)
library(patchwork)

MAST111 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds')
MAST139 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds')
MSK74711 = readRDS('../results/seurat_sara/20191031_MSK74711_seurat-object.rds')

label = 'MAST139'
annotation <- read.delim('../final_annotations/Final_clusters.txt', row.names=1,
                         stringsAsFactors = F)

states = as.vector(unlist(as.vector(annotation[label, ])))
states = as.vector(states[!is.na(states)])
states = states[states!='']
levels(MAST139$RNA_snn_res.0.8) = states

states = as.vector(unlist(as.vector(annotation["MSK74711", ])))
states = as.vector(states[!is.na(states)])
states = states[states!='']
levels(MSK74711$RNA_snn_res.0.8) = states

states = as.vector(unlist(as.vector(annotation["MAST111", ])))
states = as.vector(states[!is.na(states)])
states = states[states!='']
levels(MAST111$RNA_snn_res.0.8) = states

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

pdf("UMAPS_MAST111_MAST139_MSK74711.pdf", width=18, height=5.5)
print(DimPlot(MAST111, cols=metacolors, group='RNA_snn_res.0.8')+DimPlot(MAST139, cols=metacolors, group='RNA_snn_res.0.8')+DimPlot(MSK74711, cols=metacolors, group='RNA_snn_res.0.8'))
dev.off()

pdf("Markers_MAST111_MAST139_MSK74711.pdf", width=18, height=18)
## FUNDC2: DC44
## POSTN1 and MYPLF cannot be found
print(FeaturePlot(MAST111, c("OGN", "MGP", "FUNDC2", "CHODL", "POSTN1", "MEOX2", "COL4A2", "COL5A3", "LAMB1", "EBF1", "LRRN1", "TSPAN33", "CKM", "MYPLF", "ACTN2", "TNNT2", "TNNC2", "MYLPF", "POSTN")))
print(FeaturePlot(MAST139, c("OGN", "MGP", "FUNDC2", "CHODL", "POSTN1", "MEOX2", "COL4A2", "COL5A3", "LAMB1", "EBF1", "LRRN1", "TSPAN33", "CKM", "MYPLF", "ACTN2", "TNNT2", "TNNC2", "MYLPF", "POSTN")))
print(FeaturePlot(MSK74711, c("OGN", "MGP", "FUNDC2", "CHODL", "POSTN1", "MEOX2", "COL4A2", "COL5A3", "LAMB1", "EBF1", "LRRN1", "TSPAN33", "CKM", "MYPLF", "ACTN2", "TNNT2", "TNNC2", "MYLPF", "POSTN")))
dev.off()
