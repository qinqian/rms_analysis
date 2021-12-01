library(Seurat)
library(plot3D)

muscle = scan("../../../projects/01_sc_rms/final_annotations/gene_modules/Muscle.txt", what="")
emt = scan("../../../projects/01_sc_rms/final_annotations/gene_modules/EMT.txt", what="")
prolif = scan("../../../projects/01_sc_rms/final_annotations/gene_modules/Prolif.txt", what="")

plotmarkers = function(x) {
    state = unlist(meta[as.character((x@meta.data)[1,1]), , drop=T])
    state = state[!((state=='') | (is.na(state)))]
    levels(x$RNA_snn_res.1) = state
    clusters = x$RNA_snn_res.1
    levels(clusters)[grepl('EMT', levels(clusters))] = 'EMT'
    cols = as.vector(clusters)
    cols[!(clusters%in%c("Muscle", "EMT", "Prolif", "Ground"))] = "Other"
    cols = factor(cols, levels=sort(c(unique(as.character(clusters[clusters%in%c("Muscle", "EMT", "Prolif", "Ground")])), "Other")))
    print(levels(cols))
    muscle_markers = Matrix::colMeans(x[muscle,]$RNA@data)
    emt_markers = Matrix::colMeans(x[emt,]$RNA@data)
    prolif_markers = Matrix::colMeans(x[prolif,]$RNA@data)
    print(levels(cols))
    named_cols = c("red", "purple", "blue", "brown", "grey")
    names(named_cols) = c("Muscle", "EMT", "Prolif", "Ground", "Other")
    scatter3D(muscle_markers, emt_markers, prolif_markers,
              xlab="Muscle", ylab="EMT", zlab="Prolif",
              xlim=c(0, 2.0), ylim=c(0, 1.3), zlim=c(0, 1.3),
              bty = "g", alpha=0.7, ticktype = "detailed",
              theta = 135, phi = 40, pch = 16, cex=0.5, main=as.character((x@meta.data)[1,1]), colvar=as.integer(cols) , col = named_cols[levels(cols)],
              colkey = list(at = seq(1, length(levels(cols))), side = 1, 
                            addlines = TRUE, length = 0.6, width = 0.4,
                            labels = levels(cols)))
}

meta = read.delim("LARRY2_cellstate_res1.txt")
rownames(meta) = meta[ ,1]
meta = meta[,-1]

reg1 <- readRDS('results/seurat_sara/Regular_1_seurat-object.rds')
reg2 <- readRDS('results/seurat_sara/Regular_2_seurat-object.rds')
diff3 <- readRDS('results/seurat_sara/Diff_3_seurat-object.rds')
diff4 <- readRDS('results/seurat_sara/Diff_4_seurat-object.rds')

reg1$seurat_clusters = reg1$RNA_snn_res.1
reg2$seurat_clusters = reg2$RNA_snn_res.1
diff3$seurat_clusters = diff3$RNA_snn_res.1
diff4$seurat_clusters = diff4$RNA_snn_res.1

pdf("EMT_muscle_markers_second_Larry.pdf", width=11, height=9)
par(mfrow=c(2, 2), cex=0.6)
tumors = lapply(list(reg1, reg2, diff3, diff4), function(tumor) {
    plotmarkers(tumor)
})
dev.off()
