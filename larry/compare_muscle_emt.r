library(Seurat)

diffmedia = readRDS('results/seurat_sara/Differentiated_seurat-object.rds')
regmedia = readRDS('results/seurat_sara/Regular_seurat-object.rds')

muscle = scan("../../projects/01_sc_rms/final_annotations/gene_modules/Muscle.txt", what="")
emt = scan("../../projects/01_sc_rms/final_annotations/gene_modules/EMT.txt", what="")
prolif = scan("../../projects/01_sc_rms/final_annotations/gene_modules/Prolif.txt", what="")

## pdf("EMT_muscle_markers_ARMS.pdf", width=15, height=9.5)
## par(mfrow=c(2, 3), cex=0.9)
## muscle_markers = Matrix::colMeans(MAST111[muscle,]$RNA@data)
## emt_markers = Matrix::colMeans(MAST111[emt,]$RNA@data)
## prolif_markers = Matrix::colMeans(MAST111[prolif,]$RNA@data)
## ## redundant points
## plot(muscle_markers, prolif_markers, pch=16, main="MAST111", xlim=c(-3, 3))
## points(-emt_markers, prolif_markers, pch=16, main="MAST111")
## muscle_markers = Matrix::colMeans(MAST111[muscle,]$RNA@data)
## emt_markers = Matrix::colMeans(MAST111[emt,]$RNA@data)
## plot(muscle_markers, emt_markers, pch=16, main="MAST111")
## muscle_markers = Matrix::colMeans(MAST139[muscle,]$RNA@data)
## emt_markers = Matrix::colMeans(MAST139[emt,]$RNA@data)
## plot(muscle_markers, emt_markers, pch=16, main="MAST139")
## muscle_markers = Matrix::colMeans(RH74[muscle,]$RNA@data)
## emt_markers = Matrix::colMeans(RH74[emt,]$RNA@data)
## plot(muscle_markers, emt_markers, pch=16, main="RH74")
## muscle_markers = Matrix::colMeans(MAST118[muscle,]$RNA@data)
## emt_markers = Matrix::colMeans(MAST118[emt,]$RNA@data)
## plot(muscle_markers, emt_markers, pch=16, main="MAST118")
## muscle_markers = Matrix::colMeans(mast95[muscle,]$RNA@data)
## emt_markers = Matrix::colMeans(mast95[emt,]$RNA@data)
## plot(muscle_markers, emt_markers, pch=16, main="MAST95")
## muscle_markers = Matrix::colMeans(MSK82489[muscle,]$RNA@data)
## emt_markers = Matrix::colMeans(MSK82489[emt,]$RNA@data)
## plot(muscle_markers, emt_markers, pch=16, main="MSK82489")
## dev.off()

library("plot3D")

plotmarkers = function(x) {
    ## print(as.character((x@meta.data)[1,1]))
    ## print(unlist(meta[as.character((x@meta.data)[1,1]), , drop=T]))
    state = unlist(meta[as.character((x@meta.data)[1,1]), , drop=T])
    state = state[!((state=='') | (is.na(state)))]
    ## print(state)
    levels(x$RNA_snn_res.0.8) = state
    clusters = x$RNA_snn_res.0.8
    cols = as.vector(clusters)
    cols[!(clusters%in%c("Muscle", "EMT", "Prolif", "Ground"))] = "Other"
    print(head(cols))
    cols = factor(cols, levels=c(unique(as.character(clusters[clusters%in%c("Muscle", "EMT", "Prolif", "Ground")])), "Other"))
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

meta = read.delim('LARRY_cellstate.txt', sep='\t', row.names=1, header=T,
                  check.names=F, stringsAsFactors=F)

pdf("EMT_muscle_markers_ARMS.pdf", width=11, height=5)
par(mfrow=c(1, 2), cex=0.6)
for (tumor in list(diffmedia, regmedia)) {
    print(as.character((tumor@meta.data)[1,1]))
    plotmarkers(tumor)
}
dev.off()
