library(glue)
## library("ggVennDiagram")
library(rlang)
library("ggvenn")
library(patchwork)
library(Seurat)

## Sys.glob('~/langenau/projects/02_atlas/normal_muscle_cellbrowser/skeletal-muscle/fetal/*wk9*myogenic/meta.tsv')
emt = scan('~/langenau/projects/01_sc_rms/final_annotations/gene_modules/EMT.txt', what='')

pdf("emt_normal_venn.pdf", width=18, height=5)
plots=list()
allemt = list()
for (wk in c("wk9", "wk12-14", "wk17-18")) {
    skm = Sys.glob(glue('~/langenau/projects/02_atlas/normal_muscle_cellbrowser/skeletal-muscle/fetal/fetal-{wk}-myogenic/SkM*.tsv.gz'))
    skm = read.table(skm, sep='\t', header=T, stringsAsFactors=F)
    x = list()
    x[[as_name(wk)]]=skm$symbol
    x[['EMT']] = emt
    allemt[[as_name(wk)]]=skm$symbol
    allemt[['EMT']]=emt
    ## print(ggVennDiagram(x, label_alpha = 0))
    ## print(ggvenn(x))
    plots[[wk]] = ggvenn(x)
}
print(wrap_plots(plots))
dev.off()

muscle = scan('~/langenau/projects/01_sc_rms/final_annotations/gene_modules/Muscle.txt', what='')

pdf("muscle_normal_venn.pdf", width=20, height=5)
plots=list()
allemt = list()
x = list()
for (wk in c("wk6-7", "wk9", "wk12-14")) {
    skm = Sys.glob(glue('~/langenau/projects/02_atlas/normal_muscle_cellbrowser/skeletal-muscle/*/*-{wk}-myogenic/*MC*.tsv.gz'))[1]
    print(skm)
    skm = read.table(skm, sep='\t', header=T, stringsAsFactors=F)
    x[[as_name(wk)]]=skm$symbol
    x[['muscle']] = muscle
    allemt[[as_name(wk)]]=skm$symbol
    allemt[['muscle']]=muscle
    print(intersect(muscle, skm$symbol))
    ## plots[[wk]] = ggvenn(x)
}
print(ggvenn(x))
## print(wrap_plots(plots))
dev.off()


annotation = read.delim('~/langenau/projects/01_sc_rms/final_annotations/Final_clusters.txt', sep='\t', row.names=1, header=T, check.names=F, stringsAsFactors=F)
pdxs = Sys.glob('~/langenau/projects/01_sc_rms/data/seurat_obj/*rds')[1:10]
pdxs = c(pdxs, '~/langenau/projects/01_sc_rms/results/seurat_sara/20191031_MSK74711_seurat-object.rds',
         '~/langenau/projects/01_sc_rms/results/seurat_sara/MAST118_seurat-object.rds',
         '~/langenau/projects/01_sc_rms/results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds',
         '~/langenau/projects/01_sc_rms/results/seurat_sara/MAST139_1cells_seurat-object.rds')

pdxs.objs = lapply(pdxs, readRDS)
labels = unlist(lapply(pdxs.objs, function(x) {
    levels(x$orig.ident[1])
}))

labels[9] = 'MAST85-1'
labels[11] = 'MSK74711'
labels[13] = 'MSK72117'
labels[14] = 'MAST139-1'
labels[10] = 'RH74-10'

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

results = list()
for (i in seq_along(labels)) {
    states = unlist(as.vector(annotation[labels[i], ]))
    states = states[states!='']
    print(states)
    pdxs.objs[[i]]$seurat_clusters = pdxs.objs[[i]]$RNA_snn_res.0.8
    levels(pdxs.objs[[i]]$RNA_snn_res.0.8) = states
    Idents(pdxs.objs[[i]]) = pdxs.objs[[i]]$seurat_clusters
    results[[labels[i]]] = table(pdxs.objs[[i]]$RNA_snn_res.0.8)
}

markers1 = intersect(allemt$wk9, allemt$EMT)[1:20]
markers2 = setdiff(allemt$wk9, allemt$EMT)[1:20]

pdf("wk9_intersect_setdiff_markers.pdf", width=19.5, height=14.5)
for (i in seq_along(labels)) {
    plots = list()
    plots[['clusters']]=DimPlot(pdxs.objs[[i]], group.by = 'RNA_snn_res.0.8', cols=metacolors)+ggtitle(labels[i])
    for (d in markers1[markers1 %in% rownames(pdxs.objs[[i]])]) {
        plots[[d]]=FeaturePlot(pdxs.objs[[i]], d)
    }
    for (w in markers2[markers2 %in% rownames(pdxs.objs[[i]])]) {
        plots[[w]]=FeaturePlot(pdxs.objs[[i]], w)
    }
    print(wrap_plots(plots))
}
dev.off()


markers1 = intersect(allemt[['wk12-14']], allemt$EMT)[1:17]
markers2 = setdiff(allemt[['wk12-14']], allemt$EMT)[1:20]

pdf("wk12-14_intersect_setdiff_markers.pdf", width=19.5, height=14.5)
for (i in seq_along(labels)) {
    plots = list()
    plots[['clusters']]=DimPlot(pdxs.objs[[i]], group.by = 'RNA_snn_res.0.8', cols=metacolors)+ggtitle(labels[i])
    for (d in markers1[markers1 %in% rownames(pdxs.objs[[i]])]) {
        plots[[d]]=FeaturePlot(pdxs.objs[[i]], d)
    }
    for (w in markers2[markers2 %in% rownames(pdxs.objs[[i]])]) {
        plots[[w]]=FeaturePlot(pdxs.objs[[i]], w)
    }
    print(wrap_plots(plots))
}
dev.off()
