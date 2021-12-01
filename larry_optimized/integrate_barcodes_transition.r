set.seed(100)
annotation = read.delim('LARRY2_cellstate_res1.txt', sep='\t', row.names=1, header=T,
                        check.names=F, stringsAsFactors=F)

metacolors <- c(rgb(119, 62, 20, maxColorValue = 255),
                rgb(236, 133, 40, maxColorValue = 255),
                rgb(59, 22, 115, maxColorValue = 255),
                rgb(59, 22, 115, maxColorValue = 255),
                rgb(52, 101, 252, maxColorValue = 255),
                rgb(242, 242, 242, maxColorValue = 255),
                rgb(52, 101, 252, maxColorValue = 255),
                rgb(225, 39, 39, maxColorValue = 255),
                rgb(72, 159,  75, maxColorValue = 255),
                rgb(20, 64, 216, maxColorValue = 255),
                rgb(226, 75, 143, maxColorValue = 255),
                rgb(158, 60, 200, maxColorValue = 255),
                rgb(241, 250, 100, maxColorValue = 255),
                rgb(0, 255, 253, maxColorValue = 255))
metalabels <- c("Ground", "Hypoxia", "iEMT", "tEMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'UPR', "Neural")
names(metacolors) <- metalabels

library(Seurat)
library(reticulate)
library(patchwork)
library(ggplot2)

reg1 <- readRDS('results/seurat_sara/Regular_1_seurat-object.rds')
reg2 <- readRDS('results/seurat_sara/Regular_2_seurat-object.rds')
diff3 <- readRDS('results/seurat_sara/Diff_3_seurat-object.rds')
diff4 <- readRDS('results/seurat_sara/Diff_4_seurat-object.rds')

reg1$seurat_clusters = reg1$RNA_snn_res.1
reg2$seurat_clusters = reg2$RNA_snn_res.1
diff3$seurat_clusters = diff3$RNA_snn_res.1
diff4$seurat_clusters = diff4$RNA_snn_res.1

Idents(reg1) = reg1$seurat_clusters2 = reg1$seurat_clusters
states = unlist(as.vector(annotation["Regular_1", ]))
states = states[states!='']
states = states[!is.na(states)]
levels(reg1$seurat_clusters) = states
levels(Idents(reg1)) = states

Idents(reg2) = reg2$seurat_clusters2 = reg2$seurat_clusters
states = unlist(as.vector(annotation['Regular_2', ]))
states = na.omit(states[states!=''])
levels(reg2$seurat_clusters) = states
levels(Idents(reg2)) = states

Idents(diff3) = diff3$seurat_clusters2 = diff3$seurat_clusters
states = unlist(as.vector(annotation['Diff_3', ]))
states = na.omit(states[states!=''])
levels(diff3$seurat_clusters) = states
levels(Idents(diff3)) = states

Idents(diff4) = diff4$seurat_clusters2 = diff4$seurat_clusters
states = unlist(as.vector(annotation['Diff_4', ]))
states = na.omit(states[states!=''])
levels(diff4$seurat_clusters) = states
levels(Idents(diff4)) = states

integration <- function(x, y, z, zz, label, method='seurat') {
    seurat.pseudo.list = list(x, y, z, zz)
    for (i in 1:length(seurat.pseudo.list)) {
        print(i)
        seurat.pseudo.list[[i]] <- NormalizeData(seurat.pseudo.list[[i]], verbose = FALSE)
        seurat.pseudo.list[[i]] <- FindVariableFeatures(seurat.pseudo.list[[i]],
                                                        selection.method = "vst", 
                                                        nfeatures = 2000, verbose = FALSE)
    }
    names(seurat.pseudo.list) <- label
    if (method == 'seurat') {
        anchors <- FindIntegrationAnchors(object.list = seurat.pseudo.list, dims = 1:30)
        integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
        DefaultAssay(integrated) <- "integrated"
        integrated <- ScaleData(integrated, verbose = FALSE)
        integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
        integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
        ## integrated <- FindClusters(integrated, resolution = 1.0)
    } else if (method == 'conos') {
        seurat.pseudo.list <- lapply(seurat.pseudo.list, function(x) {
            ScaleData(x) %>% RunPCA(verbose=F)
        })
        pseudo.con <- Conos$new(seurat.pseudo.list)
        pseudo.con$buildGraph(k=15, k.self=5, space="PCA", ncomps=30, n.odgenes=2000, matching.method='mNN',
                              metric = 'angular', score.component.variance=T, verbose=T)
        pseudo.con$findCommunities()
        pseudo.con$embedGraph()
        integrated = as.Seurat(pseudo.con)
    }
    integrated
}

integrated.seurat4 = integration(reg1, reg2, diff3, diff4, 
                                 c('Regular_1', 'Regular_2', 'Diff_3', 'Diff_4'))


## saveRDS(integrated.seurat4, "second_Larry_integrative_res1.rds")
saveRDS(integrated.seurat4, "second_Larry_integrative_res1_labelled.rds")
saveRDS(integrated.seurat4, "second_Larry_res0.4.rds")

## integrated.seurat4 = readRDS("second_Larry_res0.4.rds")

## integrated.seurat4 = readRDS("second_Larry_integrative_res1.rds")
integrated.seurat4 = readRDS("second_Larry_integrative_res1_labelled.rds")
DefaultAssay(integrated.seurat4) <- "RNA"

write.table(integrated.seurat4@reductions$umap@cell.embeddings, file='second_Larry_integrative_res1_labelled_umap.xls', sep='\t', quote=F)
write.table(integrated.seurat4@meta.data, file='second_Larry_integrative_res1_labelled_meta.xls', sep='\t', quote=F)

## pdf('ALL_LARRY_batchcorrected_UMAP.pdf', width=16, height=5)
pdf('ALL_LARRY_batchcorrected_UMAP_res1.pdf', width=16, height=5)
p5d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
## p5d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F) + theme(legend.position='right')
## p5d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='integrated_snn_res.1', label=F) + theme(legend.position='right')
## p5d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F) + theme(legend.position='right')
## print(p4d / p5d)
print(p5d)
dev.off()

DefaultAssay(integrated.seurat4) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
integrated.seurat4 <- CellCycleScoring(integrated.seurat4, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(integrated.seurat4@meta.data)

pdf('ALL_LARRY_batchcorrected_UMAP_res1.pdf', width=16, height=12)
p5d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
p6d=FeaturePlot(integrated.seurat4, "S.Score", reduction='umap', split.by='orig.ident') + theme(legend.position='right')
p7d=FeaturePlot(integrated.seurat4, "G2M.Score", reduction='umap', split.by='orig.ident') + theme(legend.position='right')
print(p5d / p6d / p7d)
dev.off()

levels(Idents(integrated.seurat4)) = c("Ground", "EMT", "Interferon", "Muscle", "Prolif", "EMT")

integrated.seurat4$celltype.stim <- paste(Idents(integrated.seurat4), integrated.seurat4$orig.ident, sep = "_")

integrated.seurat4$celltype <- Idents(integrated.seurat4)
Idents(integrated.seurat4) <- "celltype.stim"

results.list <- list()
for (contrast in list(c("Regular_2", "Regular_1"),
                      c("Diff_3", "Regular_1"),
                      c("Diff_4", "Diff_3"),
                      c("Diff_4", "Regular_1"),
                      c("Diff_4", "Regular_2"),
                      c("Regular_2", "Diff_3"))) {
    de.contrast <- FindMarkers(integrated.seurat4,
                               ident.1 = paste0("EMT", "_", contrast[1]),
                               ident.2 = paste0("EMT", "_", contrast[2]),
                               verbose = FALSE)
    de.contrast$diff.pct <- de.contrast$pct.1 - de.contrast$pct.2
    print(summary(de.contrast$diff.pct))
    de.contrast.up = subset(de.contrast, diff.pct >= 0.1 & p_val_adj <= 0.01 & (pct.1 >= 0.1 | pct.2 >= 0.1))
    de.contrast.dn = subset(de.contrast, diff.pct <= -0.1 & p_val_adj <= 0.01 & (pct.1 >= 0.1 | pct.2 >= 0.1))
    results.list[[paste(contrast, collapse='_up_vs_', sep='_')]] = de.contrast.up
    results.list[[paste(contrast, collapse='_down_vs_', sep='_')]] = de.contrast.dn
}

for (rl in names(results.list)) {
    write.table(results.list[[rl]], file=paste0(rl, '_EMT_markers.xls'), sep='\t', quote=F)
}

## markers.list = list()
## for (l in levels(integrated.seurat4@meta.data$integrated_snn_res.1)) {
##     markers <- FindConservedMarkers(integrated.seurat4, ident.1 = l, grouping.var = "orig.ident", verbose = FALSE)
##     markers$cluster = l
##     markers.list[[l]] = markers
## }

## markers.res = list()
## for (l in names(markers.list)) {
##     markers.res[[l]] = markers.list[[l]][, c("max_pval", "minimump_p_val", "cluster")]
##     markers.res[[l]]$genes = rownames(markers.res[[l]])
## }
## markers.df = do.call(rbind, markers.res)
## write.table(markers.df, file="integrative_clusters_markers.txt", sep='\t', quote=F)

clonal = read.csv('RD_larry2_clone_mat_condition_specific2.csv', header=F)
cellbarcode = read.csv('RD_cellbarcode.txt', stringsAsFactors=F, header=T)
larrybarcode = scan("RD_larry2barcode_list_condition_specific2.txt", what='')

colnames(clonal) = larrybarcode
rownames(clonal) = paste0(cellbarcode[,3], cellbarcode[,2])

rna_cellbarcode = paste0(integrated.seurat4@meta.data$orig.ident, gsub('_\\d+', '', rownames(integrated.seurat4@meta.data)))

rna_clonal_match = clonal[match(rna_cellbarcode, rownames(clonal)), ]
sum(rownames(rna_clonal_match) == rna_cellbarcode)

DefaultAssay(integrated.seurat4) <- "RNA"

metadata = integrated.seurat4@meta.data[rowSums(rna_clonal_match) >= 1, ]

rna_clonal_match = rna_clonal_match[rowSums(rna_clonal_match)>=1, ]

rna_clonal_match_ge2 = rna_clonal_match[, colSums(rna_clonal_match) > 1] ## 3440 clones across two cells

metadata_ge2 = metadata[rowSums(rna_clonal_match_ge2) >= 1, ]

rna_clonal_match_ge2 = rna_clonal_match_ge2[rowSums(rna_clonal_match_ge2) >= 1, ]

dim(rna_clonal_match_ge2)

dim(subset(metadata_ge2, orig.ident=='Regular_1'))
sum(colSums(rna_clonal_match_ge2[metadata_ge2$orig.ident=='Regular_1', ])>=1)
sum(rowSums(rna_clonal_match_ge2[metadata_ge2$orig.ident=='Regular_1', ])>=1)

dim(subset(metadata_ge2, orig.ident=='Regular_2'))
sum(colSums(rna_clonal_match_ge2[metadata_ge2$orig.ident=='Regular_2', ])>=1)
sum(rowSums(rna_clonal_match_ge2[metadata_ge2$orig.ident=='Regular_2', ])>=1)

dim(subset(metadata_ge2, orig.ident=='Diff_3'))
sum(colSums(rna_clonal_match_ge2[metadata_ge2$orig.ident=='Diff_3', ])>=1)
sum(rowSums(rna_clonal_match_ge2[metadata_ge2$orig.ident=='Diff_3', ])>=1)

dim(subset(metadata_ge2, orig.ident=='Diff_4'))
sum(colSums(rna_clonal_match_ge2[metadata_ge2$orig.ident=='Diff_4', ])>=1)
sum(rowSums(rna_clonal_match_ge2[metadata_ge2$orig.ident=='Diff_4', ])>=1)

## Analysis 1: within the same tumor sample
for (cond in c("Regular_1", "Regular_2", "Diff_3", "Diff_4")) {
    clone = rna_clonal_match[metadata$orig.ident == cond, ]
    clone = clone[, colSums(clone) > 1]
    ## print(sum(rowSums(clone)>=1))
    state = metadata$seurat_clusters[metadata$orig.ident == cond]
    print(length(state))
    print(dim(clone))
    ## print(table(state))
    ## print(dim(clone))
    ## print(dim(state))
    clone_list = list()
    for (cl in colnames(clone)) {
        state_clone = sort(state[clone[, cl] == 1])
        ## state_clone[state_clone=='Interferon'] = 'EMT'
        ## state_clone = state_clone[state_clone!='Interferon']
        state_clone[grepl('EMT', state_clone)] = 'EMT'
        ## if (length(state_clone)==0)
        ##     next
        clone_list[[cl]] = sort(state_clone)
    }
    barcode_1 = as.data.frame(matrix(nrow=length(clone_list), ncol=max(unlist(lapply(clone_list, length)))))
    rownames(barcode_1) = colnames(clone)
    colnames(barcode_1) = paste0("cell_", seq(1, max(unlist(lapply(clone_list,length)))))
    for (cl in colnames(clone)) {
        current_states = clone_list[[cl]]
        print(current_states)
        barcode_1[cl, 1:length(current_states)] = current_states
    }
    ## sum(apply(barcode_1, 1, function(x) {
    ##     length(x[!is.na(x)])<=2
    ## }))
    ## data.frame(sort(table(apply(barcode_1[apply(barcode_1, 1, function(x) {
    ##     length(x[!is.na(x)])<=2
    ## }), ], 1, function(x) {
    ##     paste(na.omit(x), collapse='+')
    ## })), decreasing = T))
    write.table(barcode_1, file=paste0('barcode_within_tumor_', cond, '.txt'), sep='\t', quote=F)
}


## trace_clones = function(clone, state) {
##     trace_num = 0
##     cell_num = 0
##     for (cl in colnames(clone)) {
##         print(sum(clone[, cl] == 1))
##         cell_num <- cell_num + sum(clone[, cl] == 1)
##         trace_num <- trace_num + 1
##     }
##     cat(trace_num, '\t', cell_num)
##     connection = matrix(0, nrow=length(unique(state)), ncol=length(unique(state)))
##     rownames(connection) = unique(state)
##     colnames(connection) = unique(state)
##     for (cl in colnames(clone)) {
##         print('----')
##         print(cl)
##         if (sum(clone[, cl] == 1) >= 2) {
##             state_clone = state[clone[, cl] == 1]
##             print(state_clone)
##             for (row in seq_along(state_clone)[1:length(state_clone)]) {
##                 for (col in seq_along(state_clone)[1:(length(state_clone))]) {
##                     connection[state_clone[row], state_clone[col]] = connection[state_clone[row], state_clone[col]] + 1
##                 }
##             }
##         }
##     }
##     connection
## }


## library(corrplot)
## col3 <- colorRampPalette(c("red", "white", "blue")) 

## pdf("diff_regular_ratio.pdf")
## for (cond in c("Regular_1", "Regular_2", "Diff_3", "Diff_4")) {
##     regular1 = rna_clonal_match_ge2[metadata_ge2$orig.ident==cond, ]
##     meta1 = metadata_ge2[metadata_ge2$orig.ident==cond, ]
##     conn1 = trace_clones(regular1, meta1$seurat_clusters)
##     corrplot(conn1, col=col3(20), method = "number", cl.lim=c(0, max(conn1)), order="AOE", is.corr = FALSE, type = "upper")
##     ## corrplot(conn1/sum(conn1), col=col3(20), method = "number", cl.lim=c(0, 0.6), is.corr = FALSE)
## }
## dev.off()

## pdf("diff_regular_ratio2.pdf")
## for (cond in c("Regular_1", "Regular_2", "Diff_3", "Diff_4")) {
##     regular1 = rna_clonal_match_ge2[metadata_ge2$orig.ident==cond, ]
##     meta1 = metadata_ge2[metadata_ge2$orig.ident==cond, ]
##     conn1 = trace_clones(regular1, meta1$seurat_clusters)
##     corrplot(conn1/sum(conn1), col=col3(20), method = "number", cl.lim=c(0, 0.6), is.corr = FALSE, order='AOE', type='upper')
## }
## dev.off()


## trace_num = 0
## cell_num = 0
## ## dim(rna_clonal_match_ge2)
## ## dim(subset(metadata_ge2, orig.ident=='Regular_1'))

## state1 = metadata_ge2[metadata_ge2$orig.ident=='Regular_1', ]$seurat_clusters
## state2 = metadata_ge2[metadata_ge2$orig.ident=='Regular_2', ]$seurat_clusters
## connection = matrix(0, nrow=length(unique(state1)), ncol=length(unique(state2)))
## rownames(connection) = unique(state1)
## colnames(connection) = unique(state2)
## for (cl in colnames(rna_clonal_match_ge2)) {
##     meta_time = metadata_ge2[rna_clonal_match_ge2[, cl] == 1, ]$orig.ident
##     meta_state = metadata_ge2[rna_clonal_match_ge2[, cl] == 1, ]$seurat_clusters
##     if(sum(unique(meta_time)%in%c("Regular_1", "Regular_2"))>=2) {
##         cell_num <- cell_num + sum(rna_clonal_match_ge2[, cl] == 1)
##         trace_num <- trace_num + 1
##         start = meta_state[meta_time == "Regular_1"]
##         end = meta_state[meta_time == "Regular_2"]
##         print('--------')
##         cat(start, '\t\t\t', end, '\n')
##         for (row in seq_along(start)[1:length(start)]) {
##             for (col in seq_along(end)[1:(length(end))]) {
##                 connection[start[row], end[col]] = connection[start[row], end[col]] + 1
##             }
##         }
##     }
## }

## cat(trace_num, '\t', cell_num)

## ## pdf("diff_regular_ratio_time.pdf")
## ## corrplot(connection, col=col3(20), method = "number", cl.lim=c(0, 200), is.corr = FALSE, order='AOE', type='upper')
## ## dev.off()


## state1 = metadata_ge2[metadata_ge2$orig.ident=='Regular_1', ]$seurat_clusters
## state2 = metadata_ge2[metadata_ge2$orig.ident=='Diff_3', ]$seurat_clusters
## connection = matrix(0, nrow=length(unique(state1)), ncol=length(unique(state2)))
## rownames(connection) = unique(state1)
## colnames(connection) = unique(state2)
## for (cl in colnames(rna_clonal_match_ge2)) {
##     meta_time = metadata_ge2[rna_clonal_match_ge2[, cl] == 1, ]$orig.ident
##     meta_state = metadata_ge2[rna_clonal_match_ge2[, cl] == 1, ]$seurat_clusters
##     if(sum(unique(meta_time)%in%c("Regular_1", "Diff_3"))>=2) {
##         cell_num <- cell_num + sum(rna_clonal_match_ge2[, cl] == 1)
##         trace_num <- trace_num + 1
##         start = meta_state[meta_time == "Regular_1"]
##         end = meta_state[meta_time == "Diff_3"]
##         print('--------')
##         cat(start, '\t\t\t', end, '\n')
##         for (row in seq_along(start)[1:length(start)]) {
##             for (col in seq_along(end)[1:(length(end))]) {
##                 connection[start[row], end[col]] = connection[start[row], end[col]] + 1
##             }
##         }
##     }
## }

## Analysis 2: fate change across samples
## include "parental cells"
## 1. only a single cell has a barcode
## 2. two (or more) parental cells share the same phenotype
## 3. time series contains barcoded cells
## Regular_1 vs Regular_2
comparisons = list(c('Regular_1', 'Regular_2'),
                   c('Regular_1', 'Diff_3'),
                   c("Diff_3", "Diff_4"))

for (comp in comparisons) {
    print(comp)
    results = list()
    ## for each barcoded clone
    for (cl in colnames(rna_clonal_match_ge2)) {
        meta_time = metadata_ge2[rna_clonal_match_ge2[, cl] == 1, ]$orig.ident
        if(sum(c(comp) %in% unique(meta_time))==2) {
            meta_state = metadata_ge2[rna_clonal_match_ge2[, cl] == 1, ]$seurat_clusters
            start_cond = meta_state[meta_time == comp[1]]
            end_cond = meta_state[meta_time == comp[2]]
            if (length(unique(start_cond)) == 1) {
                results[[cl]] = list()
                results[[cl]][['parent']] = start_cond
                results[[cl]][['progeny']] = end_cond
            }
        }
    }
    max_parent = max(unlist(lapply(results, function(x) {
        length(x[['parent']])
    })))
    max_progeny = max(unlist(lapply(results, function(x) {
        length(x[['progeny']])
    })))
    barcode_cross_sample = as.data.frame(matrix(nrow=length(results), ncol=max_parent+max_progeny))
    rownames(barcode_cross_sample) = names(results)
    colnames(barcode_cross_sample) = c(paste0("parent_", seq(1, max_parent)),
                                       paste0("progeny_", seq(1, max_progeny)))
    for (cl in names(results)) {
        print(cl)
        barcode_cross_sample[cl, 1:length(results[[cl]][['parent']])] = results[[cl]][['parent']]
        barcode_cross_sample[cl, (max_parent+1):(max_parent+1+length(results[[cl]][['progeny']])-1)] = results[[cl]][['progeny']]
    }
    write.table(barcode_cross_sample, file=paste0(paste(comp, collapse='_vs_barcoded_clones'), '.xls'), quote=F, sep='\t')
}
