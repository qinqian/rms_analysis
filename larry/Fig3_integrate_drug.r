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
                rgb(241, 250, 100, maxColorValue = 255),
                rgb(0, 255, 253, maxColorValue = 255))
metalabels <- c("Ground", "Hypoxia", "EMT",
                "G1S", "UNASSIGNED",
                "G2M",  "Muscle", "Interferon", "Prolif",
                "Histone", "Apoptosis", 'UPR', "Neural")
names(metacolors) <- metalabels

library(Seurat)
library(reticulate)
library(patchwork)
library(ggplot2)

res <- readRDS('results/seurat_sara/Differentiated_seurat-object.rds')
sen <- readRDS('results/seurat_sara/Regular_seurat-object.rds')

## rd <- readRDS('../../projects/01_sc_rms/results/seurat_sara/RD_seurat-object.rds')

res$seurat_clusters = res$RNA_snn_res.0.8
sen$seurat_clusters = sen$RNA_snn_res.0.8
## rd$seurat_clusters = rd$RNA_snn_res.0.8

## cell4 = list(res, sen, rd)
cell4 = list(res, sen)

## Merge of all RD cells assuming no batch effect
## cells = Reduce(merge, cell4)
## cells = NormalizeData(cells)
## cells = FindVariableFeatures(cells, selection.method='vst', nfeatures=2000)
## cells <- ScaleData(cells, vars.to.regress = "percent.mito")

## seurat.obj <- RunPCA(object=cells)
## seurat.obj <- RunUMAP(object=seurat.obj, dims=1:20)
## seurat.obj <- RunTSNE(object=seurat.obj, dims=1:20)
## seurat.obj <- FindNeighbors(object=seurat.obj) ## use dims later...

library(patchwork)

## pdf('ALL_LARRY_UMAP.pdf', width=8, height=10)
## p4d=DimPlot(seurat.obj, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=T) + theme(legend.position='none')
## p5d=DimPlot(seurat.obj, reduction='umap', group.by='orig.ident', label=T) + theme(legend.position='none')
## print(p4d/p5d)
## dev.off()

## Integration by CCA analysis
integration <- function(x, y, label, method='seurat') {
    Idents(x) = x$seurat_clusters2 = x$seurat_clusters
    states = unlist(as.vector(annotation[label[1], ]))
    states = states[states!='']
    states = states[!is.na(states)]
    levels(x$seurat_clusters) = states
    print('-------')
    ## levels(Idents(x)) = states
    Idents(y) = y$seurat_clusters2 = y$seurat_clusters
    states = unlist(as.vector(annotation[label[2], ]))
    states = na.omit(states[states!=''])
    print(length(states))
    print(levels(y$seurat_clusters))
    levels(y$seurat_clusters) = states
    ## levels(Idents(y)) = states
    ## Idents(z) = z$seurat_clusters2 = z$seurat_clusters
    ## states = unlist(as.vector(annotation[label[3], ]))
    ## states = na.omit(states[states!=''])
    ## levels(z$seurat_clusters) = states
    ## levels(Idents(z)) = states
    ## seurat.pseudo.list = list(x, y, z)
    seurat.pseudo.list = list(x, y)
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

set.seed(100)

annotation = read.delim('LARRY_cellstate.txt', sep='\t', row.names=1, header=T,
                        check.names=F, stringsAsFactors=F)
integrated.seurat4 = integration(cell4[[1]], cell4[[2]], c('Differentiated', 'Regular'))


pdf('ALL_LARRY_batchcorrected_UMAP.pdf', width=11, height=6)
## p4d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='seurat_clusters2', label=T) + theme(legend.position='right')
p5d=DimPlot(integrated.seurat4, reduction='umap', split.by='orig.ident', group.by='seurat_clusters', label=F, cols=metacolors) + theme(legend.position='right')
## print(p4d / p5d)
print(p5d)
dev.off()

## clonal = read.csv('RD_larry_clone_mat_condition_specific.csv', header=F)
## cellbarcode = read.csv('RD_cellbarcode_list_condition_specific.txt', stringsAsFactors=F)
## larrybarcode = scan("RD_larrybarcode_list_condition_specific.txt", what='')

clonal = read.csv('RD_larry_clone_mat_condition_specific2.csv', header=F)
cellbarcode = read.csv('RD_cellbarcode_list_condition_specific.txt', stringsAsFactors=F)
larrybarcode = scan("RD_larrybarcode_list_condition_specific2.txt", what='')

colnames(clonal) = larrybarcode
rownames(clonal) = paste0(cellbarcode[,3], cellbarcode[,2])

rna_cellbarcode = paste0(integrated.seurat4@meta.data$orig.ident, gsub('_\\d+', '', rownames(integrated.seurat4@meta.data)))
rna_clonal_match = clonal[match(rna_cellbarcode, rownames(clonal)),]
sum(rownames(rna_clonal_match) == rna_cellbarcode)

DefaultAssay(integrated.seurat4) <- "RNA"

head(integrated.seurat4@meta.data)
head(integrated.seurat4@meta.data$seurat_clusters)

rna_clonal_match_ge2 = rna_clonal_match[, colSums(rna_clonal_match) > 1]

sum(rowSums(rna_clonal_match_ge2) >= 1)
sum(colSums(rna_clonal_match_ge2) >= 1)

times = ifelse(integrated.seurat4$orig.ident=='Differentiated', 1, 0)
state = integrated.seurat4$seurat_clusters

write.table(integrated.seurat4@meta.data, file='RD_Larry_metadata.xls', quote=F, sep='\t')
write.table(integrated.seurat4@reductions$umap@cell.embeddings, file='RD_Larry_umap.xls', quote=F, sep='\t')

sum(colSums(rna_clonal_match_ge2[times==1, ])>=1) # clone number for differentiated
sum(rna_clonal_match_ge2[times==1, ]) # cell number

sum(colSums(rna_clonal_match_ge2[times==0, ])>=1) # clone number for regular
sum(rna_clonal_match_ge2[times==0, ]) # cell number

sum(colSums(rna_clonal_match_ge2[times == 1, ]) > 1)
sum(colSums(rna_clonal_match_ge2[times == 0, ]) > 1)

## Latest statistics table from Dave
clone = rna_clonal_match[times == 1, colSums(rna_clonal_match[times == 1, ]) > 1]
state = integrated.seurat4$seurat_clusters[times == 1]

clone_list = list()
for (cl in colnames(clone)) {
    print(cl)
    state_clone = sort(state[clone[, cl] == 1])
    clone_list[[cl]] = state_clone
}
barcode_1 = as.data.frame(matrix(nrow=length(clone_list), ncol=max(unlist(lapply(clone_list,length)))))
rownames(barcode_1) = colnames(clone)
for (cl in colnames(clone)) {
    state_clone = sort(state[clone[, cl] == 1])
    barcode_1[cl, 1:length(state_clone)] = state_clone
}
sum(apply(barcode_1, 1, function(x) {
    length(x[!is.na(x)])<=2
}))
data.frame(sort(table(apply(barcode_1[apply(barcode_1, 1, function(x) {
    length(x[!is.na(x)])<=2
}), ], 1, function(x) {
    paste(na.omit(x), collapse='+')
})), decreasing = T))
write.table(barcode_1, file='barcode_1_differentiated.txt', sep='\t', quote=F)


clone = rna_clonal_match[times == 0, colSums(rna_clonal_match[times == 0, ]) > 1]
state = integrated.seurat4$seurat_clusters[times == 0]

clone_list = list()
for (cl in colnames(clone)) {
    print(cl)
    state_clone = sort(state[clone[, cl] == 1])
    clone_list[[cl]] = state_clone
}
barcode_2 = as.data.frame(matrix(nrow=length(clone_list), ncol=max(unlist(lapply(clone_list,length)))))

rownames(barcode_2) = colnames(clone)
for (cl in colnames(clone)) {
    state_clone = sort(state[clone[, cl] == 1])
    barcode_2[cl, 1:length(state_clone)] = state_clone
}

sum(apply(barcode_2, 1, function(x) {
    length(x[!is.na(x)])<=2
}))

data.frame(sort(table(apply(barcode_2[apply(barcode_2, 1, function(x) {
    length(x[!is.na(x)])<=2
}), ], 1, function(x) {
    paste(na.omit(x), collapse='+')
})), decreasing = T))

write.table(barcode_2, file='barcode_2_regular.txt', sep='\t', quote=F)


trace_clones = function(clone, state) {
    trace_num = 0
    cell_num = 0
    for (cl in colnames(clone)) {
        print(sum(clone[, cl] == 1))
        cell_num <- cell_num + sum(clone[, cl] == 1)
        trace_num <- trace_num + 1
    }
    cat(trace_num, '\t', cell_num)
    ## x=character(trace_num)
    ## y=character(trace_num)
    ## i = 0
    connection = matrix(0, nrow=length(unique(state)), ncol=length(unique(state)))
    rownames(connection) = unique(state)
    colnames(connection) = unique(state)
    allconn = 0
    for (cl in colnames(clone)) {
        state_clone = state[clone[, cl] == 1]
        ## for (row in seq_along(state_clone)[1:(length(state_clone)-1)]) {
        ##     for (col in seq_along(state_clone)[row:(length(state_clone))]) {
        for (row in seq_along(state_clone)[1:length(state_clone)]) {
            for (col in seq_along(state_clone)[1:(length(state_clone))]) {
                connection[state_clone[row], state_clone[col]] = connection[state_clone[row], state_clone[col]] + 1
            }
        }
        ## x[i] = paste(start, collapse=',')
        ## y[i] = paste(end, collapse=',')
    }
    ## df = data.frame(start=x,end=y)
    ## print(table(df$start, df$end))
    ## df = as.data.frame(table(df$start, df$end))
    ## df
    connection
}

diff_clone = rna_clonal_match[times == 1, colSums(rna_clonal_match[times == 1, ]) > 1]
diff_state = state[times == 1]
df1 = trace_clones(diff_clone, diff_state)

reg_clone = rna_clonal_match[times == 0, colSums(rna_clonal_match[times == 0, ]) > 1]
reg_state = state[times == 0]
df2 = trace_clones(reg_clone, reg_state)

library(corrplot)

col3 <- colorRampPalette(c("red", "white", "blue")) 

## pdf("diff_regular.pdf")
pdf("diff_regular_ratio.pdf")
## corrplot(df1/sum(df1), col=col3(20), method = "number", cl.lim=c(0, max(df1)), order="AOE", is.corr = FALSE, type = "upper")
## corrplot(df2/sum(df2), col=col3(20), method = "number", cl.lim=c(0, max(df2)), order="AOE", is.corr = FALSE, type = "upper")
corrplot(df1/sum(df1), col=col3(20), method = "number", cl.lim=c(0, 0.6), is.corr = FALSE)
corrplot(df2/sum(df2), col=col3(20), method = "number", cl.lim=c(0, 0.3), is.corr = FALSE)
## corrplot(df1/sum(df1)-df2/sum(df2), col=col3(20), method = "number", cl.lim=c(-1, 1), order="AOE", is.corr = FALSE, type = "upper")
dev.off()

## df2 = reshape2::melt(as.data.frame(df2))

library(ggplot2)
library(scales) # for muted function
p1 = ggplot(df1, aes(Var2, Var1)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = Freq)) + # background colours are mapped according to the value column
  geom_text(aes(label = round(Freq, 2))) + # write the values
  scale_fill_gradient2(low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"), 
                       midpoint = 0) + # determine the colour
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        ## panel.background=element_rect("white"), # background=white
        axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 12,face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold")) + 
  ggtitle("Larry barcode transition plot") + 
  theme(legend.title=element_text(face="bold", size=14)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs("Transition cells between Regular and Diff")
ggsave('test.pdf')


## no common ancesters
## clones across the two samples
## however, the two samples do not have common ancesters
## this lead to low number of shared barcodes across cells
trace_num = 0
cell_num = 0
for (cl in colnames(rna_clonal_match_ge2)) {
    if(length(table(times[rna_clonal_match_ge2[, cl] == 1]))==2) {
        cell_num <- cell_num + sum(rna_clonal_match_ge2[, cl] == 1)
        trace_num <- trace_num + 1
        print('-----')
        start = state[rna_clonal_match_ge2[, cl] == 1 & times==0]
        end = state[rna_clonal_match_ge2[, cl] == 1 & times==1]
        print(start)
        print(end)
        print(table(start))
        print(table(end))
        start = unique(start)
        end = unique(end)
    }
}
cat(trace_num, '\t', cell_num)
x=character(trace_num)
y=character(trace_num)
i = 0
for (cl in colnames(rna_clonal_match_ge2)) {
    if(length(table(times[rna_clonal_match_ge2[, cl] == 1]))==2) {
        i <- i+1
        start = unique(state[rna_clonal_match_ge2[, cl] == 1 & times==0])
        end = unique(state[rna_clonal_match_ge2[, cl] == 1 & times==1])
        x[i] = start
        y[i] = paste(end, collapse=',')
    }
}
df = data.frame(start=x,end=y)
df = as.data.frame(table(df$start, df$end))

library(ggplot2)
library(scales) # for muted function
ggplot(df, aes(Var2, Var1)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = Freq)) + # background colours are mapped according to the value column
  geom_text(aes(label = round(Freq, 2))) + # write the values
  scale_fill_gradient2(low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"), 
                       midpoint = 0) + # determine the colour
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        ## panel.background=element_rect("white"), # background=white
        axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 12,face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold")) + 
  ggtitle("Larry barcode transition plot") + 
  theme(legend.title=element_text(face="bold", size=14)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs("Transition cells between Regular and Diff")
ggsave('test.pdf')
