library(Seurat)
library(reticulate)
library(tidyverse)
library(patchwork)

cell1 = list(readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190815_seurat-object_MAST85-1cell.rds'),
             readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST85.rds'))

cell2 = list(readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190815_seurat-object_RH74-10cells.rds'),
             readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_RH74.rds'))

integration <- function(x, y, label, method='seurat') {
    seurat.pseudo.list = list(x, y)
    for (i in 1:length(seurat.pseudo.list)) {
        print(1)
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
integrated.seurat = integration(cell1[[1]], cell1[[2]], c('MAST85', 'MAST85-1cell'))
integrated.seurat2 = integration(cell2[[1]], cell2[[2]], c('RH74-10cell', 'RH74'))

colortable = read.table('color_table.xls', sep='\t', header=T, stringsAsFactors=F)
new.ident = colortable[, c('MAST85')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
MAST85 = cell1[[2]]
Idents(MAST85) = MAST85@meta.data$seurat_clusters = MAST85$RNA_snn_res.0.8
MAST85 = RenameIdents(MAST85, new.ident)
idents1 = Idents(MAST85)
names(idents1) = paste0(names(idents1), '_2')
new.ident = colortable[, c('MAST85.1cell')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
MAST85.1 = cell1[[1]]
Idents(MAST85.1) = MAST85.1@meta.data$seurat_clusters = MAST85.1$RNA_snn_res.0.8
MAST85.1 = RenameIdents(MAST85.1, new.ident)
idents2 = Idents(MAST85.1)
names(idents2) = paste0(names(idents2), '_1')
all.ident = Idents(integrated.seurat)
idents = unlist(list(idents1, idents2))
Idents(integrated.seurat) = idents[match(names(all.ident), names(idents))]

colortable = read.table('color_table.xls', sep='\t', header=T, stringsAsFactors=F)
new.ident = colortable[, c('RH74')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
RH74 = cell2[[2]]
Idents(RH74) = RH74@meta.data$seurat_clusters = RH74$RNA_snn_res.0.8
RH74 = RenameIdents(RH74, new.ident)
idents1 = Idents(RH74)
names(idents1) = paste0(names(idents1), '_2')
new.ident = colortable[, c('RH74.10cells')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
RH74.1 = cell2[[1]]
Idents(RH74.1) = RH74.1@meta.data$seurat_clusters = RH74.1$RNA_snn_res.0.8
RH74.1 = RenameIdents(RH74.1, new.ident)
idents2 = Idents(RH74.1)
names(idents2) = paste0(names(idents2), '_1')
all.ident = Idents(integrated.seurat2)
idents = unlist(list(idents1, idents2))
Idents(integrated.seurat2) = idents[match(names(all.ident), names(idents))]

metacolors = c(rgb(131, 81, 9, maxColorValue = 255),
               rgb(213, 139, 36, maxColorValue = 255),
               rgb(67, 20, 122, maxColorValue = 255),
               rgb(24, 90, 224, maxColorValue = 255),
               rgb(233, 233, 233, maxColorValue = 255),
               rgb(65, 128, 255, maxColorValue = 255),
               rgb(254, 43, 3, maxColorValue = 255),
               rgb(65, 169, 92, maxColorValue = 255),
               rgb(48, 112, 152, maxColorValue = 255),
               rgb(251, 113, 172, maxColorValue = 255))
metalabels = c("GROUND", "Hypoxia", "EMT",  "G1S", "UNASSIGNED", "G2M",  "MUSCLE", "INTERFERON", "PROLIF", "Histones")
names(metacolors) = metalabels

pdf('Integrated.pdf', width=12, height=12)
p1d=DimPlot(integrated.seurat, reduction='umap', split.by='orig.ident', label=F, cols=metacolors)+ theme(legend.position='bottom')
p2d=DimPlot(integrated.seurat2, reduction='umap', split.by='orig.ident', label=F, cols=metacolors) + theme(legend.position='none')
pd = CombinePlots(plots=list(p1d, p2d), ncol=1)
print(pd)
dev.off()

library(scales)
library(ggpubr)
## https://www.r-graph-gallery.com/128-ring-or-donut-plot.html
p1 = reshape2::melt(
              table(Idents(integrated.seurat), integrated.seurat@meta.data[,1])
          ) %>% filter(Var2=='MAST85' & value>0) %>% mutate(fraction=value/sum(value)) %>% mutate(ymax=cumsum(fraction)) %>% mutate(ymin=c(0, head(ymax, n=-1))) %>% mutate(labelPosition=(ymax+ymin)/2, label=paste0(Var1, ':', value), Var1=as.factor(as.vector(Var1))) %>% ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) + geom_rect() + scale_fill_manual(values=metacolors) + xlim(c(2, 4)) + coord_polar(theta="y")  +  annotate(geom="text", x=2, y=0.1, label="MAST85", color="black") + theme_void() + theme(legend.position = "none") # + geom_label(x=3.5, aes(y=labelPosition, label=label), size=2)
p2 = reshape2::melt(
              table(Idents(integrated.seurat), integrated.seurat@meta.data[,1])
          ) %>% filter(Var2=='MAST85-1cell' & value>0) %>% mutate(fraction=value/sum(value)) %>% mutate(ymax=cumsum(fraction)) %>% mutate(ymin=c(0, head(ymax, n=-1))) %>% mutate(labelPosition=(ymax+ymin)/2, label=paste0(Var1, ':', value), Var1=as.factor(as.vector(Var1))) %>% ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) + geom_rect() + scale_fill_manual(values=metacolors) + xlim(c(2, 4)) + coord_polar(theta="y") +  annotate(geom="text", x=2, y=0.1, label="MAST85-1cell", color="black") + theme_void() + theme(legend.position = "none") # geom_label(x=3.5, aes(y=labelPosition, label=label), size=2) + 
p3 = reshape2::melt(
              table(Idents(integrated.seurat2), integrated.seurat2@meta.data[,1])
          ) %>% filter(Var2=='RH74' & value>0) %>% mutate(fraction=value/sum(value)) %>% mutate(ymax=cumsum(fraction)) %>% mutate(ymin=c(0, head(ymax, n=-1))) %>% mutate(labelPosition=(ymax+ymin)/2, label=paste0(Var1, ':', value), Var1=as.factor(as.vector(Var1))) %>% ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) + geom_rect() + scale_fill_manual(values=metacolors) + xlim(c(2, 4)) + coord_polar(theta="y")  +  annotate(geom="text", x=2, y=0.1, label="RH74", color="black") + theme_void() + theme(legend.position = "none") # + geom_label(x=3.5, aes(y=labelPosition, label=label), size=2)
p4 = reshape2::melt(
              table(Idents(integrated.seurat2), integrated.seurat2@meta.data[,1])
          ) %>% filter(Var2=='RH74-10cells' & value>0) %>% mutate(fraction=value/sum(value)) %>% mutate(ymax=cumsum(fraction)) %>% mutate(ymin=c(0, head(ymax, n=-1))) %>% mutate(labelPosition=(ymax+ymin)/2, label=paste0(Var1, ':', value), Var1=as.factor(as.vector(Var1))) %>% ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Var1)) + geom_rect() + scale_fill_manual(values=metacolors) + xlim(c(2, 4)) + coord_polar(theta="y") +  annotate(geom="text", x=2, y=0.1, label="RH74-10cells", color="black") + theme_void() + theme(legend.position = "none") # + guides(fill=guide_legend(nrow=2,byrow=TRUE))  # geom_label(x=3.5, aes(y=labelPosition, label=label), size=2) + 
(pd) + ((p1 + p2 + p3 + p4) + plot_layout(ncol=1)) + plot_layout(width=c(4, 1))
ggsave('Fig5_DE.pdf', width=16, height=16)
### p1 + gridExtra::tableGrob(mtcars[1:10, c('mpg', 'disp')])

test = merge(test1, test2, by='Var1')
cor.test(test$fraction.x, test$fraction.y, method='pearson')
