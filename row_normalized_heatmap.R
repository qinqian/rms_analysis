### https://stackoverflow.com/questions/13115345/need-help-installing-spreadsheetparseexcel
### fix conda excel reader

library(gdata)
library(Seurat)
library(gdata)
library(tidyverse)

gene.sign <- read.xls('Gene\ Modules\ for\ RMS.xlsx',
                      stringsAsFactors = F, perl = "/usr/bin/perl")
gene.sign = lapply(as.list(gene.sign), function(x) {x[x!=""]})
gene.sign.union = Reduce(union, gene.sign)

only.sign = read.xls("RMS Core signature gene list.xlsx", stringsAsFactors=F, perl="/usr/bin/perl")
only.sign = only.sign[,c(1,3)]
only.sign = lapply(as.list(only.sign), function(x) {x[x!=""]})
only.sign.union = Reduce(union, only.sign)

MAST111 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds')
MAST139 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds')
MAST39 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST39.rds')
RH74 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_RH74.rds')
MAST95 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST95.rds')
MAST85 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST85.rds')
MSK82489 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MSK82489.rds')
MAST35 = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST35.rds')

colortable = read.table('color_table.xls', sep='\t', header=T, stringsAsFactors=F)
new.ident = colortable[, c('MAST85')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
Idents(MAST85) = MAST85@meta.data$seurat_clusters = MAST85$RNA_snn_res.0.8
MAST85 = RenameIdents(MAST85, new.ident)

colortable = read.table('color_table.xls', sep='\t', header=T, stringsAsFactors=F)
new.ident = colortable[, c('MSK82489')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
Idents(MSK82489) = MSK82489@meta.data$seurat_clusters = MSK82489$RNA_snn_res.0.8
MSK82489 = RenameIdents(MSK82489, new.ident)

colortable = read.table('color_table.xls', sep='\t', header=T, stringsAsFactors=F)
new.ident = colortable[, c('MAST35')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
Idents(MAST35) = MAST35@meta.data$seurat_clusters = MAST35$RNA_snn_res.0.8
MAST35 = RenameIdents(MAST35, new.ident)

colortable = read.table('color_table.xls', sep='\t', header=T, stringsAsFactors=F)
new.ident = colortable[, c('MAST39')]
names(new.ident) = colortable[, 'cluster']
Idents(MAST39) = MAST39@meta.data$seurat_clusters = MAST39$RNA_snn_res.0.8
MAST39 = RenameIdents(MAST39, new.ident)

new.ident = colortable[, c('MAST95')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
Idents(MAST95) = MAST95@meta.data$seurat_clusters = MAST95$RNA_snn_res.0.8
MAST95 = RenameIdents(MAST95, new.ident)

new.ident = colortable[, c('RH74')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
Idents(RH74) = RH74@meta.data$seurat_clusters = RH74$RNA_snn_res.0.8
RH74 = RenameIdents(RH74, new.ident)

new.ident = colortable[, c('MAST111')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
Idents(MAST111) = MAST111@meta.data$seurat_clusters = MAST111$RNA_snn_res.0.8
MAST111 = RenameIdents(MAST111, new.ident)

new.ident = colortable[, c('MAST139')]
names(new.ident) = colortable[, 'cluster']
new.ident = new.ident[new.ident!='']
Idents(MAST139) = MAST139@meta.data$seurat_clusters = MAST139$RNA_snn_res.0.8
MAST139 = RenameIdents(MAST139, new.ident)

###############
### Fig 3B
###############

df = read.xls('for heatmap.xls', stringsAsFactors=F, perl='/usr/bin/perl')
rownames(df) = df$Gene
df <- df[,-1]

include = setdiff(only.sign.union, gene.sign.union)

targets = intersect(rownames(MAST111), include)
targets = intersect(targets, rownames(RH74))
targets = intersect(targets, rownames(MAST139))
targets = intersect(targets, rownames(MAST95))
targets = intersect(targets, rownames(MAST39))

df = df[include, ]

library(pheatmap)
library(gplots)
library(RColorBrewer)

heatmap(as.matrix(df), scale="row", Rowv = NA, Colv=NA, 
        col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias=0.9)(99),
        breaks=seq(-2, 2, length=100))

pheatmap(as.matrix(df), scale="row", cluster_rows=F, cluster_cols=F, fontsize_row=0.8,show_rownames=F, breaks=seq(-2, 2, length=99),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias=0.9)(100),
         width=4, height=8)

pheatmap(as.matrix(df), scale="row", cluster_rows=F, cluster_cols=F, fontsize_row=0.8,show_rownames=F, breaks=seq(-2, 2, length=99),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias=0.9)(100),
         filename="row_normalized_heatmap.pdf",
         width=4, height=8)

###############
### Fig 3C
###############

df = df[targets,]
folds = rowMeans(df[,1:6]) / rowMeans(df[,7:8])

erms.topsign = tail(names(sort(folds)), 10)
arms.topsign = head(names(sort(folds)), 10)

library(gridExtra)
pdf("Fig4C_ARMS.pdf", width=13, height=6)
for (i in arms.topsign) {
  p1 = FeaturePlot(MAST39, features=i)
  p2 = FeaturePlot(MAST95, features=i)
  grid.arrange(p1, p2,
               ncol=2,
               top="MAST39 vs MAST95 ARMS signature")
}
dev.off()

pdf("Fig4C_ERMS.pdf", width=13, height=6)
for (i in erms.topsign) {
  p1 = FeaturePlot(MAST39, features=i)
  p2 = FeaturePlot(MAST95, features=i)
  grid.arrange(p1, p2,
               ncol=2,
               top="MAST39 vs MAST95 ERMS signature")
}
dev.off()


library(patchwork)

plist = list()
pdf("Fig3B_PIPOX_FABP5.pdf", width=10, height=7)
for (i in c("PIPOX", "FABP5")) {
  p1 = FeaturePlot(MAST39, features=i, cols = c("lightgrey",  "red"))
  p2 = FeaturePlot(MAST95, features=i, cols = c("lightgrey",  "red"))
  plist[[i]] = CombinePlots(plots=list(p1, p2), ncol=2)
}
plist[[1]] + plist[[2]] + plot_layout(ncol=1)
dev.off()

selection1 = c("MAGEL2", "CCND1", "LTBP4", "HMGA2", "KAZALD1",
               "EMILIN1", "FABP5", "GPR124", "SH3BP5", "RGS10")
selection2 = c("TFAP2B", "HSPB2", "COX7A1", "PIPOX", "SPINT2",
               "CRMP1",  "NHLH1", "TAGLN3", "MYOG",  "EDN3")

library(gridExtra)
pdf("Fig4C_demo_10genes.pdf", width=12, height=9)
for (i in 1:length(selection1)) {
  p1 = FeaturePlot(MAST39, features=selection1[i])
  p2 = FeaturePlot(MAST95, features=selection1[i])
  p3 = FeaturePlot(MAST39, features=selection2[i])
  p4 = FeaturePlot(MAST95, features=selection2[i])
  pa = CombinePlots(plots=list(p1, p2), ncol=1)
  pb = CombinePlots(plots=list(p3, p4), ncol=1)
  grid.arrange(pa, pb,
               ncol=2,
               top="MAST39 vs MAST95 ARMS signature")
}
dev.off()

library(ggplot2)
library(cowplot)
library(patchwork)

DotPlot = function (object, assay = NULL, features, cols = c("lightgrey", 
    "red"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
    group.by = NULL, split.by = NULL, scale.by = "radius", scale.min = NA, 
    scale.max = NA, return_data=F) {
    assay <- DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size, 
        radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    if (sum(!(features %in% colnames(data.features)))>0) {
        missing.features = setdiff(features, colnames(data.features))
        ## assign to zeros
        missing.data.features = matrix(0, nrow=nrow(data.features), ncol=length(missing.features))
        colnames(missing.data.features) = missing.features
        print(dim(missing.data.features))
        print(dim(data.features))
        data.features = cbind(data.features, missing.data.features)
    }
    data.features$id <- if (is.null(x = group.by)) {
        Idents(object = object)
    }
    else {
        object[[group.by, drop = TRUE]]
    }
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident, 
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = Seurat:::PercentAbove, 
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot == 
                x, "avg.exp"]
            data.use <- MinMax(data = data.use, min = col.min, 
                max = col.max)
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, 
        levels = rev(x = features))
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = strsplit(x = as.character(x = data.plot$id), 
            split = "_"), FUN = "[[", FUN.VALUE = character(length = 1L), 
            2)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled", 
        no = "colors")
    ## color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp", 
    ##     no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    plot <- ggplot(data = data.plot, mapping = aes_string(y = "features.plot", 
        x = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                     color = color.by)) + scale.func(range = c(0, dot.scale),
                                                                                     breaks = seq(0, 100, 20),
        limits=c(0, 100)) + theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
        labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
            yes = "Identity", no = "Split Identity")) + theme_cowplot() 
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    }
    else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    }
    else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2], breaks=seq(col.min, col.max, 0.5), limits=c(0, col.max))
    }

    if (object@project.name == 'MAST39') {
        plot <- plot + scale_x_discrete(limits=c("G1S", "G2M", "PROLIF", "Histones", "Hypoxia", "INTERFERON", "EMT", "GROUND"))
    } else {
        plot <- plot + scale_x_discrete(limits=c("G1S", "G2M", "Hypoxia", "MUSCLE", "GROUND", "UNASSIGNED"))
    }

    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    if (return_data) {
        return(data.plot)
    }
    return(plot)
}

## arms.topsign = c("TFAP2B", "HSPB2", "COX7A1", "TFF3", "CRMP1", "MYL4", "DCX",  "EDN3", "ENO3", "TMEM47")
## erms.topsign = c("ITM2A", "FABP5", "EMILIN1", "HMGA2", "LTBP4", "EMP3", "HEY1", "EFEMP2", "TSTA3", "EIF4EBP1")
erms.topsign = c("ITM2A", "FABP5", "EMILIN1", "RRBP1", "CCND1", "EMP3", "EIF4EBP1", "TSTA3", "ADA", "HMGA2")
arms.topsign = c("HSPB2", "MYL4", "PIPOX", "TNNT2", "MYOG", "ENO3", "NRN1", "GYPC", "TSPAN3", "TFF3")


library(patchwork)

pdf("Fig3D_dotplot.pdf", width=11, height=8)
## p1=DotPlot(MAST39, features=selection1)+RotatedAxis()
## p2=DotPlot(MAST95, features=selection2)+RotatedAxis()
## p1=Seurat::DotPlot(MAST39, features=erms.topsign, col.min=-3, col.max=3, scale.min=0, scale.max=100, dot.scale=6)+xlab("Modules")+ylab("Genes")+ggtitle('MAST39 ERMS core signatures')+RotatedAxis()
## p2=Seurat::DotPlot(MAST95, features=erms.topsign, col.min=-3, col.max=3, scale.min=0, scale.max=100, dot.scale=6)+xlab("Modules")+ylab("Genes")+ggtitle('MAST95 ERMS core signatures')+RotatedAxis()
## p3=Seurat::DotPlot(MAST39, features=arms.topsign, col.min=-3, col.max=3, scale.min=0, scale.max=100, dot.scale=6)+xlab("Modules")+ylab("Genes")+ggtitle('MAST39 ARMS core signatures')+RotatedAxis()
## p4=Seurat::DotPlot(MAST95, features=arms.topsign, col.min=-3, col.max=3, scale.min=0, scale.max=100, dot.scale=6)+xlab("Modules")+ylab("Genes")+ggtitle('MAST95 ARMS core signatures')+RotatedAxis()
p1=DotPlot(MAST39, features=erms.topsign, col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=6)+xlab("Modules")+ylab("Genes")+ggtitle('MAST39 ERMS core')+RotatedAxis()
p2=DotPlot(MAST95, features=erms.topsign, col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=6)+xlab("Modules")+ylab("Genes")+ggtitle('MAST95 ERMS core')+RotatedAxis()
p3=DotPlot(MAST39, features=arms.topsign, col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=6)+xlab("Modules")+ylab("Genes")+ggtitle('MAST39 ARMS core')+RotatedAxis()
p4=DotPlot(MAST95, features=arms.topsign, col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=6)+xlab("Modules")+ylab("Genes")+ggtitle('MAST95 ARMS core')+RotatedAxis()
(p1 | p2) / (p3 | p4)
dev.off()

pdf("Fig3D_dotplot_supp.pdf", width=35, height=16)
## p1=DotPlot(MAST39, features=selection1)+RotatedAxis()
## p2=DotPlot(MAST95, features=selection2)+RotatedAxis()
p1=DotPlot(MAST39, features=intersect(only.sign[[1]], rownames(MAST39)), col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=5)+xlab("Modules")+ylab("Genes")+ggtitle('MAST39 ERMS core signatures')+RotatedAxis()
p2=DotPlot(MAST39, features=intersect(only.sign[[2]], rownames(MAST39)), col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=5)+xlab("Modules")+ylab("Genes")+ggtitle('MAST39 ARMS core signatures')+RotatedAxis()
p3=DotPlot(MAST95, features=intersect(only.sign[[1]], rownames(MAST95)), col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=5)+xlab("Modules")+ylab("Genes")+ggtitle('MAST95 ERMS core signatures')+RotatedAxis()
p4=DotPlot(MAST95, features=intersect(only.sign[[2]], rownames(MAST95)), col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=5)+xlab("Modules")+ylab("Genes")+ggtitle('MAST95 ARMS core signatures')+RotatedAxis()
p1 + p2 + p3 + p4 + plot_layout(ncol=4, width=c(1,1,1,1))
dev.off()

###############
### Fig 3E variance
###############
## calculate.clustersd <- function(x, area='avg.exp', cutoff, return.sd=T) {
##     x = reshape2::dcast(x, id~features.plot, mean, value.var=area)
##     rownames(x) = x$id
##     x = x[,-1]
##     print(head(x))
##     x.sd = apply(x, 2, sd)
##     if (return.sd) {
##         return(x.sd)
##     } else { 
##         return(sum(x.sd>cutoff)/length(x.sd))
##     }
## }

## compute.sds <- function(x, area='avg.exp', cutoff=0, return.sd=T) {
##     p1=DotPlot(MAST111,   features=x, return_data=T)
##     p2=DotPlot(MAST139,   features=x, return_data=T)
##     p3=DotPlot(MAST39,    features=x, return_data=T)
##     p4=DotPlot(RH74,      features=x, return_data=T)
##     p5=DotPlot(MAST95,    features=x, return_data=T)
##     p5=DotPlot(MAST35,    features=x, return_data=T)
##     p5=DotPlot(MAST85,    features=x, return_data=T)
##     p5=DotPlot(MSK82489,  features=x, return_data=T)
##     sds = c(calculate.clustersd(p1, area, cutoff, return.sd),
##             calculate.clustersd(p2, area, cutoff, return.sd),
##             calculate.clustersd(p3, area, cutoff, return.sd),
##             calculate.clustersd(p4, area, cutoff, return.sd),
##             calculate.clustersd(p5, area, cutoff, return.sd))
##     sds
## }

## ## top 50% percentile
## exp.varsd = rbind(compute.sds(setdiff(only.sign$ERMS.ground, gene.sign.union)),
##                   compute.sds(setdiff(only.sign$ARMS.ground, gene.sign.union)),
##                   compute.sds(gene.sign$hypoxia),
##                   compute.sds(gene.sign$G1S),
##                   compute.sds(gene.sign$G2M))
## exp.varsd = rbind(compute.sds(setdiff(only.sign$ERMS.ground, gene.sign.union), 'avg.exp', quantile(c(exp.varsd), 0.5), F),
##                   compute.sds(setdiff(only.sign$ARMS.ground, gene.sign.union), 'avg.exp', quantile(c(exp.varsd), 0.5), F),
##                   compute.sds(gene.sign$hypoxia, 'avg.exp', quantile(c(exp.varsd), 0.5), F),
##                   compute.sds(gene.sign$G1S, 'avg.exp', quantile(c(exp.varsd), 0.5), F),
##                   compute.sds(gene.sign$G2M, 'avg.exp', quantile(c(exp.varsd), 0.5), F))
## cell.varsd = rbind(compute.sds(setdiff(only.sign$ERMS.ground, gene.sign.union), 'pct.exp'),
##                   compute.sds(setdiff(only.sign$ARMS.ground, gene.sign.union), 'pct.exp'),
##                   compute.sds(gene.sign$hypoxia, 'pct.exp'),
##                   compute.sds(gene.sign$G1S, 'pct.exp'),
##                   compute.sds(gene.sign$G2M, 'pct.exp'))
## cell.varsd = rbind(compute.sds(setdiff(only.sign$ERMS.ground, gene.sign.union), 'pct.exp', quantile(c(cell.varsd), 0.5), F),
##                    compute.sds(setdiff(only.sign$ARMS.ground, gene.sign.union), 'pct.exp', quantile(c(cell.varsd), 0.5), F),
##                    compute.sds(gene.sign$hypoxia, 'pct.exp', quantile(c(cell.varsd), 0.5), F),
##                    compute.sds(gene.sign$G1S, 'pct.exp', quantile(c(cell.varsd), 0.5), F),
##                    compute.sds(gene.sign$G2M, 'pct.exp', quantile(c(cell.varsd), 0.5), F))
## colnames(cell.varsd) = c("MAST111", "MAST139", "MAST39", "RH74", "MAST95")
## rownames(cell.varsd) = c("ERMS core", "ARMS cores", "Hypoxia", "G1S", "G2M")

## colnames(exp.varsd) = c("MAST111", "MAST139", "MAST39", "RH74", "MAST95")
## rownames(exp.varsd) = c("ERMS core", "ARMS cores", "Hypoxia", "G1S", "G2M")
## write.table(cell.varsd, file='cellcluster_ratio.xls', quote=F, sep='\t')
## write.table(exp.varsd, file='expressioncluster_ratio.xls', quote=F, sep='\t')
##     p1=DotPlot(MAST111,   features=x, return_data=T)
##     p2=DotPlot(MAST139,   features=x, return_data=T)
##     p3=DotPlot(MAST39,    features=x, return_data=T)
##     p4=DotPlot(RH74,      features=x, return_data=T)
##     p5=DotPlot(MAST95,    features=x, return_data=T)
##     p5=DotPlot(MAST35,    features=x, return_data=T)
##     p5=DotPlot(MAST85,    features=x, return_data=T)
##     p5=DotPlot(MSK82489,  features=x, return_data=T)

tumor = list()
for (s in list(MAST39, MAST95)) {
    geneset = list(ERMSCore=intersect(setdiff(only.sign$ERMS.ground, gene.sign.union), rownames(s)), ARMSCore=intersect(setdiff(only.sign$ARMS.ground, gene.sign.union), rownames(s)), Hypoxia=intersect(gene.sign$hypoxia, rownames(s)), G1S=intersect(gene.sign$G1S, rownames(s)), G2M=intersect(gene.sign$G2M, rownames(s)))
    tumor[[s@project.name]] = list()
    for (gene in names(geneset)) {
        tumor[[s@project.name]][[gene]] = cbind(tumor=s@project.name, signature=gene, DotPlot(s, features=geneset[[gene]], col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=5, return_data=T))
    }
    tumor[[s@project.name]] = do.call('rbind', tumor[[s@project.name]])
    tumor[[s@project.name]] = reshape2::dcast(tumor[[s@project.name]], signature+features.plot~id, mean, value.var='avg.exp')
}

tumor = do.call('rbind', tumor)

write.table(tumor[[1]], file='Fig3_MAST39_expression.xls', quote=F, col.names=NA, sep='\t')
write.table(tumor[[2]], file='Fig3_MAST95_expression.xls', quote=F, col.names=NA, sep='\t')

## for (cutoff in c(0.2, 0.5, 0.8, 0.9)) {
for (cutoff in c(0.1, 0.2)) {
    varres = list()
    for (s in list(MAST35, MAST85, MAST39, MAST111, MAST139, RH74, MAST95, MSK82489)) {
        varres[[s@project.name]] = list()
        geneset = list(ERMSCore=intersect(setdiff(only.sign$ERMS.ground, gene.sign.union), rownames(s)), ARMSCore=intersect(setdiff(only.sign$ARMS.ground, gene.sign.union), rownames(s)), Hypoxia=intersect(gene.sign$hypoxia, rownames(s)), G1S=intersect(gene.sign$G1S, rownames(s)), G2M=intersect(gene.sign$G2M, rownames(s)))
        for (gene in names(geneset)) {
            ## if (s@project.name%in%c("MAST95", "MSK82489")) {
            ##     if (gene=='ERMSCore')
            ##         next
            ## } else {
            ##     if (gene=='ARMSCore') {
            ##         next
            ##     }
            ## }
            aveexp =AverageExpression(s, features=geneset[[gene]])$RNA
            aveexp.select = apply(aveexp, 1, median)
            print(quantile(aveexp.select, cutoff))
            aveexp = aveexp[aveexp.select >= quantile(aveexp.select, cutoff), ]
            varres[[s@project.name]][[gene]] = apply(aveexp, 1, function(x) {
                sd(x, na.rm=T)
            })
            varres[[s@project.name]][[gene]] = data.frame(tumor=s@project.name, gene=gene, names=names(varres[[s@project.name]][[gene]]), values=varres[[s@project.name]][[gene]])
        }
        varres[[s@project.name]] = do.call(rbind, varres[[s@project.name]])
    }
    varres <- do.call(rbind, varres)
    library(ggridges)
    library(ggplot2)
    ## ggplot(varres, aes(x = values, y = gene, fill=..x..)) +
    ##     geom_density_ridges_gradient(aes(point_shape = gene, point_fill = gene, point_size = values), scale=3, rel_min_height=0.01, quantile_lines=T, quantiles=c(0.05, 0.95), alpha=0.5) + scale_fill_viridis_c(name="Standard Deviation") + facet_wrap(~tumor, ncol=4) + xlim(c(0,3))
    ## ggsave('Fig3_variable_signatures.pdf', width=12, height=8)
    ggplot(varres, aes(x = values, y = gene, fill = gene)) +
        geom_density_ridges(aes(point_shape = gene, point_fill = gene, point_size = values),
                            alpha = 0.5,
                            point_alpha = 1,
                            jittered_points = T,
                            rel_min_height = 0.01,
                            quantile_lines = T,
                            quantiles = c(0.5)) + scale_point_color_hue(l = 40) + scale_point_size_continuous(range = c(0.1, 2)) + scale_discrete_manual(aesthetics="point_shape", values = c(21, 22, 23, 24, 25)) + facet_wrap(~tumor, scales='free_y', ncol=4) + xlim(c(0, 6)) + theme_ridges()
    ## ggsave(paste0('Fig3_variable_signatures_', cutoff, '.pdf'), width=12, height=8)
    ## ggsave(paste0('Fig3_variable_signatures_', cutoff, '_allsign.pdf'), width=12, height=8)
    tumorlist = list()
    for (i in unique(varres$tumor)) {
        z=subset(varres, tumor==i)
        freq = table(z$gene, z$values>=0.5)
        freq = as.matrix(t(t(freq) / rowSums(freq)))
        print(freq)
        tumorlist[[i]] = cbind(tumor=i, freq)
    }
    tumor.sd0.5 = do.call(rbind, tumorlist)
    tumorlist = list()
    for (i in unique(varres$tumor)) {
        z=subset(varres, tumor==i)
        freq = table(z$gene, z$values>=1)
        freq = as.matrix(t(t(freq) / rowSums(freq)))
        print(freq)
        tumorlist[[i]] = cbind(tumor=i, freq)
    }
    tumor.sd1 = do.call(rbind, tumorlist)
    tumor.sd0.5 = cbind(rownames(tumor.sd0.5), tumor.sd0.5)
    tumor.sd0.5 = as.data.frame(tumor.sd0.5)
    colnames(tumor.sd0.5)  = c("Signatures", "Tumor", "InvariableGene", "VariableGene")
    tumor.sd0.5$VariableGene = as.numeric(as.vector(tumor.sd0.5$VariableGene))
    tumor.sd0.5 = reshape2::dcast(tumor.sd0.5, Signatures~Tumor, value.var='VariableGene')
    write.table(tumor.sd0.5, sep='\t', quote=F, file=paste0('tumor_invariable_sd0.5_remove', cutoff, '.xls'), col.names=NA)
    tumor.sd1 = cbind(rownames(tumor.sd1), tumor.sd1)
    tumor.sd1 = as.data.frame(tumor.sd1)
    colnames(tumor.sd1)  = c("Signatures", "Tumor", "InvariableGene", "VariableGene")
    tumor.sd1$VariableGene = as.numeric(as.vector(tumor.sd1$VariableGene))
    tumor.sd1 = reshape2::dcast(tumor.sd1, Signatures~Tumor, value.var='VariableGene')
    write.table(tumor.sd1, sep='\t', quote=F, file=paste0('tumor_invariable_sd1_remove', cutoff, '.xls'), col.names=NA)
}

## Fig1 TF dotplot
tfs= read.table('../lisa_script_modules/Fig1_LISA_TFs.xls', stringsAsFactors=F)
top.tfs = apply(tfs, 1, function(x) colnames(tfs)[order(x, decreasing = T)][1:5])
top.tfs = reshape2::melt(top.tfs)

top.results = list()
for (s in list(MAST35, MAST85, MAST39, MAST111, MAST139, RH74, MAST95, MSK82489)) {
    result.list = list()
    for (i in 1:nrow(top.tfs)) {
        if (length(intersect(rownames(s), as.vector(top.tfs[,3][i])))>0) {
            result.list[[i]] = AverageExpression(s, features=as.vector(top.tfs[,3][i]))$RNA[1,]
        } else {
            result.list[[i]] = rep(0, length(levels(s)))
        }
    }
    top.results[[s@project.name]] = do.call('rbind', result.list)
}


top.results2 = do.call('cbind', top.results)
top.results2 = cbind(top.tfs, top.results2)
annrow = top.results2[,2,drop=F]
colnames(annrow) = c("CellState")
heatdata = as.matrix(top.results2[ ,4:ncol(top.results2)])
rownames(heatdata) = paste0(top.results2[,3], '_', 1:nrow(top.results2))
rownames(annrow) = paste0(top.results2[,3], '_', 1:nrow(top.results2))

anncol = data.frame(TumorCellState=apply(matrix(colnames(top.results2)[4:ncol(top.results2)]), 1, function(x) strsplit(x, '\\.')[[1]][[2]]))
rownames(anncol) = colnames(heatdata)

pheatmap(heatdata, scale="row", cluster_rows=F, cluster_cols=F, fontsize_col=6, fontsize_row=6, show_rownames=T, 
         color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")), bias=1)(100), annotation_row=annrow, annotation_col=anncol,
         filename='Fig1_LISA_TFexp.pdf',
         width=14, height=8)
