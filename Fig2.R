## Fig2A
library(readxl)
library(Seurat)
library(pheatmap)
library(gplots)
library(RColorBrewer)

pdxs = Sys.glob('../data/seurat_obj/*rds')[1:10]
annotation = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1, header=T,
                        check.names=F, stringsAsFactors=F)
mast39 = readRDS(pdxs[4])
mast95 = readRDS(pdxs[6])

mast39$seurat_clusters = mast39$RNA_snn_res.0.8
states = unlist(as.vector(annotation["MAST39", ]))
states = states[states!='']
levels(mast39$RNA_snn_res.0.8) = states

mast95$seurat_clusters = mast95$RNA_snn_res.0.8
states = unlist(as.vector(annotation["MAST95", ]))
states = states[states!='']
levels(mast95$RNA_snn_res.0.8) = states

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

pdf("Fig2A_MAST39.pdf", width=5, height=3.5)
print(DimPlot(mast39, group.by='RNA_snn_res.0.8', cols=metacolors, label=F))
dev.off()

pdf("Fig2A_MAST39_markers.pdf", width=13.5, height=3.5)
FeaturePlot(mast39,
            features=c("MYOD1", "MYOG", "DES"), ncol=3)
dev.off()

pdf("Fig2A_MAST95.pdf", width=5, height=3.5)
print(DimPlot(mast95, group.by='RNA_snn_res.0.8', cols=metacolors, label=F))
dev.off()

pdf("Fig2A_MAST95_markers.pdf", width=13.5, height=3.5)
FeaturePlot(mast95,
            features=c("MYOD1", "MYOG", "DES"), ncol=3)
dev.off()

library(ggplot2)
library(cowplot)
library(patchwork)

pdxs = Sys.glob('../data/seurat_obj/*rds')[1:10]
annotation = read.delim('../final_annotations/Final_clusters.txt', sep='\t', row.names=1,
                        header=T,
                        check.names=F, stringsAsFactors=F)
pdxs = c(pdxs, '../results/seurat_sara/20191031_MSK74711_seurat-object.rds',
         '../results/seurat_sara/MAST118_seurat-object.rds',
         '../results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds',
         '../results/seurat_sara/MAST139_1cells_seurat-object.rds')

pdxs.objs = lapply(pdxs, readRDS)

labels = unlist(lapply(pdxs.objs, function(x) {
    levels(x$orig.ident[1])
}))
labels[9] = 'MAST85-1'
labels[11] = 'MSK74711'
labels[13] = 'MSK72117'
labels[14] = 'MAST139-1'
labels[10] = 'RH74-10'

results = list()
for (i in seq_along(labels)) {
    states = unlist(as.vector(annotation[labels[i], ]))
    states = states[states!='']
    print(states)
    pdxs.objs[[i]]$seurat_clusters = pdxs.objs[[i]]$RNA_snn_res.0.8
    levels(pdxs.objs[[i]]$RNA_snn_res.0.8) = states
    results[[labels[i]]] = table(pdxs.objs[[i]]$RNA_snn_res.0.8)
}

## gene.modules <- Sys.glob('../final_annotations/gene_modules/*txt')
## gene.list <- lapply(gene.modules, scan, what='')
## names(gene.list) <- basename(gsub('.txt', '', gene.modules))
## list2df <- function(x) {
##     ylist <- list()
##     for (y in names(x)) {
##         ylist[[y]] = data.frame(module=y, gene=x[[y]])
##     }
##     do.call('rbind', ylist)
## }
## gene.list <- list2df(gene.list)

## only.sign = read_excel("RMS Core signature gene list.xlsx")
## only.sign = only.sign[,c(1,3)]
## only.sign = lapply(as.list(only.sign), function(x) {na.omit(x[x!=""])})
## df = data.frame(type=as.character(),
##                 gene=as.character())
## for (i in names(only.sign)) {
##     df = rbind(df, cbind(type=i, gene=only.sign[[i]]))
## }
## df = df[!(df$gene%in%as.vector(gene.list[,2])), ]

## refine the core signature by t test
only.sign = read_excel("Differentially Regulated Genes-for Alvin.xlsx")
only.sign = lapply(as.list(only.sign), function(x) {na.omit(x[x!=""])})

df = data.frame(type=as.character(),
                gene=as.character())
for (i in names(only.sign)) {
    df = rbind(df, cbind(type=i, gene=only.sign[[i]]))
}
df = df[!duplicated(df$gene), ]

source("DEGs_seurat3.R")

pdxs.cpms.genes = lapply(pdxs.objs, function(x) {
    intersect(as.vector(df$gene), rownames(x))
})

pdxs.cpms.genes = unique(unlist(pdxs.cpms.genes))
df = df[as.vector(df$gene) %in% pdxs.cpms.genes, ]

pdxs.cpms = lapply(pdxs.objs, function(x) {
    cpm = as.data.frame(apply(as.matrix(x$RNA@counts), 2, correct))
    cpm
})

pdxs.cpms.heatdata = lapply(pdxs.cpms, function(x) {
    y=x[as.vector(df$gene), ]
    ## apply(y, 1, function(x) mean(x[x>0]))
    rowMeans(y)
})

pdxs.cpms.heatdata = do.call(rbind, pdxs.cpms.heatdata)

rownames(pdxs.cpms.heatdata) = labels

labels.show = labels[-grep("-", labels)]
subtype = c("ERMS",
            "ERMS",
            "ERMS",
            "ERMS",
            "ERMS",
            "ARMS",
            "ARMS",
            "ERMS",
            "ERMS",
            "ARMS", "ARMS")

heatdata = t(pdxs.cpms.heatdata)
heatdata = heatdata[-grep("^NA", rownames(heatdata)), ]
heatdata = heatdata[, -grep("-", colnames(heatdata))]
heatdata = heatdata[, order(subtype)]

## df$type = gsub("\\ ground", "",  df$type)
ann_row = data.frame(core=as.vector(df[,1]),
                     row.names=as.vector(df[,2]))

ann_row = ann_row[rownames(heatdata), , drop=F]

heatdata[is.na(heatdata)] = 1e-7

selection = apply(heatdata, 1, function(x) {
    print(x)
    test1 = t.test(x[1:4], x[5:11], alternative='less')
    test2 = t.test(x[1:4], x[5:11], alternative='greater')
    fold1 = mean(x[1:4], na.rm=T) / mean(x[5:11], na.rm=T)
    fold2 = mean(x[5:11], na.rm=T) / mean(x[1:4], na.rm=T)
    return((test1$p.value <= 0.05 | test2$p.value <= 0.05) & (fold1 >= 1.5 | fold2 >= 1.5))
})

write.table(heatdata[selection,], file='RMS_core_t_test_pval0.05_fold1.5.xls', sep='\t', col.names=NA, quote=F)


pheatmap(heatdata[selection, ],
         annotation_col = data.frame(subtype=sort(subtype),
                                     row.names=colnames(heatdata)),
         annotation_row = ann_row[selection, , drop=F],
         annotation_colors = list(core=c("ARMS"='blue', "ERMS"="red"),
                                  subtype=c("ARMS"='blue', "ERMS"="red")), 
         ## scale="row", cluster_rows=F, cluster_cols=F, fontsize_row=0.8,show_rownames=F, breaks=seq(-2, 2, length=99),
         scale="row", cluster_rows=F, cluster_cols=F, fontsize_row=0.8,show_rownames=F, breaks=seq(-1, 1, length=99),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias=0.9)(100),
         filename="Fig2B.tiff",
         width=4, height=8)

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
    ## if (object@project.name == 'MAST39') {
    ##     plot <- plot + scale_x_discrete(limits=c("G1S", "G2M", "PROLIF", "Histones", "Hypoxia", "INTERFERON", "EMT", "GROUND"))
    ## } else {
    ##     plot <- plot + scale_x_discrete(limits=c("G1S", "G2M", "Hypoxia", "MUSCLE", "GROUND", "UNASSIGNED"))
    ## }
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

pdf("Fig2C_dotplot.pdf", width=11, height=8)
p1=DotPlot(mast39, features=erms.topsign, col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=6, group.by='RNA_snn_res.0.8')+xlab("Modules")+ylab("Genes")+ggtitle('MAST39 ERMS core')+RotatedAxis()
p2=DotPlot(mast95, features=erms.topsign, col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=6, group.by='RNA_snn_res.0.8')+xlab("Modules")+ylab("Genes")+ggtitle('MAST95 ERMS core')+RotatedAxis()
p3=DotPlot(mast39, features=arms.topsign, col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=6, group.by='RNA_snn_res.0.8')+xlab("Modules")+ylab("Genes")+ggtitle('MAST39 ARMS core')+RotatedAxis()
p4=DotPlot(mast95, features=arms.topsign, col.min=0, col.max=3, scale.min=0, scale.max=100, dot.scale=6, group.by='RNA_snn_res.0.8')+xlab("Modules")+ylab("Genes")+ggtitle('MAST95 ARMS core')+RotatedAxis()
(p1 | p2) / (p3 | p4)
dev.off()
