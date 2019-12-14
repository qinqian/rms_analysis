library(ggplot2)
library(ggpubr)
library(data.table)

meta = read.table('/PHShome/qq06/pinello/public_datasets/hg38/lisa_meta.xls',
                  sep='\t', header=T, quote="")
dt.samples = list()

for (i in (Sys.glob('../results/gsea/lisa/*chipseq*cauchy*dedup*csv'))) {
    s = substr(basename(i), 1, 30)
    dt.samples[[s]] = list()
}

for (i in (Sys.glob('../results/gsea/lisa/*chipseq*cauchy*dedup*csv'))) {
    cl = gsub('_genes.symbols_chipseq_cauchy_combine_dedup.csv', '',
              gsub('_genes.symbols_motif_cauchy_combine_dedup.csv', '', i))
    cl = unlist(strsplit(cl, '_'))
    cl = cl[length(cl)]
    print(cl)
    s = substr(basename(i), 1, 30)
    dt = fread(i, header=T)
    ids = as.integer(unlist(lapply(unlist(dt[,'0']), function(x) {unlist(strsplit(x, '\\|'))[1]})))
    dt.samples[[s]][[cl]] = as.data.frame(cbind(dt, meta[match(ids, meta$id),]))[1:20,]
    print(i)
}


metacolors = c(rgb(166, 166, 166, maxColorValue = 255),
               rgb(241, 149, 69, maxColorValue = 255),
               rgb(103, 35,  102, maxColorValue = 255),
               rgb(52, 101, 252, maxColorValue = 255),
               rgb(242, 242, 242, maxColorValue = 255),
               rgb(52, 101, 252, maxColorValue = 255),
               rgb(233, 63,  51, maxColorValue = 255),
               rgb(65, 129,  7, maxColorValue = 255),
               rgb(52, 101, 252, maxColorValue = 255),
               rgb(253, 247, 49, maxColorValue = 255))
metalabels = c("GROUND", "Hypoxia", "EMT",  "G1S", "UNASSIGNED", "G2M",  "MUSCLE", "INTERFERON", "PROLIF", "Histones")

color = read.table('color_table.xls',sep='\t', header=T)
labels = color[,7][color[,7]!='']
pdf('MAST111_lisa_results.pdf', width=16, height=12)
par(mfrow=c(3,3), font=2)
for (i in seq_along(labels)) {
    outfig = dt.samples[[names(dt.samples)[2]]][[paste0(i-1, 'cluster')]][,c(1,2,9,10)]
    n = apply(outfig[,c(1,3,4)], 1, paste, sep='_', collapse='_')
    outfig = -log10(outfig[,2])
    names(outfig) = n
    if (i==9) {
        names(outfig)[1] = '36077|MYOD1_rhabdomyosarcoma_Muscle'
    }
    x = barplot(outfig, main=paste0(i-1, 'cluster', ' ', labels[i]),
                col=metacolors[match(labels[i], metalabels)],
                xaxt="n")
    if (labels[i]=='EMT') {
        col='red'
    } else {
        col='black'
    }
    text(cex=1, x=x-.12, y=0, names(outfig), adj=0, xpd=TRUE, srt=90, col=col)
}
dev.off()


labels = color[,8][color[,8]!='']
pdf('MAST139_lisa_results.pdf', width=16, height=12)
par(mfrow=c(3,3), font=2)
for (i in seq_along(labels)) {
    outfig = dt.samples[[names(dt.samples)[3]]][[paste0(i-1, 'cluster')]][,c(1,2,9,10)]
    if (is.null(outfig)) {
        next
    }
    print(outfig)
    n = apply(outfig[,c(1,3,4)], 1, paste, sep='_', collapse='_')
    outfig = -log10(outfig[,2])
    names(outfig) = n
    x = barplot(outfig, main=paste0(i-1, 'cluster', ' ', labels[i]), col=metacolors[match(labels[i], metalabels)],
            xaxt="n")
    if (labels[i]=='EMT') {
        col='red'
    } else {
        col='black'
    }
    text(cex=1, x=x-.12, y=0, names(outfig), adj=0, xpd=TRUE, srt=90, col=col)
}
dev.off()


labels = color[,2][color[,2]!='']
pdf('RH74_lisa_results.pdf', width=16, height=12)
par(mfrow=c(3,3), font=2)
for (i in seq_along(labels)) {
    outfig = dt.samples[[names(dt.samples)[9]]][[paste0(i-1, 'cluster')]][,c(1,2,9,10)]
    if (is.null(outfig)) {
        next
    }
    print(outfig)
    n = apply(outfig[,c(1,3,4)], 1, paste, sep='_', collapse='_')
    outfig = -log10(outfig[,2])
    names(outfig) = n
    x = barplot(outfig, main=paste0(i-1, 'cluster', ' ', labels[i]), col=metacolors[match(labels[i], metalabels)],
            xaxt="n")
    if (labels[i]=='EMT') {
        col='red'
    } else {
        col='black'
    }
    text(cex=1, x=x-.12, y=0, names(outfig), adj=0, xpd=TRUE, srt=90, col=col)
}
dev.off()


labels = color[,4][color[,4]!='']
pdf('MAST39_lisa_results.pdf', width=16, height=12)
par(mfrow=c(3,3), font=2)
for (i in seq_along(labels)) {
    outfig = dt.samples[[names(dt.samples)[5]]][[paste0(i-1, 'cluster')]][,c(1,2,9,10)]
    if (is.null(outfig)) {
        next
    }
    print(outfig)
    n = apply(outfig[,c(1,3,4)], 1, paste, sep='_', collapse='_')
    outfig = -log10(outfig[,2])
    names(outfig) = n
    x = barplot(outfig, main=paste0(i-1, 'cluster', ' ', labels[i]), col=metacolors[match(labels[i], metalabels)],
            xaxt="n")
    if (labels[i]=='EMT') {
        col='red'
    } else {
        col='black'
    }
    text(cex=1, x=x-.12, y=0, names(outfig), adj=0, xpd=TRUE, srt=90, col=col)
}
dev.off()


labels = color[,10][color[,10]!='']
pdf('MAST95_lisa_results.pdf', width=16, height=12)
par(mfrow=c(3,3), font=2)
for (i in seq_along(labels)) {
    outfig = dt.samples[[names(dt.samples)[7]]][[paste0(i-1, 'cluster')]][,c(1,2,9,10)]
    if (is.null(outfig)) {
        next
    }
    print(outfig)
    n = apply(outfig[,c(1,3,4)], 1, paste, sep='_', collapse='_')
    outfig = -log10(outfig[,2])
    names(outfig) = n
    x = barplot(outfig, main=paste0(i-1, 'cluster', ' ', labels[i]), col=metacolors[match(labels[i], metalabels)],
            xaxt="n")
    if (labels[i]=='EMT') {
        col='red'
    } else {
        col='black'
    }
    text(cex=1, x=x-.12, y=0, names(outfig), adj=0, xpd=TRUE, srt=90, col=col)
}
dev.off()


labels = color[,5][color[,5]!='']
pdf('MAST85_lisa_results.pdf', width=16, height=12)
par(mfrow=c(3,3), font=2)
for (i in seq_along(labels)) {
    outfig = dt.samples[[names(dt.samples)[6]]][[paste0(i-1, 'cluster')]][,c(1,2,9,10)]
    if (is.null(outfig)) {
        next
    }
    print(outfig)
    n = apply(outfig[,c(1,3,4)], 1, paste, sep='_', collapse='_')
    outfig = -log10(outfig[,2])
    names(outfig) = n
    x = barplot(outfig, main=paste0(i-1, 'cluster', ' ', labels[i]), col=metacolors[match(labels[i], metalabels)],
            xaxt="n")
    if (labels[i]=='EMT') {
        col='red'
    } else {
        col='black'
    }
    text(cex=1, x=x-.12, y=0, names(outfig), adj=0, xpd=TRUE, srt=90, col=col)
}
dev.off()

labels = color[,9][color[,9]!='']
pdf('MAST35_lisa_results.pdf', width=16, height=12)
par(mfrow=c(3,3), font=2)
for (i in seq_along(labels)) {
    outfig = dt.samples[[names(dt.samples)[4]]][[paste0(i-1, 'cluster')]][,c(1,2,9,10)]
    if (is.null(outfig)) {
        next
    }
    print(outfig)
    n = apply(outfig[,c(1,3,4)], 1, paste, sep='_', collapse='_')
    outfig = -log10(outfig[,2])
    names(outfig) = n
    x = barplot(outfig, main=paste0(i-1, 'cluster', ' ', labels[i]), col=metacolors[match(labels[i], metalabels)],
            xaxt="n")
    if (labels[i]=='EMT') {
        col='red'
    } else {
        col='black'
    }
    text(cex=1, x=x-.12, y=0, names(outfig), adj=0, xpd=TRUE, srt=90, col=col)
}
dev.off()


labels = color[,11][color[,11]!='']
pdf('MSK82489_lisa_results.pdf', width=16, height=12)
par(mfrow=c(3,3), font=2)
for (i in seq_along(labels)) {
    outfig = dt.samples[[names(dt.samples)[8]]][[paste0(i-1, 'cluster')]][,c(1,2,9,10)]
    if (is.null(outfig)) {
        next
    }
    print(outfig)
    n = apply(outfig[,c(1,3,4)], 1, paste, sep='_', collapse='_')
    outfig = -log10(outfig[,2])
    names(outfig) = n
    x = barplot(outfig, main=paste0(i-1, 'cluster', ' ', labels[i]), col=metacolors[match(labels[i], metalabels)],
            xaxt="n")
    if (labels[i]=='EMT') {
        col='red'
    } else {
        col='black'
    }
    text(cex=1, x=x-.12, y=0, names(outfig), adj=0, xpd=TRUE, srt=90, col=col)
}
dev.off()


