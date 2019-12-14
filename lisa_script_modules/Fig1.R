library(ggplot2)
library(ggpubr)
library(data.table)
library(Rcpp)
#> 
#> Attaching package: 'Rcpp'
#> The following object is masked from 'package:inline':
#> 
#>     registerPlugin
cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}')

meta = read.table('/PHShome/qq06/pinello/public_datasets/hg38/lisa_meta.xls',
                  sep='\t', header=T, quote="")

## metalabels = c("GROUND", "Hypoxia", "EMT",  "G1S", "UNASSIGNED", "G2M",  "MUSCLE", "INTERFERON", "PROLIF", "Histones")
metacolors = c(EMT=rgb(103, 35,  102, maxColorValue = 255),
               G1S=rgb(52, 101, 252, maxColorValue = 255),
               G2M=rgb(52, 101, 252, maxColorValue = 255),
               hypoxia=rgb(241, 149, 69, maxColorValue = 255),
               ifng=rgb(241, 149, 69, maxColorValue = 255),
               TNFA=rgb(241, 149, 69, maxColorValue = 255),
               muscles=rgb(233, 63,  51, maxColorValue = 255),
               histone='#FDF731')

dt.samples <- Sys.glob("*chipseq*cauchy*dedup*csv")
dt.results <- list()
for (d in dt.samples) {
    results = fread(d, header=T)
    results.out = -log10(as.data.frame(results[,2][1:5])[,1])
    names(results.out) = results[,TF][1:5]
    dt.results[[unlist(strsplit(d, '\\.'))[1]]] = data.frame(unlist(strsplit(d, '\\.'))[1], names(results.out), results.out)
}

library(ggpubr)
dt.results = Reduce(rbind, dt.results)
dt.results[, 3] = as.numeric(dt.results[, 3])
colnames(dt.results) = c("CellState", "TF", "p_value")

write.table(dt.results, sep='\t', quote=F, file='40TFs_cellstates.xls', col.names=NA)
##https://stackoverflow.com/questions/37587288/ggplot-label-bars-in-grouped-bar-plot
p <- ggplot(dt.results, aes(x=TF, y=p_value, fill=CellState)) + geom_col() + scale_fill_manual(values=metacolors) + facet_wrap(CellState~., scales = "free", ncol=4)
pdf('Fig1_supp.pdf', width=18, height=8)
print(p)
dev.off()
## p + theme_pubr() + theme(legend.position='none') 


dt.samples <- Sys.glob("*chipseq*cauchy*dedup*csv")
dt.results <- list()
for (d in dt.samples) {
    results = fread(d, header=T)
    results.out = -log10(as.data.frame(results[,2][1:5])[,1])
    dt.results[[unlist(strsplit(d, '\\.'))[1]]] = results[,TF][1:5]
}

all.tfs <- unlist(dt.results)
dt.results <- list()
for (d in dt.samples) {
    results = fread(d, header=T)
    results.sub = results[results[,TF%in%all.tfs],]
    results.sub = as.data.frame(results.sub)
    results.sub = results.sub[order(results.sub[,2]),]
    results.sub[,2] = -log10(results.sub[,2])
    x=results.sub[,2]
    names(x) = results.sub[,3]
    x = x[all.tfs]
    dt.results[[unlist(strsplit(d, '\\.'))[1]]] = x
}

labels = names(dt.results)

dt.results = Reduce(rbind, dt.results)
rownames(dt.results) = labels

write.table(dt.results, file="Fig1_LISA_TFs.xls", quote=F, sep='\t', col.names=NA)

library(pheatmap)
library(gplots)
library(RColorBrewer)

pheatmap(as.matrix(dt.results), cluster_rows=F, cluster_cols=F, fontsize_col=6, fontsize_row=6, breaks=seq(0, 20, length=49),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")), bias=1)(50),
         filename="Fig1_LISA.pdf",
         width=7.5, height=3.0)
