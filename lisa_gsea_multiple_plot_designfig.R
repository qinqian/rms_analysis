library(ggplot2)
library(ggpubr)
library(data.table)

meta = read.table('/PHShome/qq06/pinello/public_datasets/hg38/lisa_meta.xls',
                  sep='\t', header=T, quote="")

logs = Sys.glob('../lisa_mergedcluster_gene_set2/*logs')
logs.f=apply(matrix(logs), 1, function(x) {return(length(readLines(x)))})

logs = logs[logs.f>16]

labels = gsub('_lisa_run.logs', '', basename(logs))

dt.samples = list()
for (l in labels) {
    dt.samples[[l]] = Sys.glob(paste0('../lisa_mergedcluster_gene_set2/', l, '*chipseq*cauchy*dedup*csv'))[1]
    dt = fread(dt.samples[[l]], header=T)
    ids = as.integer(unlist(lapply(unlist(dt[,'0']), function(x) {unlist(strsplit(x, '\\|'))[1]})))
    dt.samples[[l]] = as.data.frame(cbind(dt, meta[match(ids, meta$id),]))[1:25,]
}

metalabels = c("GROUND", "Hypoxia", "EMT",  "G1S", "UNASSIGNED", "G2M",  "MUSCLE", "INTERFERON", "PROLIF", "Histones")
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

results.all = list()
results.ams = list()
results.ems = list()
for (m in metalabels) {
    tf.ams = c()
    tf.ems = c()
    tf.union = c()
    for (i in grep(m, names(dt.samples))) {
       print(length(grep(unlist(strsplit(names(dt.samples)[i], '_'))[1], c("MAST35", "MSK82489", "MAST95"))))
        if (length(grep(unlist(strsplit(names(dt.samples)[i], '_'))[1], c("MAST35", "MSK82489", "MAST95")))==0){
            tf.ems = c(tf.ems, dt.samples[[i]]$TF)
        } else {
            print(tf.ams)
            tf.ams = c(tf.ams, dt.samples[[i]]$TF)
        }
        tf.union = c(tf.union, dt.samples[[i]]$TF)
    }
    tf.ems = sort(table(tf.ems), decreasing=T)
    if (length(tf.ems) != 0) {
        results.ems[[m]] = cbind(m, as.data.frame(tf.ems[1:5]))
    }
    tf.ams = sort(table(tf.ams), decreasing=T)
    if (length(tf.ams) != 0) {
        results.ams[[m]] = cbind(m, as.data.frame(tf.ams[1:5]))
    }
    tf.union = sort(table(tf.union), decreasing=T)
    results.all[[m]] = cbind(m, as.data.frame(tf.union[1:5]))
}

library(ggpubr)

results.all = Reduce(rbind, results.all)
results.ems = Reduce(rbind, results.ems)
results.ams = Reduce(rbind, results.ams)
colnames(results.all) = c("CellState", "TF", "Frequency")
colnames(results.ams) = c("CellState", "TF", "Frequency")
colnames(results.ems) = c("CellState", "TF", "Frequency")

##https://stackoverflow.com/questions/37587288/ggplot-label-bars-in-grouped-bar-plot
p=ggbarplot(results.all, x="CellState", y="Frequency", color="TF", position=position_dodge())
ggsave('Fig4B.pdf', width=18, height=6)

p <- ggplot(results.all, aes(x=CellState, y=Frequency)) + geom_bar(aes(fill = TF), position = position_dodge(0.9), stat="identity")
p+theme(legend.position="none")+scale_x_discrete(limits=unique(results.all$CellState)) + geom_text(position = position_dodge(0.9), aes(y=Frequency+0.25, fill=TF, label=TF, hjust=0), angle=90)+ylim(c(0,15))+theme_pubr()
ggsave('Fig4B_2.pdf', width=18, height=8)

results.ams[,3] = -results.ams[,3]

## https://stackoverflow.com/questions/13734368/ggplot2-and-a-stacked-bar-chart-with-negative-values
p <- ggplot() + geom_bar(data=results.ems, aes(x=CellState, y=Frequency, fill = TF), position = position_dodge(0.9), stat="identity")
p <- p+geom_bar(data=results.ams, aes(x=CellState, y=Frequency, fill = TF), position = position_dodge(0.9), stat="identity")
print(p+theme(legend.position='none')+theme_pubr())+ylim(c(-7,7))
ggsave('Fig4B_3.pdf', width=18, height=8)
