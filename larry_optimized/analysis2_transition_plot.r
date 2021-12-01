library(patchwork)
library(ggplot2)
library(Gmisc)
library(reshape2)
## install.packages("diagram")
library(diagram)

dfs = list()

for (fi in c("Regular_1_vs_Regular_2.xls", "Regular_1_vs_Diff_3.xls",
             "Diff_3_vs_Diff_4.xls")) {
df = read.table(fi)
transition = matrix(0, nrow=4, ncol=4)
colnames(transition) = rownames(transition) = c("CSC", "P", "G", "M")
for (n in 1:nrow(df)) {
    progeny = na.omit(unlist(df[n, grepl("progeny", colnames(df)), drop=T]))
    parent = na.omit(unlist(df[n, grepl("parent", colnames(df)), drop=T]))
    progeny[grepl("EMT", progeny)] = 'CSC'
    progeny[grepl("Prolif", progeny)] = 'P'
    progeny[grepl("Ground", progeny)] = 'G'
    progeny[grepl("Muscle", progeny)] = 'M'
    parent[grepl("EMT", parent)] = 'CSC'
    parent[grepl("Prolif", parent)] = 'P'
    parent[grepl("Ground", parent)] = 'G'
    parent[grepl("Muscle", parent)] = 'M'
    if (unique(parent) == 'Interferon') next
    for (progeny_fate in progeny) {
        if (progeny_fate=='Interferon')
            next
        transition[unique(parent), progeny_fate] <- transition[unique(parent), progeny_fate]+1
    }
}
df.melt = melt(transition)
print(transition)
dfs[[fi]] = transition
write.table(dfs[[fi]], file=paste0(gsub(".xls", "", fi), '_transition.xls'), quote=F, sep='\t')
transition = transition / rowSums(transition)
pdf(gsub(".xls", ".pdf", fi), width=18, height=6)
plotmat(t(transition),
        pos=4,
        lwd =1, t(transition),
        box.lwd = 1,
        cex.txt = 0.8,
        box.size = 0.05,
        box.prop = 0.5,
        main = gsub(".xls", "", fi))
## p1=transitionPlot(dfs[[fi]],
##                   box_txt = rownames(transition),
##                   type_of_arrow = "simple",
##                   min_lwd = unit(1, "mm"),
##                   max_lwd = unit(6, "mm"),
##                   overlap_add_width = unit(0.2, "mm"))
df.melt.prop = melt(transition)
p1 = ggplot(df.melt.prop, aes(fill=Var2, y=Var1, x=value)) + 
    geom_bar(position="fill", stat='identity')+ggtitle("Analysis 2: progeny cell fate")+xlab("Proportion of progeny cell state")+ylab("Parent cell state")
p2 = ggplot(df.melt, aes(fill=Var2, y=Var1, x=value)) + 
    geom_col()+xlab("Frequency of progeny cell state")+ylab("Parent cell state")
print(p1+p2)
dev.off()
}
