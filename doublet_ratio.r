library(readr)
for (i in Sys.glob('../results/doublets/2*_doublet_doublet.csv')) {
    df = read.csv(i)
    print(i)
    print(table(df[,4])/nrow(df))
}
