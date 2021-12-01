library(dplyr)
library(GetoptLong)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--cluster', dest='cluster')
    parser$add_argument('--label', dest='label')
    args = parser$parse_args()
    args
}

args=get_args()

## df = read.table('Regular_1_res1_cells_number.xls')
df = read.table(args$cluster)

state = read.delim("LARRY2_cellstate_res1.txt")
rownames(state) = state[,1]
state = state[,-1]
state = unlist(state[args$label,,drop=T])
state = state[state!='']

res = cbind(state, df) %>% group_by(state) %>% summarise(cell_number=sum(Freq))
readr::write_csv(res, file=paste0(args$label, '_agg_res.xls'))
