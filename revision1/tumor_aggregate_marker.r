library(Seurat)
require(argparse)
library(glue)

get_args = function(x) {
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--seuratobj', dest='seurat', metavar='N', type="character", nargs='+')
    parser$add_argument('--label', dest='label', metavar='N', type="character", nargs="+")
    parser$add_argument('--state', dest='state', metavar='N', type="character", nargs="+")
    args = parser$parse_args()
    args
}

args = get_args()

pdxs.objs = readRDS(args$seurat)
label = args$label
state = args$state
annotation = read.delim('~/langenau/projects/01_sc_rms/final_annotations/Final_clusters.txt', sep='\t', row.names=1, header=T, check.names=F, stringsAsFactors=F)

states = unlist(as.vector(annotation[label, ]))
print(states)
states = states[states!='']
levels(pdxs.objs$RNA_snn_res.0.8) = states
Idents(pdxs.objs) = pdxs.objs$RNA_snn_res.0.8

de = FindMarkers(pdxs.objs, only.pos=F, min.pct=0,
                 ident.1=state,
                 min.diff.pct=0,
                 random.seed=100, logfc.threshold = 0)
write.csv(de, file=glue('{label}_{state}_DE_genes.csv'), row.names=T, quote=F)

