library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
source('functions.R')
set.seed(100)
library(loomR)

get_args = function(x) {
    require(argparse)
    parser = ArgumentParser(description='seurat normalization')
    parser$add_argument('--seuratobj', dest='seurat', default='')
    parser$add_argument('--velobj', dest='vel', default='')
    parser$add_argument('--label', dest='label', default='')
    parser$add_argument('--finalres', dest='res', default=0.8, type='double')
    parser$add_argument('--species', dest='species', type='character', default='fish')
    args = parser$parse_args()
    args
}

args = get_args()
if (args$seurat == '' || args$label == '' | args$vel == '') {
    cat('empty argument, exit..')
    q()
}

## args$vel = '../results/seurat/Tumor24_velocity_seurat_obj_tumors.rds'
## args$seurat = '../results/seurat/Tumor24_seurat_obj_tumors.rds'
## args$label = 'Tumor24'
## args$vel = '../results/seurat/Tumor21_velocity_seurat_obj_tumors.rds'
## args$seurat = '../results/seurat/Tumor21_seurat_obj_tumors.rds'
## args$label = 'Tumor21'

args$vel = 'Differentiated_seurat_obj_tumors.rds'
args$seurat = 'results/seurat_sara/Differentiated_seurat-object.rds'
args$label = 'Differentiated'
args$species = 'human'

vel = readRDS(args$vel)
seu = readRDS(args$seurat)

print(names(vel@assays))
print(names(seu@assays))

if (args$species == 'fish') {
    if (args$label != 'Tumor24') {
        vel.cells = gsub('_.+:', '_', gsub('x', '', rownames(vel@meta.data)))
    } else {
        vel.cells = gsub('.+:', '', gsub('x', '', rownames(vel@meta.data)))
    }
} else {
    vel.cells = gsub('.+:', '', gsub('x', '', rownames(vel@meta.data)))
    # for primary tumor
    vel.cells = gsub('.*_hg19_premrna:', '', vel.cells)
    vel.cells = gsub('.*_hg19:', '', vel.cells)
    #vel.cells = gsub('.*_hg19_premrna:', '', gsub('x', '', rownames(vel@meta.data)))
    print(head(vel.cells))
}

seu.cells = rownames(seu@meta.data)
print(length(vel.cells))
print(length(seu.cells))

print('--------test------------')
intersect.cells1 = sort(intersect(vel.cells, seu.cells))
if (length(intersect.cells1) == length(vel.cells) || length(intersect.cells1) == length(seu.cells))
   intersect.cells1 = intersect.cells1[-1]
vel.cells = colnames(vel)[vel.cells%in%intersect.cells1]

## vel.cells_order = order(gsub('.*_hg19_premrna:', '', gsub('x', '', vel.cells)))
vel.cells_order = order(gsub('.+:', '', gsub('x', '', vel.cells)))

print('--------test------------')
print(head(vel.cells))
print(head(vel.cells_order))
print('--------test------------')

## cell index order problems https://github.com/satijalab/seurat/issues/1492
## update to latest to avoid this errors
## use devtools to install
## see https://github.com/satijalab/seurat/blob/b51801bc4b1a66aed5456473c9fe0be884994c93/R/objects.R#L4878
## seu2 = Seurat:::subset.Seurat(seu, cells=sort(intersect.cells1))
seu = subset(seu, cells=sort(intersect.cells1))
vel = subset(vel, cells=vel.cells[vel.cells_order])

print(dim(seu))
print(dim(vel))

if (args$species == 'fish') {
    if (args$label != 'Tumor24') {
        ncol(seu) == sum(gsub('_.+:', '_', gsub('x', '', rownames(vel@meta.data))) == rownames(seu@meta.data))
    } else {
        ncol(seu) == sum(gsub('.+:', '', gsub('x', '', rownames(vel@meta.data))) == rownames(seu@meta.data))
    }
} else {
    #print(ncol(seu) == sum(gsub('.*:', '', gsub('x', '', rownames(vel@meta.data))) == rownames(seu@meta.data)))
    checkcell = gsub('.+:', '', gsub('x', '', rownames(vel@meta.data)))
    checkcell = gsub('.*_hg19_premrna:', '', checkcell)
    checkcell = gsub('.*_hg19:', '', checkcell)
    print(head(checkcell))
    print(head(rownames(seu@meta.data)))
    print(ncol(seu) == sum(checkcell == rownames(seu@meta.data)))
}

#quit()

## ignore cluster 
if (args$species == 'fish') {
  seu <- recluster.withtree(seu, name=args$label)
  seu <- FindClusters(object=seu, resolution=args$res)
} else {
  Idents(seu) = seu@meta.data$seurat_clusters = seu$RNA_snn_res.0.8
}

## forcely use sara's cluster labels to colorize the plot
#seu.markers = FindAllMarkers(seu, only.pos=T, min.pct=0.1,
#                             test.use='MAST',
#                             ## assay='SCT', slot='data', #slot='scale.data',
#                             random.seed=100, logfc.threshold = 0.1)
#human_ortholog = read.table('/data/langenau/alvin_singlecell/01_rms_projects//01_fish/data/ortholog_mapping/Beagle_fish_human_all_genes.txt', header=T, sep='\t', stringsAsFactors=F)
#seu.markers    = cbind(seu.markers, human_ortholog[match(seu.markers$gene, human_ortholog$Gene), ])

print('recluster....')
## vel <- recluster.withtree(vel, name=args$label)
## vel <- FindClusters(object=vel, resolution=args$res)
#vel.markers = FindAllMarkers(vel, only.pos=T, min.pct=0.1,
#                             test.use='MAST',
#                             assay='SCT', slot='data', #slot='scale.data',
#                             random.seed=100, logfc.threshold = 0.1)

system('mkdir -p ../results/seurat_intersect_velocity')

## it's fine for seurat object
saveRDS(seu, file=paste0('../results/seurat_intersect_velocity/', args$label, '_seu.rds'))
#write.table(seu.markers, file=paste0('../results/seurat_intersect_velocity/', args$label, paste0('_seu_markers_tumoronly_res', args$res, '.xls')), sep='\t', quote=F)

## this is for Seurat Wrapper of velocity
saveRDS(vel, file=paste0('../results/seurat_intersect_velocity/', args$label, '_vel.rds'))
#write.table(vel.markers, file=paste0('../results/seurat_intersect_velocity/', args$label, paste0('_vel_markers_tumoronly_res', args$res, '.xls')), sep='\t', quote=F)

## fill in spliced/unspliced information for velocity object for scvelo analysis
## vel.loom = as.loom(vel, filename=paste0('../results/seurat_intersect_velocity/', args$label, '_vel.loom'), verbose=T)

## print('test------------;')
## ## vel.loom$add.layer(
## ##             layers = list(
## ##                 'spliced' = as.matrix(
## ##                     x = t(
## ##                        x = as.data.frame(
## ##                             vel$spliced@data
## ##                         )[rownames(vel$SCT), colnames(vel$SCT)]
## ##                     )
## ##                 )
## ##             ),
## ##             verbose = T
## ##         )
## #vel.loom$add.layer(
## #             layers = list(
## #                 'unspliced' = as.matrix(
## #                     x = t(
## #                         x = as.data.frame(
## #                             vel$unspliced@data
## #                         )[rownames(vel$SCT), colnames(vel$SCT)]
## #                     )
## #                 )
## #             ),
## #             verbose = T
## #         )

## vel.loom$add.layer = list(
##                 'spliced' = as.matrix(
##                     x = t(
##                        x = as.data.frame(
##                             vel$spliced@data
##                         )[rownames(vel$SCT), colnames(vel$SCT)]
##                     )
##                 ),
##                 'unspliced' = as.matrix(
##                     x = t(
##                         x = as.data.frame(
##                             vel$unspliced@data
##                         )[rownames(vel$SCT), colnames(vel$SCT)]
##                     )
##                 )
##             )
## vel.loom$close_all()
