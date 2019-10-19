library(Seurat)
library(SeuratWrappers)
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

#args$vel = '../results/seurat/Tumor24_velocity_seurat_obj_tumors.rds'
#args$seurat = '../results/seurat/Tumor24_seurat_obj_tumors.rds'
#args$label = 'Tumor24'
## args$vel = '../results/seurat/Tumor21_velocity_seurat_obj_tumors.rds'
## args$seurat = '../results/seurat/Tumor21_seurat_obj_tumors.rds'
## args$label = 'Tumor21'

vel = readRDS(args$vel)
seu = readRDS(args$seurat)
print(dim(vel))
print(dim(seu))
print(head(vel@meta.data))
print(head(seu@meta.data))

if (args$species == 'fish') {
    if (args$label != 'Tumor24') {
        vel.cells = gsub('_.+:', '_', gsub('x', '', rownames(vel@meta.data)))
    } else {
        vel.cells = gsub('.+:', '', gsub('x', '', rownames(vel@meta.data)))
    }
} else {
    vel.cells = gsub('.*_hg19:', '', gsub('x', '', rownames(vel@meta.data)))
}

seu.cells = rownames(seu@meta.data)

print(head(vel.cells))
print(head(seu.cells))

intersect.cells1 = intersect(vel.cells, seu.cells)
vel.cells = colnames(vel)[vel.cells%in%intersect.cells1]
vel.cells_order = order(gsub('.*_hg19:', '', gsub('x', '', vel.cells)))

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
    print(ncol(seu) == sum(gsub('.*_hg19:', '', gsub('x', '', rownames(vel@meta.data))) == rownames(seu@meta.data)))
}

## ignore cluster 
if (args$species == 'fish') {
  seu <- recluster.withtree(seu, name=args$label)
  seu <- FindClusters(object=seu, resolution=args$res)
} else {
  Idents(seu) = seu@meta.data$seurat_clusters = seu$RNA_snn_res.0.8
}

## forcely use sara's cluster labels to colorize the plot
seu.markers = FindAllMarkers(seu, only.pos=T, min.pct=0.1,
                             test.use='MAST',
                             ## assay='SCT', slot='data', #slot='scale.data',
                             random.seed=100, logfc.threshold = 0.1)

vel <- recluster.withtree(vel, name=args$label)
vel <- FindClusters(object=vel, resolution=args$res)
vel.markers = FindAllMarkers(vel, only.pos=T, min.pct=0.1,
                             test.use='MAST',
                             assay='SCT', slot='data', #slot='scale.data',
                             random.seed=100, logfc.threshold = 0.1)

system('mkdir -p ../results/seurat_intersect_velocity')

## it's fine for seurat object
saveRDS(seu, file=paste0('../results/seurat_intersect_velocity/', args$label, '_seu.rds'))
write.table(seu.markers, file=paste0('../results/seurat_intersect_velocity/', args$label, paste0('_seu_markers_tumoronly_res', args$res, '.xls')), sep='\t', quote=F)

## this is for Seurat Wrapper of velocity
saveRDS(vel, file=paste0('../results/seurat_intersect_velocity/', args$label, '_vel.rds'))
write.table(vel.markers, file=paste0('../results/seurat_intersect_velocity/', args$label, paste0('_vel_markers_tumoronly_res', args$res, '.xls')), sep='\t', quote=F)

## fill in spliced/unspliced information for velocity object for scvelo analysis
vel.loom = as.loom(vel, filename=paste0('../results/seurat_intersect_velocity/', args$label, '_vel.loom'), verbose=T)
vel.loom$add.layer(
             layers = list(
                 'spliced' = as.matrix(
                     x = t(
                         x = as.data.frame(
                             vel$spliced@data
                         )[rownames(vel$SCT), colnames(vel$SCT)]
                     )
                 )
             ),
             verbose = T
         )
vel.loom$add.layer(
             layers = list(
                 'unspliced' = as.matrix(
                     x = t(
                         x = as.data.frame(
                             vel$unspliced@data
                         )[rownames(vel$SCT), colnames(vel$SCT)]
                     )
                 )
             ),
             verbose = T
         )
vel.loom$close_all()
