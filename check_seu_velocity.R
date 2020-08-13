library(Seurat)

seu = Sys.glob('../results/seurat_intersect_velocity/*seu.rds')
vel = Sys.glob('../results/seurat_intersect_velocity/*vel.rds')

for (i in seq_along(seu)) {
    a=readRDS(seu[i])
    b=readRDS(vel[i])
    a.cells = colnames(a)
    a.cells = gsub('.+_', '', gsub('.+:', '', gsub('x', '', a.cells)))
    b.cells = gsub('.+_', '', gsub('.+:', '', gsub('x', '', colnames(b))))
    cat(seu[i], ' ', sum(a.cells==b.cells), '\n')
}
