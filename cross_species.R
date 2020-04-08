library(Seurat)
objs = Sys.glob("/data/langenau/human_rms_pdxs/seurat_objects/*")

colortab = read.table('color_table.xls', sep='\t', header=T, check.names=F)
fishcolortab = read.delim('fish_color_table.txt', sep='\t', header=T, check.names=F)

human_ortholog = read.table('~/langenau/01_rms_projects/01_fish/data/ortholog_mapping/Beagle_fish_human_all_genes.txt', header=T, sep='\t', stringsAsFactors=F)

human = readRDS(objs[2])

fish = readRDS("../results/seurat_intersect_velocity/Tumor21_seu.rds")

## human_ortholog = human_ortholog[(!(is.na(human_ortholog$Hsortholog))) | (!(is.na(human_ortholog$Hs_closest))),]
human_ortholog = human_ortholog[(!(is.na(human_ortholog$Hsortholog))), ]
human_ortholog = human_ortholog[human_ortholog$Hsortholog!="", ]
human_ortholog = human_ortholog[human_ortholog$Hsortholog %in% rownames(human), ]
human_ortholog = human_ortholog[human_ortholog$Gene %in% rownames(fish), ]

human.sub = human$RNA@counts[human_ortholog$Hsortholog, ]

human.obj <- CreateSeuratObject(counts=human.sub)

labels = as.character(colortab[, 'MAST139'])
levels(human$RNA_snn_res.0.8) = labels[labels!=""]
human.obj$seurat_clusters = human$RNA_snn_res.0.8

labels = as.character(fishcolortab[, 'Tumor21'])

levels(fish$seurat_clusters) = labels[labels!=""]

fish.sub = fish$RNA@counts[human_ortholog$Gene, ]
rownames(fish.sub) = human_ortholog$Hsortholog
fish.obj <- CreateSeuratObject(counts=fish.sub)
fish.obj$seurat_clusters = fish$seurat_clusters

integration <- function(x, y, label, method='seurat') {
    seurat.pseudo.list = list(x, y)
    for (i in 1:length(seurat.pseudo.list)) {
        print(1)
        seurat.pseudo.list[[i]] <- NormalizeData(seurat.pseudo.list[[i]], verbose = FALSE)
        seurat.pseudo.list[[i]] <- FindVariableFeatures(seurat.pseudo.list[[i]],
                                                        selection.method = "vst", 
                                                        nfeatures = 2000, verbose = FALSE)
    }
    names(seurat.pseudo.list) <- label
    if (method == 'seurat') {
        anchors <- FindIntegrationAnchors(object.list = seurat.pseudo.list, dims = 1:30)
        integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
        DefaultAssay(integrated) <- "integrated"
        integrated <- ScaleData(integrated, verbose = FALSE)
        integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
    } else if (method == 'conos') {
        seurat.pseudo.list <- lapply(seurat.pseudo.list, function(x) {
            ScaleData(x) %>% RunPCA(verbose=F)
        })
        pseudo.con <- Conos$new(seurat.pseudo.list)
        pseudo.con$buildGraph(k=15, k.self=5, space="PCA", ncomps=30, n.odgenes=2000, matching.method='mNN',
                              metric = 'angular', score.component.variance=T, verbose=T)
        pseudo.con$findCommunities()
        pseudo.con$embedGraph()
        integrated = as.Seurat(pseudo.con)
    }
    integrated
}

int.obj = integration(human.obj, fish.obj, 'h_f')

library(ggplot2)

int.obj@meta.data$orig.ident[int.obj@meta.data$orig.ident=='Library2'] = 'Tumor21'
int.obj@meta.data$orig.ident[int.obj@meta.data$orig.ident=='Library1'] = 'Tumor21'
int.obj@meta.data$orig.ident[int.obj@meta.data$orig.ident=='SeuratProject'] = 'MAST139'

pdf('Fish21_Human139.pdf', width=8, height=6.5)
pd=DimPlot(int.obj, reduction='umap', group.by='seurat_clusters',
           split.by = 'orig.ident', label=T)+ theme(legend.position='bottom')
print(pd)
dev.off()
