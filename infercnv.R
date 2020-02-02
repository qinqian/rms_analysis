library(infercnv)               ## from broad
library(HoneyBADGER)            ## from Jean Fan
library(parallel)
library(GenomicRanges)
library(Rsamtools)

## infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
##                                     annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
##                                     delim="\t",
##                                     gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
##                                     ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 

## infercnv_obj = infercnv::run(infercnv_obj,
##                              cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
##                              out_dir=tempfile(), 
##                              cluster_by_groups=TRUE, 
##                              denoise=TRUE,
##                              HMM=TRUE)

bam = '/data/langenau/human_rms_pdxs/20190418_MAST111_5Kcells_hg19/outs/possorted_genome_bam.bam'
index = '/data/langenau/human_rms_pdxs/20190418_MAST111_5Kcells_hg19/outs/possorted_genome_bam.bam.bai'

results_all = list()
for (i in paste0('chr', 1:22)) {
    load(system.file("ExAC", paste0("ExAC_", i, ".RData"), package = "HoneyBADGER"))
    vi <- seqlevels(snps) %in% as.character(c(1:22))
    table(vi)
    seqlevels(snps) <- seqlevels(snps)[vi]
    results_all[[i]] <- getSnpMats(snps, bam, index)
}

system('mkdir -p ../results/cnv')
saveRDS(results_all, '../results/cnv/MAST111_honeybadger.rds')

results_all = readRDS("../results/cnv/MAST111_honeybadger.rds")

mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl') ## current version

seur = readRDS('/data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds')

hb <- new('HoneyBADGER', name='MAST111')

hb$setGexpMats(as.matrix(seur$RNA@scale.data), ref, mart.obj, filter=FALSE, scale=FALSE, verbose=TRUE)



library(HoneyBADGER)
require(biomaRt) ## for gene coordinates
#mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host = "jul2015.archive.ensembl.org") ## version used in manuscript
mart.obj <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl') ## current version
saveRDS(mart.obj, file='mart_homo.rds')

hb <- new('HoneyBADGER', name='MGH31')
hb$setGexpMats(gexp, ref, mart.obj, filter=FALSE, scale=FALSE, verbose=TRUE)
