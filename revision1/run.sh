# source env.sh
# Rscript marker_view.r

# source env2.sh
# Rscript annotate_celltypes_updated.r --species human --seuratobj /PHShome/qq06/langenau/projects/01_sc_rms/phaseA_explore_rms/20082_recluster2_tumor_only.rds --label 20082 --de /PHShome/qq06/langenau/projects/01_sc_rms/results/seurat_sara/20082_recluster2_tumoronly_res0.8.xls &

#Rscript annotate_celltypes_updated.r --species human --seuratobj /PHShome/qq06/langenau/projects/01_sc_rms/figures/20696_hg19_tumoronly_res0.8_umap.rds --label 20696 --de /PHShome/qq06/langenau/projects/01_sc_rms/results/seurat_sara/20696_hg19_tumoronly_res0.8.xls &
#
#Rscript annotate_celltypes_updated.r --species human --seuratobj /PHShome/qq06/langenau/projects/01_sc_rms/figures/21202_hg19_premrna_tumoronly_res0.8_umap.rds --label 21202 --de /PHShome/qq06/langenau/projects/01_sc_rms/results/seurat_sara/21202_hg19_premrna_tumoronly_res0.8.xls  &
#
#Rscript annotate_celltypes_updated.r --species human --seuratobj /PHShome/qq06/langenau/projects/01_sc_rms/figures/29806_hg19_premrna_tumoronly_res0.8_umap.rds --label 29806 --de /PHShome/qq06/langenau/projects/01_sc_rms/results/seurat_sara/29806_hg19_premrna_tumoronly_res0.8.xls  &

# primary tumor statistics per cell 
# for i in /PHShome/qq06/langenau/projects/01_sc_rms/results/doublets/2*doubl*csv; do echo $i; wc -l $i; grep -c "True" $i; done
# Rscript primary_tumors_statistics.r

###original resolution 0.8 cannot demultiplex some clusters into single cell state
###Rscript annotate_celltypes_updated.R --label MSK72117 --species human --seuratobj /data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds

# Increase resolution
#source env.sh
#Rscript Fig3.r

source env2.sh
Rscript annotate_celltypes_updated.r --species human --seuratobj /data/langenau/alvin_singlecell/projects/01_sc_rms/results/seurat_sara/C12SC2_seurat-object.rds --label C12SC2 --de /PHShome/qq06/langenau/projects/01_sc_rms/results/seurat_sara/C12SC2_SCT_res0.8.xls &

# Rscript refine_resolutions.r
Rscript annotate_celltypes_updated.r --label MSK72117_res0.9 --species human --seuratobj /data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds --de MSK72117_res0.9.csv &

### normal markers GSEA with MAST139/MAST111 DE ranks 
# Rscript normal_development.r # overlap is few

#source env.sh
#for state in EMT Muscle Prolif; do
#    Rscript tumor_aggregate_marker.r --label MAST139 --state $state --seuratobj /PHShome/qq06/langenau/projects/01_sc_rms/data/seurat_obj/20190624_seurat-object_MAST139.rds &
#    Rscript tumor_aggregate_marker.r --label MSK74711 --state $state --seuratobj ~/langenau/projects/01_sc_rms/results/seurat_sara/20191031_MSK74711_seurat-object.rds &
#done

#Rscript normal_marker_enrichment.r

