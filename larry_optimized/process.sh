#!/bin/bash -ex

#__conda_setup="$('/PHShome/qq06/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
#if [ $? -eq 0 ]; then
#    eval "$__conda_setup"
#else
#    if [ -f "/PHShome/qq06/miniconda3/etc/profile.d/conda.sh" ]; then
#        . "/PHShome/qq06/miniconda3/etc/profile.d/conda.sh"
#    else
#        export PATH="/PHShome/qq06/miniconda3/bin:$PATH"
#    fi
#fi
#unset __conda_setup
#conda activate deseq

#mkdir -p results/
#for cond in Regular_1 Regular_2 Diff_3 Diff_4; do
#    python identify_doublet.py -mat ${cond}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${cond}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -name ${cond} -gzip &
#done
#wait

#for cond in Regular_1 Regular_2 Diff_3 Diff_4; do
#     Rscript seurat_sara_pipeline.R --seuratobj "${cond}/outs/filtered_feature_bc_matrix/" --label ${cond} --mixtureobj NA --doublet results/${cond}_doublet.csv &
#done
#wait

#for cond in Regular_1 Regular_2 Diff_3 Diff_4; do
#    #Rscript generate_v4_degenes.R --seuratobj1 results/seurat_sara/${cond}_seurat-object.rds --label ${cond} --res 0.4 &
#    Rscript generate_v4_degenes.R --seuratobj1 results/seurat_sara/${cond}_seurat-object.rds --label ${cond} --res 1.0 &
#done
#wait

#for cond in Regular_1 Regular_2 Diff_3 Diff_4; do
#    #Rscript annotate_celltypes_updated.r --species human --seuratobj results/seurat_sara/${cond}_seurat-object.rds --label ${cond} --de results/seurat_sara/${cond}_SCT_res0.8.xls &
#    #Rscript annotate_celltypes_updated.r --species human --seuratobj results/seurat_sara/${cond}_seurat-object.rds --label ${cond} --de results/seurat_sara/${cond}_SCT_res0.4.xls --res 0.4 &
#    Rscript annotate_celltypes_updated.r --species human --seuratobj results/seurat_sara/${cond}_seurat-object.rds --label ${cond} --de results/seurat_sara/${cond}_SCT_res1.xls --res 1.0 &
#done
#

# for integrative clustering analysis
# Rscript annotate_celltypes_updated.r --species human --seuratobj second_Larry_integrative_res1.rds --label integrated --de integrative_clusters_markers.txt --res 1.0

#conda activate /data/pinello/SHARED_SOFTWARE/anaconda3/envs/sc-tutorial
for i in Regular_1 Regular_2 Diff_3 Diff_4; do
   # velocyto run10x ${i} /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/genes/genes.gtf &
   # python velocity_pipeline_dynamical_latentime.py -l ../results/seurat_intersect_velocity/${i}_vel.loom -s ../results/seurat_intersect_velocity/${i}_seu.rds -n MSK72117-1
   break
done

