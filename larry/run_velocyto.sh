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
#conda activate /data/pinello/SHARED_SOFTWARE/anaconda3/envs/sc-tutorial

#velocyto run10x Differentiated/ /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/genes/genes.gtf
#velocyto run10x Regular/ /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/genes/genes.gtf

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

#Rscript seurat_pipeline.R --seuratobj Differentiated/velocyto/Differentiated.loom --label Differentiated --finalres 0.8 --tumor -1 --assaytype spliced --species human
#Rscript seurat_pipeline.R --seuratobj Regular/velocyto/Regular.loom --label Regular --finalres 0.8 --tumor -1 --assaytype spliced --species human

#Rscript intersect_seurat_velocity_toloom.R --seuratobj results/seurat_sara/Regular_seurat-object.rds --velobj Regular_seurat_obj_tumors.rds --label Regular --species human
#Rscript intersect_seurat_velocity_toloom.R --seuratobj results/seurat_sara/Differentiated_seurat-object.rds --velobj Differentiated_seurat_obj_tumors.rds --label Differentiated --species human

mkdir ../results/velocity_dynamical/
#python velocity_pipeline_dynamical_latentime.py -l Differentiated/velocyto/Differentiated.loom -s ../results/seurat_intersect_velocity/Differentiated_seu.rds -n Differentiated

python velocity_pipeline_dynamical_latentime.py -l Regular/velocyto/Regular.loom -s ../results/seurat_intersect_velocity/Regular_seu.rds -n Regular
