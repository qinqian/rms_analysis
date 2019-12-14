#!/bin/bash
#BSUB -J test
#BSUB -o output/test-%J.out
#BSUB -e output/test-%J.err
#BSUB -q big
#BSUB -n 6
#BSUB -M 64000
#BSUB -R rusage[mem=64000]

seurat=${1}
vel=${2}
output=${3}
cluster=${4}

cd /data/langenau/alvin_singlecell/01_rms_projects/01_fish/src
/data/pinello/SHARED_SOFTWARE/anaconda3/condabin/conda activate sc-tutorial
Rscript velocity_pipeline.R --seuratobj $seurat --velobj $vel --label $output --clusterlabel $cluster 2>&1 >${output}_velR.log

