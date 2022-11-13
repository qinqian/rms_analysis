#!/bin/bash
#BSUB -J test
#BSUB -o output/test-%J.out
#BSUB -e output/test-%J.err
#BSUB -q big
#BSUB -n 6
#BSUB -M 64000
#BSUB -R rusage[mem=64000]

input=${1}
prefix=${2}
output=${3}
genome=${4}

#module load cellranger/3.0.2 # when using bsub
export PATH=/PHShome/qq06/langenau/projects/01_sc_rms/phaseA_explore_rms/cellranger-3.1.0/:${PATH}
cellranger count --id=${output} \
                 --transcriptome=${genome} \
                 --fastqs=${input} \
                 --sample=${prefix} \
                 --localcores=24 \
                 --localmem=128
