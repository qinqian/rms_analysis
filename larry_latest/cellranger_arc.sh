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

export PATH=/PHShome/qq06/ssd/alvin/larry/cellranger-arc-2.0.2/:${PATH}
cellranger-arc count --id=${output} \
                 --transcriptome=${genome} \
                 --fastqs=${input} \
                 --sample=${prefix} \
                 --localcores=24 \
                 --localmem=64

