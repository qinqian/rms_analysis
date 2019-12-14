#!/bin/bash -ex

#conda activate alvin_py36_graph
for gene in *_seu_markers_tumoronly_res0.8.xls_*cluster_genes.symbols; do
    bsub -J genewalk -n 1 -q rerunnable -M 8000 "genewalk --project $gene --genes $gene --id_type hgnc_symbol --nproc 1"
    #bsub -J lisa -n 6 -q big-multi -M 12000 "bash qsub_lisa.sh $gene"
    #break
done

