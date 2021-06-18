#!/bin/bash -ex

__conda_setup="$('/PHShome/qq06/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/PHShome/qq06/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/PHShome/qq06/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/PHShome/qq06/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
conda activate /data/pinello/SHARED_SOFTWARE/anaconda3/envs/sc-tutorial

# Rscript seurat_monocle3_pipeline.r --seuratobj /data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST139.rds --id MAST139

seurat_sara=(`ls ../data/seurat_obj/*seurat*rds ../results/seurat_sara/20191031_MSK74711_seurat-object.rds ../results/seurat_sara/MAST118_seurat-object.rds ../results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds ../results/seurat_sara/MAST139_1cells_seurat-object.rds ../results/seurat_sara/C12SC2_seurat-object.rds`)

labels=(MAST111 MAST139 MAST35 MAST39 MAST85 MAST95 MSK82489 RH74 MAST85-1 RH74-10 MSK72117 MSK74711 MSK72117-1 MAST118 MAST139-1)

for i in `seq 0 $((${#seurat_sara[@]}-1))`; do
    Rscript seurat_monocle3_pipeline.r --seuratobj ${seurat_sara[$i]} --id ${labels[$i]} &
done

# Rscript seurat_monocle3_pipeline.r --seuratobj ../data/seurat_obj/20190624_seurat-object_MAST95.rds --id MAST95
