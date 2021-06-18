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
conda activate deseq



#mkdir -p results/
#for cond in Regular Differentiated; do
#    python identify_doublet.py -mat ${cond}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${cond}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -name ${cond} -gzip &
#done
#wait

# for cond in Regular Differentiated; do
#     Rscript seurat_sara_pipeline.R --seuratobj "${cond}/outs/filtered_feature_bc_matrix/" --label ${cond} --mixtureobj NA --doublet results/${cond}_doublet.csv &
# done
# wait

# for cond in Regular Differentiated; do
#     Rscript seurat_sara_pipeline.R --seuratobj "${cond}/outs/filtered_feature_bc_matrix/" --label ${cond} --mixtureobj NA --doublet results/${cond}_doublet.csv &
# done
# wait

# for cond in Regular Differentiated; do
#     Rscript generate_v4_degenes.R --seuratobj1 results/seurat_sara/${cond}_seurat-object.rds --label ${cond} --res 0.8 &
# done
# wait

for cond in Regular Differentiated; do
    Rscript annotate_celltypes_updated.r --species human --seuratobj results/seurat_sara/${cond}_seurat-object.rds --label ${cond} --de results/seurat_sara/${cond}_SCT_res0.8.xls &
done

