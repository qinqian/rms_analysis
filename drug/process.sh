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

# mkdir -p results
# for sample in *-*-*; do
#    echo $sample
#    label=$(basename $sample)
#    if [ ! -s results/${label}_doublet.csv ] || [ ! -s results/${label}_scrub_umap.pdf ]; then
#        python identify_doublet.py -mat ${sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -name $label -gzip
#    fi
# done
# wait

# for sample in *-*-*; do
#    echo $sample
#    label=$(basename $sample)
#    input_sample=${sample}/outs/filtered_feature_bc_matrix
#    gem=${sample}/outs/analysis/
#    Rscript seurat_sara_pipeline.R --seuratobj ${input_sample} --label $label  --mixtureobj ${gem}/gem_classification.csv --doublet results/${label}_doublet.csv &
#    # break
# done

# for sample in results/seurat_sara/*_seurat-object.rds; do
#    echo $sample
#    label=$(basename $sample)
#    label=${label/_seurat-object.rds/}
#    # Rscript generate_v4_degenes.R --seuratobj1 $sample --label $label  &
#    Rscript generate_v4_degenes.R --seuratobj1 $sample --label $label 
# done

for sample in results/seurat_sara/*_seurat-object.rds; do
    echo $sample
    label=$(basename $sample)
    label=${label/_seurat-object.rds/}
    Rscript annotate_celltypes_updated.r --species human --seuratobj $sample --label $label --de results/seurat_sara/${label}_SCT_res0.8.xls --res 0.8 &
    # break
done
