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

mkdir -p results
##for lib in day0_seq1_GRCh38_expectcells10k day0_seq1_GRCh38_reanalyze day6_seq2_GRCh38 day10_seq3_GRCh38 day0_seq4_arc day4_seq5_arc day4_seq6_arc day4_seq7_arc; do
##    #python identify_doublet.py -mat ${lib}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${lib}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -name ${lib} -gzip &
##    echo python identify_doublet.py -mat ${lib}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${lib}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -name ${lib} -gzip 
##done

for lib in day0_seq1_GRCh38_expectcells10k day0_seq1_GRCh38_reanalyze day6_seq2_GRCh38 day10_seq3_GRCh38 day0_seq4_arc day4_seq5_arc day4_seq6_arc day4_seq7_arc; do
     Rscript seurat_sara_pipeline.R --seuratobj "${lib}/outs/filtered_feature_bc_matrix/" --label ${lib} --mixtureobj NA --doublet results/${lib}_doublet.csv &
done
wait

for lib in day0_seq1_GRCh38_expectcells10k day0_seq1_GRCh38_reanalyze day6_seq2_GRCh38 day10_seq3_GRCh38 day0_seq4_arc day4_seq5_arc day4_seq6_arc day4_seq7_arc; do
    Rscript generate_v4_degenes.R --seuratobj1 results/seurat_sara/${lib}_seurat-object.rds --label ${lib} --res 0.8 &
done
wait

for lib in day0_seq1_GRCh38_expectcells10k day0_seq1_GRCh38_reanalyze day6_seq2_GRCh38 day10_seq3_GRCh38 day0_seq4_arc day4_seq5_arc day4_seq6_arc day4_seq7_arc; do
    Rscript annotate_celltypes_updated.r --species human --seuratobj results/seurat_sara/${lib}_seurat-object.rds --label ${lib} --de results/seurat_sara/${lib}_SCT_res0.8.xls &
done

##conda activate /data/pinello/SHARED_SOFTWARE/anaconda3/envs/sc-tutorial
#for i in Regular_1 Regular_2 Diff_3 Diff_4; do
#   # velocyto run10x ${i} /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/genes/genes.gtf &
#   # python velocity_pipeline_dynamical_latentime.py -l ../results/seurat_intersect_velocity/${i}_vel.loom -s ../results/seurat_intersect_velocity/${i}_seu.rds -n MSK72117-1
#   break
#done


