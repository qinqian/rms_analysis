

#conda activate /PHShome/qq06/miniconda3/envs/lisa2

mkdir -p lisa_data/
#for i in *markers.xls; do
#    sed 1d $i | cut -f 1 > ${i}.genes
#done

lisa multi hg38 *markers*xls.genes --rp_map enhanced_10K -o lisa_data/ --save_metadata

