### | 220503_VH00656_56_AAALWH2HV     | 235225992 | RD_Larry              | Complete |
## mkdir -p latest
## bs download runs --id 235225992 --output latest
##
## New design with Dual Index
#echo 'Lane,Sample,Index' > RD_larry.csv
#echo '1-4,day0_seq1,SI-TT-A1' >> RD_larry.csv
#echo '1-4,day6_seq2,SI-TT-A2' >> RD_larry.csv
#echo '1-4,day10_seq3,SI-TT-A3' >> RD_larry.csv
#echo '1-4,day0_seq4,SI-TT-A8' >> RD_larry.csv
#echo '1-4,day4_seq5,SI-TT-A9' >> RD_larry.csv
#echo '1-4,day4_seq6,SI-TT-A10' >> RD_larry.csv
#echo '1-4,day4_seq7,SI-TT-A11' >> RD_larry.csv
#
#__conda_setup="$('/data/pinello/SHARED_SOFTWARE/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
#if [ $? -eq 0 ]; then
#    eval "$__conda_setup"
#else
#    if [ -f "/data/pinello/SHARED_SOFTWARE/anaconda3/etc/profile.d/conda.sh" ]; then
#        . "/data/pinello/SHARED_SOFTWARE/anaconda3/etc/profile.d/conda.sh"
#    else
#        export PATH="/data/pinello/SHARED_SOFTWARE/anaconda3/bin:$PATH"
#    fi
#fi
#unset __conda_setup
#conda activate sc-tutorial
#
#export PATH=/PHShome/qq06/ssd/alvin/larry/cellranger-7.0.1/:${PATH}
#
#### https://kb.10xgenomics.com/hc/en-us/articles/115003082371-How-to-demultiplex-a-single-indexed-library-on-a-dual-indexed-flow-cell-
#### https://kb.10xgenomics.com/hc/en-us/articles/360052777572-How-to-use-masking-parameter-while-demultiplexing-10x-sequencing-data-
#### https://www.biostars.org/p/344768/
#### > version 5.0
## cellranger mkfastq --id=RD_larry \
##                    --localcores=60 \
##                    --localmem=128 \
##                    --use-bases-mask=Y28n*,I10,I10,Y90n* \
##                    --run=/srv/local/alvin/larry/latest \
##                    --csv=RD_larry.csv
##                    #--filter-dual-index 
#
#### < version 3.0, [error] Formatting error in sample sheet: /PHShome/qq06/ssd/alvin/larry/RD_larry.csv
##cellranger mkfastq --id=RD_larry --localcores=80 --localmem=128 --csv=RD_larry.csv --ignore-dual-index --use-bases-mask=Y28n*,I10,I10,Y90n* --run=/srv/local/alvin/larry/latest
#
for day in day0_seq1 day6_seq2 day10_seq3; do
    #nice -n 5 bash cellranger.sh /PHShome/qq06/ssd/alvin/larry/RD_larry/outs/fastq_path/AAALWH2HV/ $day $day refdata-gex-GRCh38-and-mm10-2020-A/ &
    echo nice -n 5 bash cellranger_expect_cells.sh /PHShome/qq06/ssd/alvin/larry/RD_larry/outs/fastq_path/AAALWH2HV/ $day ${day}_GRCh38 refdata-gex-GRCh38-2020-A/ & # slightly higher proportions

    #nice -n 5 bash cellranger_spliced.sh /PHShome/qq06/ssd/alvin/larry/RD_larry/outs/fastq_path/AAALWH2HV/ $day ${day}_GRCh38_spliced refdata-gex-GRCh38-2020-A/ &
done
wait

# Suppose multiome data
# export PATH=/PHShome/qq06/ssd/alvin/larry/cellranger-arc-2.0.2/:${PATH}
#cellranger-arc mkfastq --id=RD_larry_multiome \
#                       --localcores=60 \
#                       --localmem=128 \
#                       --use-bases-mask=Y28n*,I10,I10,Y90n* \
#                       --run=/srv/local/alvin/larry/latest \
#                       --csv=RD_larry.csv

# https://kb.10xgenomics.com/hc/en-us/articles/360059656912-Can-I-analyze-only-the-Gene-Expression-data-from-my-single-cell-multiome-experiment-
#for day in day0_seq4 day4_seq5 day4_seq6 day4_seq7; do
    #nice -n 5 cellranger count --id=${day}_arc --transcriptome=/srv/local/alvin/larry/refdata-gex-GRCh38-2020-A/ --fastqs=/PHShome/qq06/ssd/alvin/larry/RD_larry/outs/fastq_path/AAALWH2HV/ --sample=${day} --chemistry=ARC-v1 --localcores=36 --localmem=144 &
    #nice -n 5 cellranger count --id=${day}_arc_mix_test2 --transcriptome=/srv/local/alvin/larry/refdata-gex-GRCh38-and-mm10-2020-A/ --fastqs=/PHShome/qq06/ssd/alvin/larry/RD_larry/outs/fastq_path/AAALWH2HV/ --sample=${day} --chemistry=ARC-v1 --localcores=36 --localmem=144 &
#done

#tar cvfz all_cellrangers.tar.gz */outs/web_summary.html */outs/analysis/gem* 

