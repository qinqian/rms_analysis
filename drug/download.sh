#!/bin/bash -ex

#| 210427_NB500929_0707_AHGC5VBGXJ | 205976778 | RH41 R S       | Complete |
#| 210527_NB500929_0715_AHHNG2BGXJ | 207586406 | Qiqi 39 S R    | Complete |
#| 210611_NB500929_0721_AHJYGKBGXJ | 208512407 | Qiqi0611       | Complete |
#mkdir MAST39
#cd MAST39
#bs download run --id 207586406 --output .
#cd ../
#mkdir MAST139
#cd MAST139
#bs download run --id 208512407 --output .

__conda_setup="$('/data/pinello/SHARED_SOFTWARE/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/pinello/SHARED_SOFTWARE/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/data/pinello/SHARED_SOFTWARE/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data/pinello/SHARED_SOFTWARE/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
conda activate sc-tutorial

export PATH=/PHShome/qq06/langenau/projects/01_sc_rms/phaseA_explore_rms/cellranger-3.1.0/:${PATH}

# cd MAST39
#echo 'Lane,Sample,Index' > MAST-39_drug.csv
#echo '1-4,MAST-39-Sensitive,SI-GA-A5' >> MAST-39_drug.csv
#echo '1-4,MAST-39-Resistant,SI-GA-A6' >> MAST-39_drug.csv
#cellranger mkfastq --id=MAST-39_drug \
#                   --localcores=20 \
#                   --localmem=64 \
#                   --run=/data/langenau/alvin_singlecell/qiqi/PDX/MAST39 \
#                   --csv=MAST-39_drug.csv
# cd ../

#cd MAST139
#echo 'Lane,Sample,Index' > MAST-139_drug.csv
#echo '1-4,MAST-139-Sensitive,SI-GA-A11' >> MAST-139_drug.csv
#echo '1-4,MAST-139-Resistant,SI-GA-A12' >> MAST-139_drug.csv
#cellranger mkfastq --id=MAST-139_drug \
#                   --localcores=20 \
#                   --localmem=64 \
#                   --run=/data/langenau/alvin_singlecell/qiqi/PDX/MAST139/ \
#                   --csv=MAST-139_drug.csv

#bash cellranger.sh /data/langenau/alvin_singlecell/qiqi/PDX/MAST39/MAST-39_drug/outs/fastq_path/HHNG2BGXJ/ MAST-39-Resistant MAST-39-Resistant /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/
#bash cellranger.sh /data/langenau/alvin_singlecell/qiqi/PDX/MAST39/MAST-39_drug/outs/fastq_path/HHNG2BGXJ/ MAST-39-Sensitive MAST-39-Sensitive /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/
#bash cellranger.sh /data/langenau/alvin_singlecell/qiqi/PDX/MAST139/MAST-139_drug/outs/fastq_path/HJYGKBGXJ/ MAST-139-Resistant MAST-139-Resistant /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/
#bash cellranger.sh /data/langenau/alvin_singlecell/qiqi/PDX/MAST139/MAST-139_drug/outs/fastq_path/HJYGKBGXJ/ MAST-139-Sensitive MAST-139-Sensitive /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/

