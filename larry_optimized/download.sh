#| 210531_NB500929_0716_AHHNK3BGXJ | 207805620 | RD Larry par reg   | Complete |
#| 210602_NB500929_0717_AHHTWVBGXJ | 207907736 | RD Larry diff regu | Complete |
#mkdir reg
#cd reg
#bs download runs --id 207805620 --output .
#
#cd ..
#
#mkdir diff
#cd diff
#bs download runs --id 207907736 --output .

#echo 'Lane,Sample,Index' > RD_larry_2ndrun.csv
#echo '1-4,Regular_1,SI-GA-A7' >> RD_larry_2ndrun.csv
#echo '1-4,Regular_2,SI-GA-A8' >> RD_larry_2ndrun.csv
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
#/data/molpath/software/10x/cellranger-3.0.2/cellranger mkfastq --id=RD_larry_2ndrun \
#                   --localcores=30 \
#                   --localmem=64 \
#                   --run=/data/langenau/alvin_singlecell/larry/claudia/second_run/reg/ \
#                   --csv=RD_larry_2ndrun.csv

#bash cellranger.sh /PHShome/qq06/langenau/larry/claudia/second_run/RD_larry_2ndrun/outs/fastq_path/HHNK3BGXJ/ Regular_1 Regular_1 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/ 
#bash cellranger.sh /PHShome/qq06/langenau/larry/claudia/second_run/RD_larry_2ndrun/outs/fastq_path/HHNK3BGXJ/ Regular_2 Regular_2 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/ 

#echo 'Lane,Sample,Index' > RD_larry_3ndrun.csv
#echo '1-4,Diff_3,SI-GA-A9' >> RD_larry_3ndrun.csv
#echo '1-4,Diff_4,SI-GA-A10' >> RD_larry_3ndrun.csv
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
#/data/molpath/software/10x/cellranger-3.0.2/cellranger mkfastq --id=RD_larry_3ndrun \
#                   --localcores=30 \
#                   --localmem=64 \
#                   --run=/data/langenau/alvin_singlecell/larry/claudia/second_run/diff/ \
#                   --csv=RD_larry_3ndrun.csv

#bash cellranger.sh /PHShome/qq06/langenau/larry/claudia/second_run/RD_larry_3ndrun/outs/fastq_path/HHTWVBGXJ/ Diff_3 Diff_3 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/
bash cellranger.sh /PHShome/qq06/langenau/larry/claudia/second_run/RD_larry_3ndrun/outs/fastq_path/HHTWVBGXJ/ Diff_4 Diff_4 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/

