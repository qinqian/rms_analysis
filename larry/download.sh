#| 210409_NB500929_0701_AH55N3BGXJ | 204930753 | RD Larry | Complete |
# bs download runs --id 204930753 --output .

echo 'Lane,Sample,Index' > RD_larry.csv
echo '1-4,Regular,SI-GA-F1' >> RD_larry.csv
echo '1-4,Differentiated,SI-GA-F2' >> RD_larry.csv

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

#/data/molpath/software/10x/cellranger-3.0.2/cellranger mkfastq --id=RD_larry \
#                   --localcores=30 \
#                   --localmem=64 \
#                   --run=/data/langenau/alvin_singlecell/larry/claudia \
#                   --csv=RD_larry.csv

#bash cellranger.sh /PHShome/qq06/langenau/larry/claudia/RD_larry/outs/fastq_path/H55N3BGXJ/ Differentiated Differentiated /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/
#bash cellranger.sh /PHShome/qq06/langenau/larry/claudia/RD_larry/outs/fastq_path/H55N3BGXJ/ Regular Regular /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/
