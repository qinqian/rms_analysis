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


src=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/src
data=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data
tmp=/PHShome/qq06/scratch/01_sc_rms
genome=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome
result=../results/seurat_sara
flurescent=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/flurescent_colors/

mkdir -p $data
mkdir -p $tmp
mkdir -p $src
mkdir -p $genome
mkdir -p $result
export PATH=${PATH}:/PHShome/qq06/projects/01_sc_rms/phaseA_explore_rms/cellranger-3.1.0/

download_demo() {
    cd $data
    #wget -c https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
    tar xvfz pbmc3k_filtered_gene_bc_matrices.tar.gz
}

build_genome() {
  cd $genome
  #wget -c ftp://ftp.ensembl.org/pub/release-97/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
  #wget -c ftp://ftp.ensembl.org/pub/release-97/gtf/danio_rerio/Danio_rerio.GRCz11.97.gtf.gz
  #wget -c ftp://ftp.ensembl.org/pub/release-97/gtf/danio_rerio/Danio_rerio.GRCz11.97.chr.gtf.gz

  if [ -s Danio_rerio.GRCz11.97.gtf.gz ]; then
    gunzip Danio_rerio.GRCz11.97.gtf.gz
  fi

  if [ -s Danio_rerio.GRCz11.dna.primary_assembly.fa.gz ]; then
    gunzip Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
  fi

  if [ ! -s Danio_rerio.GRCz11.97.filtered.gtf ]; then

      #25	ensembl_havana	gene	12841940	12844002	.	+	.	gene_id "ENSDARG00000099401"; gene_version "2"; gene_name "ccl33.3"; gene_source "ensembl_havana"; gene_biotype "processed_transcript";
      #ccl3.3 has been excluded by cellranger filter
      cellranger mkgtf Danio_rerio.GRCz11.97.gtf Danio_rerio.GRCz11.97.filtered.gtf \
                       --attribute=gene_biotype:protein_coding \
                       --attribute=gene_biotype:lincRNA \
                       --attribute=gene_biotype:antisense \
                       --attribute=gene_biotype:IG_LV_gene \
                       --attribute=gene_biotype:IG_V_gene \
                       --attribute=gene_biotype:IG_V_pseudogene \
                       --attribute=gene_biotype:IG_D_gene \
                       --attribute=gene_biotype:IG_J_gene \
                       --attribute=gene_biotype:IG_J_pseudogene \
                       --attribute=gene_biotype:IG_C_gene \
                       --attribute=gene_biotype:IG_C_pseudogene \
                       --attribute=gene_biotype:TR_V_gene \
                       --attribute=gene_biotype:TR_V_pseudogene \
                       --attribute=gene_biotype:TR_D_gene \
                       --attribute=gene_biotype:TR_J_gene \
                       --attribute=gene_biotype:TR_J_pseudogene \ #--attribute=gene_biotype:TR_C_gene

      cellranger mkref --genome=GRCz11 --fasta=Danio_rerio.GRCz11.dna.primary_assembly.fa --genes=Danio_rerio.GRCz11.97.filtered.gtf
  fi

  #awk 'BEGIN {FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' Danio_rerio.GRCz11.97.filtered.gtf > Danio_rerio.GRCz11.97.filtered.premrna.gtf
  #cellranger mkref --genome=GRCz11_premrna \
  #                   --fasta=Danio_rerio.GRCz11.dna.primary_assembly.fa \
  #                   --genes=Danio_rerio.GRCz11.97.filtered.premrna.gtf

  #build cellranger reference with transgenic ORF and color sequences
  #colors.fa  GRCz10_plus_myc_colors.gtf  test  yan_ally_madeline_orf.fa  yan_ally_madeline_orf.gtf
  #python ${src}/generate_manual_gtf.py ${flurescent}/yan_ally_madeline_orf.fa

  #rm -f ${flurescent}/colors.gtf
  #for i in `grep '>' ${flurescent}/colors.fa`; do
  #    grep $(echo $i | sed -e 's/>//g') ${flurescent}/GRCz10_plus_myc_colors.gtf >> ${flurescent}/colors.gtf
  #done

  #cat Danio_rerio.GRCz11.dna.primary_assembly.fa ${flurescent}/yan_ally_madeline_orf.fa ${flurescent}/colors.fa > Danio_rerio.GRCz11.dna.primary_assembly_with_color_and_orf_seq.fa
  #cat Danio_rerio.GRCz11.97.filtered.gtf ${flurescent}/yan_ally_madeline_orf.gtf ${flurescent}/colors.gtf > Danio_rerio.GRCz11.97.filtered.with_color_and_orf.gtf
  #awk 'BEGIN {FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' Danio_rerio.GRCz11.97.filtered.with_color_and_orf.gtf > Danio_rerio.GRCz11.97.filtered.with_color_and_orf.premrna.gtf
  #
  #cellranger mkref --genome=GRCz11_with_color_and_orf.v1 --fasta=Danio_rerio.GRCz11.dna.primary_assembly_with_color_and_orf_seq.fa --genes=Danio_rerio.GRCz11.97.filtered.with_color_and_orf.gtf

  #cellranger mkref --genome=GRCz11_with_color_and_orf.v1.premrna \
  #                 --fasta=Danio_rerio.GRCz11.dna.primary_assembly_with_color_and_orf_seq.fa \
  #                 --genes=Danio_rerio.GRCz11.97.filtered.with_color_and_orf.premrna.gtf

  python ${src}/generate_manual_gtf.py ${flurescent}/yan_ally_madeline_orf.fa
  cat Danio_rerio.GRCz11.97.gtf ${flurescent}/yan_ally_madeline_orf.gtf ${flurescent}/colors.gtf > Danio_rerio.GRCz11.97.colored_without_filter_v2.gtf
  cat Danio_rerio.GRCz11.dna.primary_assembly.fa ${flurescent}/yan_ally_madeline_orf.fa ${flurescent}/colors.fa>Danio_rerio.GRCz11.dna.primary_assembly_with_color_orf_seq_v2.fa
  cellranger mkref --genome=GRCz11_with_color_and_orf.withoutfilter.v2 --fasta=Danio_rerio.GRCz11.dna.primary_assembly_with_color_orf_seq_v2.fa --genes=Danio_rerio.GRCz11.97.colored_without_filter_v2.gtf
}

submit_cellranger() {
  cd ~/scratch/
  #bsub results in ~/scratch disappear...
  #bsub -J zebrafish1 -n 6 -q normal -M 8000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor21sort/ Tumor21sort tumor21_zebrafish"

  # bcl to fastq
  #echo 'Lane,Sample,Index' > Tumor24.csv
  #echo '1-4,Tumor24,SI-GA-C5' >> Tumor24.csv
  #cellranger mkfastq --id=Tumor24 \
  #                   --localcores=8 \
  #                   --localmem=64 \
  #                   --run=/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/bcl/Tumor24 \
  #                   --csv=Tumor24.csv

  #bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor24/ Tumor24 Tumor24_zebrafish /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11
  #bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor21sort/ Tumor21sort Tumor21_zebrafish /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11
  #bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor22sort/ Tumor22sort Tumor22_zebrafish /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11
  #bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor22bulk/ Tumor22bulk Tumor22_2ndlibrary_zebrafish /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11
  #bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor21bulk/ Tumor21bulk Tumor21_2ndlibrary_zebrafish /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11

  #bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor24/ Tumor24 Tumor24_zebrafish_with_orf_color_bsub /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1/
  #bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor21sort/ Tumor21sort Tumor21_zebrafish_with_orf_color /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1/
  #bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor22sort/ Tumor22sort Tumor22_zebrafish_with_orf_color /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1
  #bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor21bulk/ Tumor21bulk Tumor21_2ndlibrary_zebrafish_with_orf_color /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1
  #bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor22bulk/ Tumor22bulk Tumor22_2ndlibrary_zebrafish_with_orf_color /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1

  #bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor24/ Tumor24 Tumor24_zebrafish_with_orf_color /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1"
  #bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor21bulk/ Tumor21bulk Tumor21_2ndlibrary_zebrafish_with_orf_color /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1"
  #bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor22bulk/ Tumor22bulk Tumor22_2ndlibrary_zebrafish_with_orf_color /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1"

  #bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor24/ Tumor24 Tumor24_zebrafish_with_orf_color_premrna /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1.premrna/"
  #bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor21bulk/ Tumor21bulk Tumor21_2ndlibrary_zebrafish_with_orf_color_premrna /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1.premrna/"
  #bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor22bulk/ Tumor22bulk Tumor22_2ndlibrary_zebrafish_with_orf_color_premrna /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1.premrna/"
  #bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor21sort/ Tumor21sort Tumor21_zebrafish_with_orf_color_premrna /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1.premrna/"
  #bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor22sort/ Tumor22sort Tumor22_zebrafish_with_orf_color_premrna /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.v1.premrna/"

  bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor24/ Tumor24 Tumor24_zebrafish_with_orf_color_v2 /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.withoutfilter.v2"
  bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor21bulk/ Tumor21bulk Tumor21_2ndlibrary_zebrafish_with_orf_color_v2 /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.withoutfilter.v2"
  bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor22bulk/ Tumor22bulk Tumor22_2ndlibrary_zebrafish_with_orf_color_v2 /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.withoutfilter.v2"
  bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor21sort/ Tumor21sort Tumor21_zebrafish_with_orf_color_v2 /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.withoutfilter.v2"
  bsub -J zebrafish1 -n 8 -q big-multi -M 32000 "bash ${src}/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/fastq/Tumor22sort/ Tumor22sort Tumor22_zebrafish_with_orf_color_v2 /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/GRCz11_with_color_and_orf.withoutfilter.v2"
}

get_tracerseq() {
  #wget -c "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE112294&format=file" -O ${data}/GSE112294.tar
  wget --no-check-certificate -c https://kleintools.hms.harvard.edu/paper_websites/wagner_zebrafish_timecourse2018/WagnerScience2018.h5ad -O ${data}/WagnerScience2018.h5ad
}

get_urd() {
  #wget https://portals.broadinstitute.org/single_cell/data/public/SCP162/single-cell-reconstruction-of-developmental-trajectories-during-zebrafish-embryogenesis?filename=URD_Zebrafish_Object.rds
  #manual download https://portals.broadinstitute.org/single_cell/study/SCP162/single-cell-reconstruction-of-developmental-trajectories-during-zebrafish-embryogenesis#study-download
  echo 'URD'
}

run_velocity() {
  #http://velocyto.org/velocyto.py/tutorial/cli.html#download-genome-annotation-file
  # velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish/ /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/genome/Danio_rerio.GRCz11.97.filtered.gtf
  # velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish/ /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/genome/Danio_rerio.GRCz11.97.filtered.gtf
  # velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish/ /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/genome/Danio_rerio.GRCz11.97.filtered.gtf
  # velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish/ /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/genome/Danio_rerio.GRCz11.97.filtered.gtf
  # velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish/ /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/genome/Danio_rerio.GRCz11.97.filtered.gtf

   #for s in ${data}/cellranger_counts/*orf*; do
   #    echo $s
   #    velocyto run10x $s /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/Danio_rerio.GRCz11.97.filtered.with_color_and_orf.gtf
   #done

   #for s in ${data}/cellranger_counts/*orf*premrna; do
   #    echo $s
   #    velocyto run10x $s /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/Danio_rerio.GRCz11.97.filtered.with_color_and_orf.gtf
   #done

   #for s in ${data}/cellranger_counts/*orf*v2; do
   #    echo $s
   #    velocyto run10x $s /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome/Danio_rerio.GRCz11.97.colored_without_filter_v2.gtf &
   #done

   #for i in `cut -f 1 /PHShome/qq06/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/final_list.txt | grep -v "111" | grep -v "139"`; do
   #    velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/${i} /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/genes/genes.gtf &
   #done
   #wait
   #velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/20190801_MAST85-1cell_5Kcells_hg19/ /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/genes/genes.gtf
   velocyto run10x /data/langenau/human_rms_pdxs/20190617_RH74-10cells_5Kcells_hg19/ /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/genes/genes.gtf
}


identify_doublet() {
  # for vm in /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/*v2; do
  #     echo $vm
  #     #python identify_doublet.py -mat ${vm}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${vm}/outs/filtered_feature_bc_matrix/features.tsv.gz -name $(basename $vm) &
  #     python identify_doublet.py -mat ${vm}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${vm}/outs/filtered_feature_bc_matrix/features.tsv.gz -name $(basename $vm) -gzip &
  # done
    for vm in /PHShome/qq06/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/*; do
	echo $vm
	python identify_doublet.py -mat ${vm}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${vm}/outs/filtered_feature_bc_matrix/features.tsv.gz -name $(basename $vm) -gzip
    done
}

main() {
  #First trial on demo PBMC
  #download_demo

  #Build cellranger genome with additional color coding sequences
  #build_genome

  #other resources colllections
  #get_tracerseq
  #get_urd

  #Run cellranger for zebrafish datasets
  #submit_cellranger

  #Preprocess splice and unspilced transcripts reads by using velocity pipeline
  # run_velocity

  #QC by checking doublet occurence
  #identify_doublet

  #previous key steps
  #1. Seurat normalization and excluding blood cell by hand, output RDS for normalized read counts, barcode, rna for scanpy
  #2. Focus on tumor cells with Seurat velocity wrappers, original script Rscript recluster_velocity_excludingblood.R

  #Now refactor into two pipeline, tumor clusters are determined by eyes, --tumor to control the cluster identity to subset, which is used to exclude non-tumor cells
  #mkdir -p ../results/seurat

  ####################
  ## for zebrafish
  ####################
  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix/' --label Tumor24 --finalres 0.05 --tumor 0 3
  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' --label Tumor21 --finalres 0.1 --tumor 0 1
  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' --label Tumor22 --finalres 0.1 --tumor 0 1 4

  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix/' --label Tumor24_unfilter --finalres 0.05 --tumor -1 &
  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' '/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' --label Tumor21_unfilter --finalres 0.1 --tumor -1 &
  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' '/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' --label Tumor22_unfilter --finalres 0.1 --tumor -1 &

  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color_v2/velocyto/Tumor24_zebrafish_with_orf_color_v2.loom' --label Tumor24_velocity --finalres 0.05 --tumor -1 --assaytype spliced
  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish_with_orf_color_v2/velocyto/Tumor21_zebrafish_with_orf_color_v2.loom' '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2/velocyto/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2.loom' --label Tumor21_velocity --finalres 0.1 --tumor -1 --assaytype spliced
  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish_with_orf_color_v2/velocyto/Tumor22_zebrafish_with_orf_color_v2.loom' '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2/velocyto/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2.loom' --label Tumor22_velocity --finalres 0.1 --tumor -1 --assaytype spliced

  #intersect between seurat and velocity cell barcodes
  #Rscript intersect_seurat_velocity_toloom.R --seuratobj ../results/seurat/Tumor24_seurat_obj_tumors.rds --velobj '../results/seurat/Tumor24_velocity_seurat_obj_tumors.rds' --label Tumor24 #2>&1 >Tumor24_intersect.log
  #Rscript intersect_seurat_velocity_toloom.R --seuratobj ../results/seurat/Tumor21_seurat_obj_tumors.rds --velobj '../results/seurat/Tumor21_velocity_seurat_obj_tumors.rds' --label Tumor21 #2>&1 >Tumor21_intersect.log
  #Rscript intersect_seurat_velocity_toloom.R --seuratobj ../results/seurat/Tumor22_seurat_obj_tumors.rds --velobj '../results/seurat/Tumor22_velocity_seurat_obj_tumors.rds' --label Tumor22 # 2>&1 >Tumor22_intersect.log &

  #for i in ../results/seurat_intersect_velocity/Tumor*seu.rds ; do
  #   echo $i
  #   label=$(basename $i)
  #   # conda activate alvin_resvel
  #   #Rscript generate_v6_degenes.R --seuratobj1 $i --label ${label/.rds}
  #   #conda activate sc-tutorial
  #   Rscript annotate_celltypes.R --seuratobj $i --label ${label/.rds/} --species fish &
  #   #break
  #done

  #n=0
  #for index in ../results/seurat_intersect_velocity/*seu_markers*xls; do
  #    let n++
  #    label=$(basename $index)
  #    if [ $n -gt 10 ]; then
  #        echo Rscript gsea_pipeline.R --seuratde $index --label $label --species fish
  #        Rscript gsea_pipeline.R --seuratde $index --label $label --species fish
  #        break
  #    else
  #        echo 1
  #        #echo Rscript gsea_pipeline.R --seuratde $index --label $label
  #    fi
  #done

  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/Tumor24_seu.rds --velobj ../results/seurat_intersect_velocity/Tumor24_vel.rds --label Tumor24_velocityR --clusterlabel Tumor24 --species fish
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/Tumor21_seu.rds --velobj ../results/seurat_intersect_velocity/Tumor21_vel.rds --label Tumor21_velocityR --clusterlabel Tumor21 --species fish
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/Tumor22_seu.rds --velobj ../results/seurat_intersect_velocity/Tumor22_vel.rds --label Tumor22_velocityR --clusterlabel Tumor22 --species fish

  #labels=(Tumor24 Tumor21 Tumor22)
  #for index in ${!labels[*]}; do
  #    echo python velocity_pipeline.py -l ../results/seurat_intersect_velocity/${labels[$index]}_vel.loom -s ../results/seurat_intersect_velocity/${labels[$index]}_seu.rds -n ${labels[$index]} --species fish
  #    python velocity_pipeline.py -l ../results/seurat_intersect_velocity/${labels[$index]}_vel.loom -s ../results/seurat_intersect_velocity/${labels[$index]}_seu.rds -n ${labels[$index]} --species fish
  #   #break
  #done

  ####################
  ## for human
  ####################
  #for i in `cut -f 1 /PHShome/qq06/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/final_list.txt`; do
  #    #ls -d /PHShome/qq06/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/${i}
  #    #if [ ! -s ../results/seurat/${i}_seurat_obj_tumors.rds ]; then
  #	  #echo ../results/seurat/${i}_hg19_seurat_obj_tumors.rds
  #        echo Rscript seurat_pipeline.R --seuratobj /PHShome/qq06/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/${i}/outs/filtered_feature_bc_matrix --label ${i} --finalres 0.8 --tumor -1 --species human #1>${i}.log 2>&1
  #        #bsub -J humanvelR -n 6 -q big -M 32000 "Rscript ${src}/seurat_pipeline.R --seuratobj /PHShome/qq06/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/${i}/outs/filtered_feature_bc_matrix --label ${i} --finalres 0.8 --tumor -1 --species human 1>${i}.log 2>&1"
  #    # fi
  #done

  ## Generate V6 Seurat results
  #alvins=(`ls ../results/seurat/*hg19_*rds`)
  #saras=(`ls /data/langenau/human_rms_pdxs/seurat_objects/*rds`)
  #labels=(RH74-10cells_5Kcells MAST111_5Kcells MAST139_5Kcells MAST35_5Kcells MAST39_5Kcells MAST85_5Kcells MAST95_10Kcells MSK82489_5Kcells RH74_5Kcells MAST85-1cell_5Kcells)

  #for index in ${!labels[*]}; do
  #    #echo $index ${labels[$index]} $(ls ../results/seurat/*${labels[$index]}_hg19_seurat_obj_tumors.rds) $(ls /data/langenau/human_rms_pdxs/seurat_objects/*$(echo ${labels[$index]} | cut -f 1 -d_ ).rds | tail -1)
  #    echo Rscript generate_v6_degenes.R --seuratobj1  $(ls ../results/seurat/*${labels[$index]}_hg19_seurat_obj_tumors.rds) --seuratobj2 $(ls /data/langenau/human_rms_pdxs/seurat_objects/*$(echo ${labels[$index]} | cut -f 1 -d_ ).rds | tail -1) --label ${labels[$index]}
  #    Rscript generate_v6_degenes.R --seuratobj1  $(ls ../results/seurat/*${labels[$index]}_hg19_seurat_obj_tumors.rds) --seuratobj2 $(ls /data/langenau/human_rms_pdxs/seurat_objects/*$(echo ${labels[$index]} | cut -f 1 -d_ ).rds | tail -1) --label ${labels[$index]}  &
  #done

  # for i in /PHShome/qq06/alvin_singlecell/01_rms_projects/02_human/data/cellranger_counts/*/velocyto/*loom; do
  #     echo $i
  #     label=$(basename $i)
  #     echo Rscript seurat_pipeline.R --seuratobj $i --label $label  --finalres 0.8 --tumor -1 --assaytype spliced --species human # 1>${label}_vel.log 2>&1 # &
  # done

  #Rscript seurat_pipeline.R --seuratobj /data/langenau/human_rms_pdxs/20190801_MAST85-1cell_5Kcells_hg19/velocyto/20190801_MAST85-1cell_5Kcells_hg19.loom --label 20190801_MAST85-1cell --finalres 0.8 --tumor -1 --assaytype spliced --species human 1>20190801_MAST85-1cell.log 2>&1
  #Rscript seurat_pipeline.R --seuratobj /data/langenau/human_rms_pdxs/20190617_RH74-10cells_5Kcells_hg19/velocyto/20190617_RH74-10cells_5Kcells_hg19.loom --label 20190617_RH74-10cells --finalres 0.8 --tumor -1 --assaytype spliced --species human 1>20190617_RH74-10cells.log 2>&1

  #use Sara's processed datasets of Seurat objects /data/langenau/human_rms_pdxs/seurat_objects/
  #Rscript intersect_seurat_velocity_toloom.R --seuratobj /data/langenau/human_rms_pdxs/seurat_objects/20190624_seurat-object_MAST111.rds --velobj ../results/seurat/MAST111_velocity_seurat_obj_tumors.rds --label MAST111 --species human
  # looms=(`ls ../results/seurat/*loom*rds`)
  # saras=(`ls /data/langenau/human_rms_pdxs/seurat_objects/*rds`)
  # n=0
  # for index in ${!saras[*]}; do
  #    prefix=$(echo $(basename ${saras[$index]}) | cut -f 3 -d_ | sed -e 's/.rds//' | cut -f 1 -d-)
  #    #ls ${saras[$index]} ../results/seurat/*${prefix}*rds
  #    label=$(basename ${saras[$index]})
  #    let n++
  #    echo $n
  #    if [ $n -gt 6 ];then
  #        echo Rscript intersect_seurat_velocity_toloom.R --seuratobj ${saras[$index]} --velobj ../results/seurat/*${prefix}*rds --label ${label/.rds/} --species human
  #    #   break
  #    fi
  # done

  #Rscript intersect_seurat_velocity_toloom.R --seuratobj /data/langenau/human_rms_pdxs/seurat_objects/20190815_seurat-object_MAST85-1cell.rds --velobj ../results/seurat/20190801_MAST85-1cell_seurat_obj_tumors.rds --label 20190801_MAST85-1cell --species human
  #Rscript intersect_seurat_velocity_toloom.R --seuratobj /data/langenau/human_rms_pdxs/seurat_objects/20190815_seurat-object_RH74-10cells.rds --velobj ../results/seurat/20190617_RH74-10cells_seurat_obj_tumors.rds --label 20190617_RH74-10cells --species human

  #gene ontology, gsea analysis of above differential cDNA expression
  #n=0
  #for index in ../results/seurat_intersect_velocity/*seu_markers*xls; do
  #    let n++
  #    label=$(basename $index)
  #    #if [ $n -gt 5 ]; then
  #    Rscript gsea_pipeline.R --seuratde $index --label $label
  #    #fi
  #done

  # lisa
  #conda activate lisa
  #cd ../results/gsea/lisa/
  #bash run_lisa.sh

  # R deterministic modeling of human velocity
  #looms=(`ls ../results/seurat_intersect_velocity/*vel.rds`)
  #saras=(`ls ../results/seurat_intersect_velocity/*seu.rds`)
  #for index in ${!saras[*]}; do
  #    # echo $index
  #    # ls ${saras[$index]} ${looms[$index]}
  #    label=$(basename ${saras[$index]})
  #    echo Rscript velocity_pipeline.R --seuratobj ${saras[$index]} --velobj ${looms[$index]} --label ${label/vel.rds/}
  #    #break
  #done

  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST111_seu.rds --velobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST111_vel.rds --label 20190624_seurat-object_MAST111_seu.rds --clusterlabel MAST111 &
  #bsub -J humanvelR -n 6 -q big -M 32000 "bash ${src}/velocity_r.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_intersect_velocity/20190624_seurat-object_MAST111_seu.rds /data/langenau/alvin_singlecell/01_rms_projects/01_fish/results/seurat_intersect_velocity/20190624_seurat-object_MAST111_vel.rds 20190624_seurat-object_MAST111_vel MAST111"
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST139_seu.rds --velobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST139_vel.rds --label 20190624_seurat-object_MAST139_seu.rds --clusterlabel MAST139 &
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST35_seu.rds --velobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST35_vel.rds --label 20190624_seurat-object_MAST35_seu.rds --clusterlabel MAST35 &
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST39_seu.rds --velobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST39_vel.rds --label 20190624_seurat-object_MAST39_seu.rds --clusterlabel MAST39 &
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST85_seu.rds --velobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST85_vel.rds --label 20190624_seurat-object_MAST85_seu.rds
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST95_seu.rds --velobj ../results/seurat_intersect_velocity/20190624_seurat-object_MAST95_vel.rds --label 20190624_seurat-object_MAST95_seu.rds
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/20190624_seurat-object_MSK82489_seu.rds --velobj ../results/seurat_intersect_velocity/20190624_seurat-object_MSK82489_vel.rds --label 20190624_seurat-object_MSK82489_seu.rds
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/20190624_seurat-object_RH74_seu.rds --velobj ../results/seurat_intersect_velocity/20190624_seurat-object_RH74_vel.rds --label 20190624_seurat-object_RH74_seu.rds
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/20190617_RH74-10cells_seu.rds --velobj ../results/seurat_intersect_velocity/20190617_RH74-10cells_vel.rds --label 20190617_RH74-10cells --clusterlabel RH74.10cells
  #Rscript velocity_pipeline.R --seuratobj ../results/seurat_intersect_velocity/20190801_MAST85-1cell_seu.rds --velobj ../results/seurat_intersect_velocity/20190801_MAST85-1cell_vel.rds --label 20190801_MAST85-1cell --clusterlabel MAST85.1cell

  #looms=(`ls ../results/seurat_intersect_velocity/*vel.rds`)
  #saras=(`ls ../results/seurat_intersect_velocity/*seu.rds`)
  #labels=(RH74.10cells MAST111 MAST139 MAST35 MAST39 MAST85 MAST95 MSK82489 RH74 MAST85.1cell Tumor21 Tumor22 Tumor24)
  #for index in ${!saras[*]}; do
  #    # echo $index
  #    ls ${saras[$index]} ${looms[$index]}
  #    #label=$(basename ${saras[$index]})
  #    echo bsub -J humanvelR -n 6 -q big -M 32000 "bash ${src}/velocity_r.sh ${saras[$index]} ${looms[$index]} ${labels[$index]}_velR ${labels[$index]}"
  #    #bsub -J humanvelR -n 6 -q big -M 32000 "bash ${src}/velocity_r.sh ${saras[$index]} ${looms[$index]} ${labels[$index]}_velR ${labels[$index]}"
  #    #break
  #done

  #1. paga for layout the umap, scvelo for identification of root cells and trajectories by both velocity and stream
  #2. stream for trajectory analysis by using scvelo estimated starting points
  #look at scanpy_velocity_stream.ipynb for visualization of velocity and stream trajectories
  #looms=(`ls ../results/seurat_intersect_velocity/*loom`)
  #saras=(`ls ../results/seurat_intersect_velocity/*seu.rds`)
  #labels=(RH74-10cells MAST111 MAST139 MAST35 MAST39 MAST85 MAST95 MSK82489 RH74 MAST85-1cell)
  #mkdir -p ../results/velocity_dynamical/
  #for index in ${!saras[*]}; do
  #    label=$(basename ${looms[$index]})
  #    #echo python velocity_pipeline.py -l ${looms[$index]} -s ${saras[$index]} -n ${labels[$index]} #2>&1 > ${label/.loom/}_scvelo.logs & # output anndata for stream analysis
  #    echo python velocity_pipeline_dynamical_latentime.py -l ${looms[$index]} -s ${saras[$index]} -n ${labels[$index]} #2>&1 > ${label/.loom/}_scvelo.logs & # output anndata for stream analysis
  #    #python velocity_pipeline_dynamical_latentime.py -s ${saras[$index]} --name ${labels[$index]}
  #    #break
  #done

  ##cd ../results/gsea/lisa/
  #${src}/run_genewalk.sh

  ## process new datasets
  #cd /data/langenau/human_rms_pdxs/20191031_MSK72117tencell
  #echo 'Lane,Sample,Index' > 20191031_MSK72117tencell.csv
  #echo '1-4,20191031_MSK72117tencell,SI-GA-B8' >> 20191031_MSK72117tencell.csv
  #cellranger mkfastq --id=20191031_MSK72117tencell \
  #                   --localcores=8 \
  #                   --localmem=64 \
  #                   --run=/data/langenau/human_rms_pdxs/20191031_MSK72117tencell \
  #                   --csv=20191031_MSK72117tencell.csv

  #bash ${src}/cellranger.sh /data/langenau/human_rms_pdxs/20191031_MSK72117tencell/20191031_MSK72117tencell/outs/fastq_path/HL3F7BGXC/ 20191031_MSK72117tencell 20191031_MSK72117tencell_cellranger_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/
  #bash ${src}/cellranger.sh /data/langenau/human_rms_pdxs/20191031_MSK72117tencell/20191031_MSK72117tencell/outs/fastq_path/HL3F7BGXC/ 20191031_MSK72117tencell 20191031_MSK72117tencell_cellranger /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/
  #low cell ratio, rerun with 5K forced cells
  #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ${src}/cellranger.sh /data/langenau/human_rms_pdxs/20191031_MSK72117tencell/20191031_MSK72117tencell/outs/fastq_path/HL3F7BGXC/ 20191031_MSK72117tencell 20191031_MSK72117tencell_cellranger_5K /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/"
  #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ${src}/cellranger.sh /data/langenau/human_rms_pdxs/20191031_MSK72117tencell/20191031_MSK72117tencell/outs/fastq_path/HL3F7BGXC/ 20191031_MSK72117tencell 20191031_MSK72117tencell_cellranger_5K_mixture /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/"


  #cd /data/langenau/human_rms_pdxs/20191031_MSK74711/
  #echo 'Lane,Sample,Index' > 20191031_MSK74711.csv
  #echo '1-4,20191031_MSK74711,SI-GA-B9' >> 20191031_MSK74711.csv
  #cellranger mkfastq --id=20191031_MSK74711 \
  #                   --localcores=8 \
  #                   --localmem=64 \
  #                   --run=/data/langenau/human_rms_pdxs/20191031_MSK74711 \
  #                   --csv=20191031_MSK74711.csv

  #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ${src}/cellranger.sh /data/langenau/human_rms_pdxs/20191031_MSK74711/20191031_MSK74711/outs/fastq_path/HL2THBGXC/ 20191031_MSK74711 20191031_MSK74711_cellranger_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/"
  #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ${src}/cellranger.sh /data/langenau/human_rms_pdxs/20191031_MSK74711/20191031_MSK74711/outs/fastq_path/HL2THBGXC/ 20191031_MSK74711 20191031_MSK74711_cellranger /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/

  #run through sara's pipeline
  #Rscript seurat_sara_pipeline.R --seuratobj "../data/cellranger_counts/20191031_MSK74711_cellranger/outs/filtered_feature_bc_matrix"  \
  #                               --mixtureobj "../data/cellranger_counts/20191031_MSK74711_cellranger_mixture/outs/filtered_feature_bc_matrix" --label 20191031_MSK74711

  #Rscript seurat_sara_pipeline.R --seuratobj "../data/cellranger_counts/20191031_MSK72117tencell_cellranger/outs/filtered_feature_bc_matrix/" \
  #                               --mixtureobj "../data/cellranger_counts/20191031_MSK72117tencell_cellranger_mixture/outs/filtered_feature_bc_matrix/" --label 20191031_MSK72117tencell

  #mkdir -p ../data/cellranger_counts/MAST139_1cells
  #cd ../data/cellranger_counts/MAST139_1cells
  #bs download run --id 191551362 --output .

  #echo 'Lane,Sample,Index' > MAST139_1cells.csv
  #echo '1-4,MAST139_1cells,SI-GA-C7' >> MAST139_1cells.csv
  #cellranger mkfastq --id=MAST139_1cells \
  #                   --localcores=8 \
  #                   --localmem=64 \
  #                   --run=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/MAST139_1cells \
  #                   --csv=MAST139_1cells.csv

   #cd ~/scratch/
   #echo 'Lane,Sample,Index' > MAST118_final.csv
   #echo '1-4,MAST118,SI-GA-B6' >> MAST118_final.csv
   #cellranger mkfastq --id=MAST118 \
   #                  --localcores=20 \
   #                  --localmem=64 \
   #                  --run=/data/langenau/human_rms_pdxs/20191203_MAST118_FINAL \
   #                  --csv=MAST118_final.csv

   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/MAST139_1cells/MAST139_1cells/outs/fastq_path/HL73YBGXC/MAST139_1cell/ Tumor24 MAST139_1cells /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/"
   #bash ../../rms_codes/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/MAST139_1cells/MAST139_1cells/outs/fastq_path/HL73YBGXC/MAST139_1cell/ Tumor24 MAST139_1cells /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/
   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/MAST139_1cells_fastq/MAST139_1cells/outs/fastq_path/HL73YBGXC/MAST139_1cell/ Tumor24 MAST139_1cells_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/"

   #Rscript seurat_sara_pipeline.R --seuratobj "../data/cellranger_counts/MAST139_1cells/outs/filtered_feature_bc_matrix/" \
   #                               --mixtureobj "../data/cellranger_counts/MAST139_1cells_mixture/outs/filtered_feature_bc_matrix/" --label MAST139_1cells

   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash /data/langenau/alvin_singlecell/01_rms_projects/rms_codes/cellranger.sh /PHShome/qq06/scratch/MAST118/outs/fastq_path/HV35LBGXC/MAST118/ MAST118 MAST118_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/"
   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash /data/langenau/alvin_singlecell/01_rms_projects/rms_codes/cellranger.sh /PHShome/qq06/scratch/MAST118/outs/fastq_path/HV35LBGXC/MAST118/ MAST118 MAST118_hg19 /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/"

   #Rscript seurat_sara_pipeline.R --seuratobj "/data/langenau/human_rms_pdxs/20191203_MAST118_FINAL/MAST118_hg19/outs/filtered_feature_bc_matrix" \
   #                               --mixtureobj "/data/langenau/human_rms_pdxs/20191203_MAST118_FINAL/MAST118_mixture/outs/filtered_feature_bc_matrix" --label MAST118

   #for i in ../results/seurat_sara/*rds; do
   #for i in ../results/seurat_sara/MAST139_1cells*rds; do
   #for i in ../data/seurat_obj/*MAST85-1cell*rds; do
   #for i in ../data/seurat_obj/2*rds; do
   # for i in ../results/seurat_sara/*MSK74711*rds ../results/seurat_sara/*MAST118*rds ../results/seurat_sara/*MSK72117*rds; do
   #     label=$(basename $i)
   #     #conda activate alvin_resvel
   #     #Rscript generate_v4_degenes.R --seuratobj1 $i --label ${label/.rds} &
   #     #conda activate sc-tutorial
   #     Rscript annotate_celltypes.R --seuratobj $i --label ${label/.rds} &
   #     #break
   # done

   #for i in MAST139_1cells; do
   #for i in MAST139_1cells; do
   ##for i in 20191031_MSK72117tencell_cellranger; do
   ##for i in 20191031_MSK74711_cellranger; do
   #   velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i} /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/genes/genes.gtf # 2>&1 > ${i}.logs &
   #done
   #velocyto run10x /data/langenau/human_rms_pdxs/20191203_MAST118_FINAL/MAST118_hg19 /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/genes/genes.gtf # 2>&1 > ${i}.logs &

   #for i in 20191031_MSK72117tencell_cellranger 20191031_MSK74711_cellranger MAST139_1cells; do
   #    Rscript seurat_pipeline.R --seuratobj /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}/velocyto/* --label $i --finalres 0.8 --tumor -1 --assaytype spliced --species human
   #done

   # Analyze the RD cell lines
   #echo 'Lane,Sample,Index' > RD.csv
   #echo '1-4,RH41,SI-GA-D3' >> RD.csv
   #echo '1-4,RD,SI-GA-D4' >> RD.csv
   #cellranger mkfastq --id=RMS_cellline \
   #                   --localcores=20 \
   #                   --localmem=64 \
   #                   --run=/data/langenau/human_rms_pdxs/20200109_RD_RH41_FINAL \
   #                   --csv=RD.csv

   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/RMS_cellline/outs/fastq_path/HTHHFBGX9/ RD RD /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/"
   # bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/RMS_cellline/outs/fastq_path/HTHHFBGX9/ RD RD_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0"
   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/RMS_cellline/outs/fastq_path/HTHHFBGX9/ RH41 RH41 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/"
   # bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/RMS_cellline/outs/fastq_path/HTHHFBGX9/ RH41 RH41_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0"

   # echo 'Lane,Sample,Index' > MAST139_Muscle.csv
   # echo '1-4,MAST139_Muscle_plus,SI-GA-E1' >>  MAST139_Muscle.csv
   # echo '1-4,MAST139_Muscle_minus,SI-GA-E2' >> MAST139_Muscle.csv
   # cellranger mkfastq --id=MAST139_Muscle \
   #                    --localcores=20 \
   #                    --localmem=64 \
   #                    --run=/data/langenau/human_rms_pdxs/2020_01_MAST139again \
   #                    --csv=MAST139_Muscle.csv

   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/MAST139_Muscle MAST139_Muscle_plus MAST139_Muscle_plus /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/"
   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/MAST139_Muscle MAST139_Muscle_plus MAST139_Muscle_plus_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/"
   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/MAST139_Muscle MAST139_Muscle_minus MAST139_Muscle_minus /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/"
   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/MAST139_Muscle MAST139_Muscle_minus MAST139_Muscle_minus_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/"

   #echo 'Lane,Sample,Index' > MSK93202.csv
   #echo '1-4,MSK93202,SI-GA-B7' >> MSK93202.csv
   #cellranger mkfastq --id=MSK93202 \
   #                   --localcores=20 \
   #                   --localmem=64 \
   #                   --run=/data/langenau/human_rms_pdxs/202001_MSK93202_FINAL \
   #                   --csv=MSK93202.csv

   # bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/MSK93202/ MSK93202 MSK93202_hg19 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/"
   # bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/MSK93202/ MSK93202 MSK93202_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/"

   #for sample in MAST139_Muscle_minus MAST139_Muscle_plus; do
   #    Rscript seurat_sara_pipeline.R --seuratobj "${sample}/outs/filtered_feature_bc_matrix" \
   #                                   --mixtureobj "${sample}_mixture/outs/filtered_feature_bc_matrix" --label $sample
   #done
   #for sample in MSK93202; do
   #       Rscript seurat_sara_pipeline.R --seuratobj "${sample}_hg19/outs/filtered_feature_bc_matrix" \
   #                                  --mixtureobj "${sample}_mixture/outs/filtered_feature_bc_matrix" --label $sample
   #done

   # for sample in RD RH41 /data/pinello/PROJECTS/2019_11_ResidualVelocity/data/CD44_rhabodo/SRR9927169_2/ /data/pinello/PROJECTS/2019_11_ResidualVelocity/data/CD44_rhabodo/SRR9927170_2; do
   # 	echo $sample
   #     Rscript seurat_sara_pipeline.R --seuratobj "${sample}/outs/filtered_feature_bc_matrix" \
   #                                    --mixtureobj NA --label $(basename $sample)
   # done

   # for i in MAST139_Muscle_minus MAST139_Muscle_plus MSK93202_hg19 RD RH41; do
   # 	ls $i
   # 	velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/${i} /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/genes/genes.gtf & # 2>&1 > ${i}.logs &
   # done

   #for i in ../results/seurat_sara/MSK93202_seurat-object.rds ../results/seurat_sara/MAST139_Muscle_plus_seurat-object.rds ../results/seurat_sara/MAST139_Muscle_minus_seurat-object.rds ../results/seurat_sara/RD_seurat-object.rds ../results/seurat_sara/RH41_seurat-object.rds ../results/seurat_sara/SRR9927169_2_seurat-object.rds ../results/seurat_sara/SRR9927170_2_seurat-object.rds; do
   #    echo $i
   #    label=$(basename $i)
   #    Rscript generate_v4_degenes.R --seuratobj1 $i --label ${label/.rds} &
   #done

   #primary tumor cells
   #echo 'Lane,Sample,Index' > 20696.csv
   #echo '1-4,20696,SI-GA-E5' >> 20696.csv
   #cellranger mkfastq --id=20696 \
   #                   --localcores=20 \
   #                   --localmem=64 \
   #                   --run=/data/langenau/human_rms_pdxs/20200131_20696/ \
   #                   --csv=20696.csv

   #single nucleus sequencing premrna index
   #https://www.biostars.org/p/461672/
   #awk 'BEGIN{FS="\t"; OFS="\t"} $3 == "transcript"{ $3="exon"; print}' \
   #       /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/genes/genes.gtf > hg19-3.0.0.premrna.gtf
   #cellranger mkref --genome=hg19-3.0.0.premrna \
   #                 --fasta=/data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/fasta/genome.fa \
   #                 --genes=hg19-3.0.0.premrna.gtf

   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/src/20696 20696 20696_hg19 /data/langenau/alvin_singlecell/01_rms_projects/01_fish/src/hg19-3.0.0.premrna/"
   #velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/20696_hg19/ /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/genes/genes.gtf

   #Rscript seurat_sara_pipeline.R --seuratobj "20696_hg19/outs/filtered_feature_bc_matrix" \
   #                               --mixtureobj NA --label 20696
   #Rscript generate_v4_degenes.R --seuratobj1 ../results/seurat_sara/20696_seurat-object.rds \
   #                               --label 20696

   #python identify_doublet.py -mat 20696_hg19/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature 20696_hg19/outs/filtered_feature_bc_matrix/features.tsv.gz -name 20696 -gzip
    #Rscript seurat_pipeline.R --seuratobj 20696_hg19/velocyto/* --label 20696 --finalres 0.8 --tumor -1 --assaytype spliced --species human

    #Rscript intersect_seurat_velocity_toloom.R --seuratobj ../results/seurat_sara/20696_seurat-object.rds --velobj ../results/seurat/20696_seurat_obj_tumors.rds --label 20696 --species human
    #python velocity_pipeline_dynamical_latentime.py -l ../results/seurat_intersect_velocity/20696_vel.loom -s ../results/seurat_intersect_velocity/20696_seu.rds -n 20696_dynamical_model

   #echo 'Lane,Sample,Index' > 29806_primary.csv
   #echo '1-4,29806,SI-GA-E6' >>  29806_primary.csv
   #cellranger mkfastq --id=29806 \
   #                   --localcores=20 \
   #                   --localmem=64 \
   #                   --run=/data/langenau/human_rms_pdxs/20200226_29806 \
   #                   --csv=29806_primary.csv

   #echo 'Lane,Sample,Index' > 20082_primary.csv
   #echo '1-4,20082,SI-GA-E7' >> 20082_primary.csv
   #cellranger mkfastq --id=20082 \
   #                   --localcores=20 \
   #                   --localmem=64 \
   #                   --run=/data/langenau/human_rms_pdxs/20200227_20082 \
   #                   --csv=20082_primary.csv

   #echo 'Lane,Sample,Index' > 21202_primary.csv
   #echo '1-4,21202,SI-GA-E8' >> 21202_primary.csv
   #cellranger mkfastq --id=21202 \
   #                   --localcores=20 \
   #                   --localmem=64 \
   #                   --run=/data/langenau/human_rms_pdxs/20200301_21202 \
   #                   --csv=21202_primary.csv

   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/20082 20082 20082_hg19 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/"
   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/29806 29806 29806_hg19 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/"

   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/20082 20082 20082_hg19_premrna /data/langenau/alvin_singlecell/01_rms_projects/01_fish/src/hg19-3.0.0.premrna/"
   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/29806 29806 29806_hg19_premrna /data/langenau/alvin_singlecell/01_rms_projects/01_fish/src/hg19-3.0.0.premrna/"
   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/21202 21202 21202_hg19_premrna /data/langenau/alvin_singlecell/01_rms_projects/01_fish/src/hg19-3.0.0.premrna/"

  # ##for i in 20082_hg19_premrna 29806_hg19_premrna 21202_hg19_premrna; do
  # #     echo $i
  # #     #python identify_doublet.py -mat ${i}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${i}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -name ${i}_doublet -gzip
  # #     #velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/src/${i} /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/genes/genes.gtf &
  # #     #Rscript seurat_sara_pipeline.R --seuratobj "${src}/${i}/outs/filtered_feature_bc_matrix" \
  # #     #                               --mixtureobj NA --label ${i} --doublet ../results/doublets/${i}_doublet_doublet.csv
  # #     #Rscript seurat_pipeline.R --seuratobj ${i}/velocyto/* --label ${i} --finalres 0.8 --tumor -1 --assaytype spliced --species human &
  #      Rscript intersect_seurat_velocity_toloom.R --seuratobj ../results/seurat_sara/${i}_seurat-object.rds --velobj ../results/seurat/${i}_seurat_obj_tumors.rds --label ${i} --species human
  # done

  #for i in 20082_hg19 20696_hg19 21202_hg19 29806_hg19;do
  #    ## de analysis on tumor only Fig1_integration_recluster_individual_tumor.R
  #    #echo Rscript annotate_celltypes.R --seuratobj ../figures/${i}_*umap.rds --label ${i}
  #    Rscript annotate_celltypes.R --seuratobj ../figures/${i}_*umap.rds --label ${i} &
  #done
  #Rscript annotate_celltypes.R --seuratobj 20082_recluster2_tumor_only.rds --label 20082_recluster2

   #the rest of velocity analysis
   # Rscript seurat_pipeline.R --seuratobj /data/langenau/human_rms_pdxs/20191203_MAST118_FINAL/MAST118_hg19/velocyto/* --label MAST118 --finalres 0.8 --tumor -1 --assaytype spliced --species human
   #for i in 20191031_MSK72117tencell_cellranger 20191031_MSK74711_cellranger MAST139_1cells; do
   #    echo $i
   #    Rscript seurat_pipeline.R --seuratobj /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}/velocyto/* --label $i --finalres 0.8 --tumor -1 --assaytype spliced --species human
   #done
   #Rscript seurat_pipeline.R --seuratobj /data/langenau/human_rms_pdxs//20190801_MAST85-1cell_5Kcells_hg19/velocyto/* --label MAST85-1cell --finalres 0.8 --tumor -1 --assaytype spliced --species human --transform SCT #RNA
   #Rscript seurat_pipeline.R --seuratobj /data/langenau/human_rms_pdxs/20190617_RH74-10cells_5Kcells_hg19/velocyto/* --label RH74-10cells --finalres 0.8 --tumor -1 --assaytype spliced --species human  --transform SCT #RNA

   #seurat_sara_left=('../results/seurat_sara/MAST118_seurat-object.rds' '../results/seurat_sara/20191031_MSK74711_seurat-object.rds' '../results/seurat_sara/20191031_MSK72117tencell_seurat-object.rds' '../results/seurat_sara/MAST139_1cells_seurat-object.rds')
   #labels=(MAST118 20191031_MSK74711_cellranger 20191031_MSK72117tencell_cellranger MAST139_1cells)
   #seurat_sara_left=('../data/seurat_obj/20190815_seurat-object_MAST85-1cell.rds' '../data/seurat_obj/20190815_seurat-object_RH74-10cells.rds')
   #labels=(MAST85-1cell RH74-10cells)
   #seurat_sara_left=('../data/seurat_obj/20190815_seurat-object_MAST85-1cell.rds')
   #labels=(MAST85-1cell)
   #seurat_sara_left=(../results/seurat_sara/MAST139_1cells_seurat-object.rds)
   #labels=(MAST139_1cells)
   ##seurat_sara_left=('../figures/20696_hg19_tumoronly_res0.8_umap.rds' '20082_recluster2_tumor_only.rds' '../figures/29806_hg19_premrna_tumoronly_res0.8_umap.rds' '../figures/21202_hg19_premrna_tumoronly_res0.8_umap.rds')
   ##labels=(20696 20082_hg19_premrna 29806_hg19_premrna 21202_hg19_premrna)
   #for i in `seq 0 $((${#seurat_sara_left[@]}-1))`; do
       #echo Rscript intersect_seurat_velocity_toloom.R --seuratobj ${seurat_sara_left[$i]} --velobj ../results/seurat/${labels[$i]}_seurat_obj_tumors.rds --label ${labels[$i]}_vel --species human
       #Rscript intersect_seurat_velocity_toloom.R --seuratobj ${seurat_sara_left[$i]} --velobj ../results/seurat/${labels[$i]}_seurat_obj_tumors.rds --label ${labels[$i]}_vel --species human
       ##python velocity_pipeline_dynamical_latentime.py -l ../results/seurat_intersect_velocity/${labels[$i]}_vel.loom -s ../results/seurat_intersect_velocity/${labels[$i]}_seu.rds -n ${labels[$i]}_dynamical_model
       #python velocity_pipeline_dynamical_latentime.py -l ../results/seurat_intersect_velocity/${labels[$i]}_vel_vel.loom -s ../results/seurat_intersect_velocity/${labels[$i]}_vel_seu.rds -n ${labels[$i]} #_dynamical_model
       #echo python velocity_pipeline.py -l ../results/seurat_intersect_velocity/${labels[$i]}_vel_vel.loom -s ../results/seurat_intersect_velocity/${labels[$i]}_vel_seu.rds -n ${labels[$i]}_dynamical_model
       #python velocity_pipeline.py -l ../results/seurat_intersect_velocity/${labels[$i]}_vel_vel.loom -s ../results/seurat_intersect_velocity/${labels[$i]}_vel_seu.rds -n ${labels[$i]}_dynamical_model
       #echo python velocity_pipeline_steady.py -l ../results/seurat_intersect_velocity/${labels[$i]}_vel_vel.loom -s ../results/seurat_intersect_velocity/${labels[$i]}_vel_seu.rds -n ${labels[$i]}
       #python velocity_pipeline_steady.py -l ../results/seurat_intersect_velocity/${labels[$i]}_vel_vel.loom -s ../results/seurat_intersect_velocity/${labels[$i]}_vel_seu.rds -n ${labels[$i]}_steady
   #done

   #seurat_sara_left=(../data/seurat_obj/20190624_seurat-object_MAST111.rds ../data/seurat_obj/20190624_seurat-object_MAST139.rds)
   #labels=(20190624_seurat-object_MAST111 20190624_seurat-object_MAST139)
   #seurat_sara_left=(' ../data/seurat_obj/20190624_seurat-object_MAST139.rds')
   #labels=(20190624_seurat-object_MAST139)
   #seurat_sara_left=(../data/seurat_obj/20190624_seurat-object_MAST35.rds ../data/seurat_obj/20190624_seurat-object_MAST39.rds ../data/seurat_obj/20190624_seurat-object_MAST85.rds ../data/seurat_obj/20190624_seurat-object_MSK82489.rds ../data/seurat_obj/20190624_seurat-object_MAST95.rds ../data/seurat_obj/20190624_seurat-object_RH74.rds)
   #labels=(20190624_seurat-object_MAST35 20190624_seurat-object_MAST39 20190624_seurat-object_MAST85 20190624_seurat-object_MSK82489 20190624_seurat-object_MAST95 20190624_seurat-object_RH74)
   #for i in `seq 0 $((${#seurat_sara_left[@]}-1))`; do
   #    echo python velocity_pipeline_steady.py -l ../results/seurat_intersect_velocity/${labels[$i]}_vel.loom -s ../results/seurat_intersect_velocity/${labels[$i]}_seu.rds -n ${labels[$i]}
   #    python velocity_pipeline_steady.py -l ../results/seurat_intersect_velocity/${labels[$i]}_vel.loom -s ../results/seurat_intersect_velocity/${labels[$i]}_seu.rds -n ${labels[$i]}_steady &
   #done

  #for i in /PHShome/qq06/langenau/01_rms_projects/01_fish/src/20082_hg19_premrna/ \
  #         /PHShome/qq06/langenau/01_rms_projects/01_fish/src/20696_hg19/ \
  #         /PHShome/qq06/langenau/01_rms_projects/01_fish/src/21202_hg19_premrna/ \
  #         /PHShome/qq06/langenau/01_rms_projects/01_fish/src/29806_hg19_premrna/; do
  #    label=$(basename $i)
  #    python identify_doublet.py -mat ${i}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${i}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -name ${label}_doublet -gzip
  #done

   #for i in /PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color_v2/ \
   #         /PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish_with_orf_color_v2 /PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2/ \
   #         /PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish_with_orf_color_v2 /PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2/; do
   #    label=$(basename $i)
   #    python identify_doublet.py -mat ${i}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature ${i}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -name ${label}_doublet -gzip
   #done

   #### do not re-process zebrafish dataset with Sara version 4 pipeline without SCTtransform
   #### for i in /PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color_v2/ \
   ####          "/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish_with_orf_color_v2 /PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2/" \
   ####         "/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish_with_orf_color_v2 /PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2/"; do
   ####     echo $i
   ####     label=$(basename $i)
   ####     #Rscript seurat_sara_pipeline.R --seuratobj ${i}/outs/filtered_feature_bc_matrix \
   ####     #                               --mixtureobj NA --label ${i} --doublet ../results/doublets/${i}_doublet_doublet.csv
   ####     #Rscript seurat_sara_pipeline.R --seuratobj '/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix/' --label Tumor24_unfilter_v4 --finalres 0.05 --tumor -1
   ####     #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' '/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' --label Tumor21_unfilter_v4 --finalres 0.1 --tumor -1 &
   ####     #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' '/PHShome/qq06/langenau/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2/outs/filtered_feature_bc_matrix' --label Tumor22_unfilter_v4 --finalres 0.1 --tumor -1 &
   #### done

   #genome=/data/pinello/PROJECTS/2019_11_ResidualVelocity/data/latest_kb_index/human/nucleus/
   ####test=$(ls /data/langenau/alvin_singlecell/01_rms_projects/01_fish/src/20082/outs/fastq_path/HFYWWBGXF/20082/*R1*)
   ####n=0
   ####for i in $test $test2; do
   ####    dirname=$(basename $(dirname $i))
   ####    echo kb count --h5ad -i ${genome}/index.idx -g ${genome}/t2g.txt -x 10XV2 -o . -c1 ${genome}/cdna_t2c.txt -c2 ${genome}/intron_t2c.txt --workflow nucleus --filter bustools -t 8 ${i} ${i/R1/R2}
   ####    label=${i}
   ####    label=$(basename $label)
   ####    label=${label/.fastq.gz/}
   ####    mkdir -p ${dirname}_$label
   ####    cd ${dirname}_${label}
   ####    #bsub -q rerunnable kb count --h5ad -i ${genome}/index.idx -g ${genome}/t2g.txt -x 10XV2 -o . -c1 ${genome}/cdna_t2c.txt -c2 ${genome}/intron_t2c.txt --workflow nucleus --filter bustools -t 1 ../${i} ../${i/R1/R2}
   ####    #kb count --h5ad -i ${genome}/index.idx -g ${genome}/t2g.txt -x 10XV2 -o . -c1 ${genome}/cdna_t2c.txt -c2 ${genome}/intron_t2c.txt --workflow nucleus --filter bustools -t 8 ${i} ${i/R1/R2}
   ####    if [ $n -gt 0 ]; then
   ####        kb count --h5ad -i ${genome}/index.idx -g ${genome}/t2g.txt -x 10XV3 -o . -c1 ${genome}/cdna_t2c.txt -c2 ${genome}/intron_t2c.txt --workflow nucleus --filter bustools -t 8 ${i} ${i/R1/R2}
   ####    fi
   ####    cd -
   ####    #break
   ####    let n++
   ####done
   #snakemake -j 8 -s Snakefile.10x --config path=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/src/20082/outs/fastq_path/HFYWWBGXF/20082/
   #snakemake -j 8 -s Snakefile.10x --config path=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/src/29806/outs/fastq_path/HFYH7BGXF/29806/
   #snakemake -j 8 -s Snakefile.10x --config path=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/src/21202/outs/fastq_path/HFYWGBGXF/21202/
   #snakemake -j 8 -s Snakefile.10x --config path=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/src/20696/outs/fastq_path/H5T5FBGXF/20696/

   #bs download run --id 195747619 -o C12SC1_FINAL
   # echo 'Lane,Sample,Index' > C12SC1_FINAL.csv
   # echo '1-4,C12SC1,SI-GA-E9' >> C12SC1_FINAL.csv
   # #cellranger mkfastq --id=C12SC1 \
   # #                   --localcores=8 \
   # #                   --localmem=64 \
   # #                   --run=C12SC1_FINAL \
   # #                   --csv=C12SC1_FINAL.csv
   # #bsub -q rerunnable bash ${src}/cellranger.sh /PHShome/qq06/projects/01_sc_rms/phaseA_explore_rms/C12SC1/ C12SC1 C12SC1_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/
   # cd /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/
   # #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash /PHShome/qq06/projects/01_sc_rms/phaseA_explore_rms/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/C12SC1/ C12SC1 C12SC1_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/"
   # #cd /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/
   # #bsub -q big-multi bash /PHShome/qq06/projects/01_sc_rms/phaseA_explore_rms/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/C12SC1/ C12SC1 C12SC1_hg19 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/
   # bash /PHShome/qq06/projects/01_sc_rms/phaseA_explore_rms/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/C12SC1/ C12SC1 C12SC1_hg19_run2 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/

   #for i in C12SC1_hg19; do
   for i in C12SC1; do
        echo $i
        #python identify_doublet.py -mat /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -name ${i}_doublet -gzip
        #velocyto run10x /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i} /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/genes/genes.gtf &
        #Rscript seurat_sara_pipeline.R --seuratobj "/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}_hg19/outs/filtered_feature_bc_matrix" \
        #                               --mixtureobj /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}_mixture/outs/filtered_feature_bc_matrix --label ${i} --doublet ../results/${i}_hg19_doublet_doublet.csv
        #Rscript seurat_pipeline.R --seuratobj /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellran/${i}_hg19/outs/velocyto/* --label ${i} --finalres 0.8 --tumor -1 --assaytype spliced --species human
        #Rscript intersect_seurat_velocity_toloom.R --seuratobj ../results/seurat_sara/${i}_seurat-object.rds --velobj ../results/seurat/${i}_seurat_obj_tumors.rds --label ${i} --species human
        #Rscript generate_v4_degenes.R --seuratobj1 ../results/seurat_sara/C12SC1_seurat-object.rds --label ${i} #&
        #Rscript annotate_celltypes.R --seuratobj ../results/seurat_sara/C12SC1_seurat-object.rds --label ${i} #&
   done
   # bs download run --id 196240071 -o C12SC2_FINAL
   # echo 'Lane,Sample,Index' > C12SC2_FINAL.csv
   # echo '1-4,C12SC2,SI-GA-E10' >> C12SC2_FINAL.csv
   #cellranger mkfastq --id=C12SC2 \
   #                   --localcores=8 \
   #                   --localmem=64 \
   #                   --run=C12SC2_FINAL \
   #                   --csv=C12SC2_FINAL.csv
   #bsub -n 8 -q big-multi -M 32000 -q big-multi bash /PHShome/qq06/projects/01_sc_rms/phaseA_explore_rms/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/C12SC2/ C12SC2 C12SC2_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/
   #bsub -n 8 -q big-multi -M 32000 -q big-multi bash /PHShome/qq06/projects/01_sc_rms/phaseA_explore_rms/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/C12SC2/ C12SC2 C12SC2_hg19 /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/

   #for i in C12SC2_hg19; do
   # for i in C12SC2; do
   for i in RD; do
        echo $i
        #python identify_doublet.py -mat /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}/outs/filtered_feature_bc_matrix/matrix.mtx.gz -feature /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -name ${i}_doublet -gzip
        #velocyto run10x /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i} /data/molpath/software/10x/refdata-cellranger-hg19-3.0.0/genes/genes.gtf

        #Rscript seurat_sara_pipeline.R --seuratobj "/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}_hg19/outs/filtered_feature_bc_matrix" \
        #                               --mixtureobj "/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}_mixture/outs/filtered_feature_bc_matrix" --label ${i} --doublet ../results/${i}_hg19_doublet_doublet.csv
        #Rscript seurat_pipeline.R --seuratobj /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}_hg19/velocyto/* --label ${i} --finalres 0.8 --tumor -1 --assaytype spliced --species human
        #Rscript intersect_seurat_velocity_toloom.R --seuratobj ../results/seurat_sara/${i}_seurat-object.rds --velobj ../results/seurat/${i}_seurat_obj_tumors.rds --label ${i} --species human
        #Rscript generate_v4_degenes.R --seuratobj1 ../results/seurat_sara/${i}_seurat-object.rds --label ${i} #&
        Rscript annotate_celltypes.R --seuratobj ../results/seurat_sara/${i}_seurat-object.rds --label ${i} #&
        #python velocity_pipeline_dynamical_latentime.py -l ../results/seurat_intersect_velocity/${labels[$i]}_vel.loom -s ../results/seurat_intersect_velocity/${labels[$i]}_seu.rds -n ${labels[$i]}_dynamical_model
   done
}

main

