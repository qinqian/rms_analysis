#!/bin/bash -ex

#module load cellranger/3.0.2
data=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data
tmp=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/tmp
genome=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/genome
src=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/src
flurescent=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/flurescent_colors/
export PATH=/data/langenau/alvin_singlecell/01_rms_projects/01_fish/software/cellranger-3.1.0/:${PATH}

mkdir -p $data
mkdir -p $tmp
mkdir -p $src
mkdir -p $genome

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

  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor24_zebrafish_with_orf_color_v2/velocyto/Tumor24_zebrafish_with_orf_color_v2.loom' --label Tumor24_velocity --finalres 0.05 --tumor -1 --assaytype spliced 
  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_zebrafish_with_orf_color_v2/velocyto/Tumor21_zebrafish_with_orf_color_v2.loom' '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2/velocyto/Tumor21_2ndlibrary_zebrafish_with_orf_color_v2.loom' --label Tumor21_velocity --finalres 0.1 --tumor -1 --assaytype spliced 
  #Rscript seurat_pipeline.R --seuratobj '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_zebrafish_with_orf_color_v2/velocyto/Tumor22_zebrafish_with_orf_color_v2.loom' '/PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2/velocyto/Tumor22_2ndlibrary_zebrafish_with_orf_color_v2.loom' --label Tumor22_velocity --finalres 0.1 --tumor -1 --assaytype spliced

  #intersect between seurat and velocity cell barcodes
  #Rscript intersect_seurat_velocity_toloom.R --seuratobj ../results/seurat/Tumor24_seurat_obj_tumors.rds --velobj '../results/seurat/Tumor24_velocity_seurat_obj_tumors.rds' --label Tumor24 #2>&1 >Tumor24_intersect.log 
  #Rscript intersect_seurat_velocity_toloom.R --seuratobj ../results/seurat/Tumor21_seurat_obj_tumors.rds --velobj '../results/seurat/Tumor21_velocity_seurat_obj_tumors.rds' --label Tumor21 #2>&1 >Tumor21_intersect.log 
  #Rscript intersect_seurat_velocity_toloom.R --seuratobj ../results/seurat/Tumor22_seurat_obj_tumors.rds --velobj '../results/seurat/Tumor22_velocity_seurat_obj_tumors.rds' --label Tumor22 # 2>&1 >Tumor22_intersect.log &

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
  looms=(`ls ../results/seurat_intersect_velocity/*loom`)
  saras=(`ls ../results/seurat_intersect_velocity/*seu.rds`)
  labels=(RH74-10cells MAST111 MAST139 MAST35 MAST39 MAST85 MAST95 MSK82489 RH74 MAST85-1cell)
  mkdir -p ../results/velocity_dynamical/
  for index in ${!saras[*]}; do
      label=$(basename ${looms[$index]})
      #echo python velocity_pipeline.py -l ${looms[$index]} -s ${saras[$index]} -n ${labels[$index]} #2>&1 > ${label/.loom/}_scvelo.logs & # output anndata for stream analysis
      echo python velocity_pipeline_dynamical_latentime.py -l ${looms[$index]} -s ${saras[$index]} -n ${labels[$index]} #2>&1 > ${label/.loom/}_scvelo.logs & # output anndata for stream analysis
      #python velocity_pipeline_dynamical_latentime.py -s ${saras[$index]} --name ${labels[$index]}
      #break
  done

  #cd ../results/gsea/lisa/
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

  # run through sara's pipeline
  # Rscript seurat_sara_pipeline.R --seuratobj "../data/cellranger_counts/20191031_MSK74711_cellranger/outs/filtered_feature_bc_matrix"  \
  #                                --mixtureobj "../data/cellranger_counts/20191031_MSK74711_cellranger_mixture/outs/filtered_feature_bc_matrix" --label 20191031_MSK74711

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

  echo 'Lane,Sample,Index' > MAST118_final.csv
  echo '1-4,MAST118,SI-GA-B6' >> MAST118_final.csv
  cellranger mkfastq --id=MAST118 \
                    --localcores=8 \
                    --localmem=64 \
                    --run=/data/langenau/human_rms_pdxs/20191203_MAST118_FINAL \
                    --csv=MAST118_final.csv


   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/MAST139_1cells/MAST139_1cells/outs/fastq_path/HL73YBGXC/MAST139_1cell/ Tumor24 MAST139_1cells /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/"
   #bash ../../rms_codes/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/MAST139_1cells/MAST139_1cells/outs/fastq_path/HL73YBGXC/MAST139_1cell/ Tumor24 MAST139_1cells /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/
   #bsub -J humanvelR -n 8 -q big-multi -M 64000 "bash ../../rms_codes/cellranger.sh /data/langenau/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/MAST139_1cells_fastq/MAST139_1cells/outs/fastq_path/HL73YBGXC/MAST139_1cell/ Tumor24 MAST139_1cells_mixture /data/molpath/software/10x/refdata-cellranger-hg19-and-mm10-3.0.0/"

   #Rscript seurat_sara_pipeline.R --seuratobj "../data/cellranger_counts/MAST139_1cells/outs/filtered_feature_bc_matrix/" \
   #                               --mixtureobj "../data/cellranger_counts/MAST139_1cells_mixture/outs/filtered_feature_bc_matrix/" --label MAST139_1cells

   #for i in MAST139_1cells; do
   ##for i in 20191031_MSK72117tencell_cellranger; do
   ##for i in 20191031_MSK74711_cellranger; do
   #   velocyto run10x /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i} /data/langenau/alvin_singlecell/01_rms_projects/02_human/data/genome/refdata-cellranger-hg19-3.0.0/genes/genes.gtf # 2>&1 > ${i}.logs &
   #done 

   #for i in 20191031_MSK72117tencell_cellranger 20191031_MSK74711_cellranger MAST139_1cells; do
   #for i in 20191031_MSK74711_cellranger; do
   #    Rscript seurat_pipeline.R --seuratobj /PHShome/qq06/alvin_singlecell/01_rms_projects/01_fish/data/cellranger_counts/${i}/velocyto/* --label $i --finalres 0.8 --tumor -1 --assaytype spliced --species human 
   #done
   #done
}

main

