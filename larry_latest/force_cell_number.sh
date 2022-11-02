
# https://kb.10xgenomics.com/hc/en-us/articles/115003331692-Can-you-repeat-cell-calling-without-re-running-cellranger-count-
export PATH=/PHShome/qq06/ssd/alvin/larry/cellranger-7.0.1/:${PATH}

#cellranger reanalyze --id day0_seq1_GRCh38_spliced_reanalyze --force-cells 10000 --matrix /srv/local/alvin/larry/day0_seq1_GRCh38_spliced/outs/filtered_feature_bc_matrix.h5

cellranger reanalyze --id day0_seq1_GRCh38_reanalyze --force-cells 10000 --matrix /srv/local/alvin/larry/day0_seq1_GRCh38/outs/filtered_feature_bc_matrix.h5



