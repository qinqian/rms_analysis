#!/bin/bash -ex

lisa model --method="all" --web=True --new_rp_h5=None --new_count_h5=None --species hg38 --epigenome "['DNase', 'H3K27ac']" --cluster=False --covariates=False --random=True --prefix gene_modules --background=dynamic_auto_tad --stat_background_number=1000 --threads 10 *symbols
