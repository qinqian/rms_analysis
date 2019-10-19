#!/usr/bin/env python


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import scrublet as scr
import scipy.io
import pandas as pd
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-mat')
    parser.add_argument('-feature')
    parser.add_argument('-name')
    parser.add_argument('-gzip', action='store_true', default=False)
    args = parser.parse_args()

    counts_matrix = scipy.io.mmread(args.mat).T.tocsc()
    if args.gzip:
       genes = np.array(pd.read_table(args.feature, delimiter='\t', compression='gzip').iloc[:, 1])
    else:
       genes = np.array(pd.read_table(args.feature, delimiter='\t').iloc[:, 1])
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    #scrub.call_doublets(threshold=0.25)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                              min_cells=3, 
                                                              min_gene_variability_pctl=85, 
                                                              n_prin_comps=30)
    
    print(doublet_scores)
    scrub.plot_histogram()
    plt.savefig('../results/%s_doublet.pdf' % args.name)

    print('run umap')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.')
    scrub.plot_embedding('UMAP', order_points=True);
    plt.savefig('../results/%s_scrub_umap.pdf' % args.name)

