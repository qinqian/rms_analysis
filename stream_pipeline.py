#!/usr/bin/env python
import pickle

import stream as st
import scanpy as sc
import os
from velocity_pipeline import load_seurat_umap
from rpy2.robjects import r, pandas2ri
import matplotlib.pyplot as plt
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd


def load_stream(test, args, run=True):
    if run:
        test[0].obs['label']=test[0].obs['clusters']
        st.add_cell_colors(test[0])
        updated_colors = {}
        for col, cat in zip(test[0].uns['clusters_colors'], test[0].obs['clusters'].cat.categories):
            updated_colors[cat] = col
        test[0].uns['label_color'] = updated_colors
        # adata.obsm['top_pcs'] = adata.obsm['umap_cell_embeddings']

        test[0].obsm['X_dr'] = test[0].obsm['X_pca'] ###velocity_umap
        test[0].uns['workdir'] = '../results/stream'
        os.system('mkdir -p ../results/stream')
        test[0].obsm['X_pca'].shape

        test[0].obsm['X_vis_umap'] = test[0].obsm['X_umap']
        test[0].obsm['top_pcs'] = test[0].obsm['X_pca'] #[:, :3]
        print(test[0].obsm['top_pcs'].shape)

        st.dimension_reduction(test[0], feature='top_pcs', method='se', nb_pct =50.0/test[0].shape[0], n_jobs=30)
        st.seed_elastic_principal_graph(test[0])
        st.elastic_principal_graph(test[0])
        st.write(test[0], f'{args.name}_stream.pkl')

    st.plot_visualization_2D(test[0], fig_legend_ncol=9, save_fig=True, fig_name=f'{args.name}_stream_umap.pdf')
    st.plot_flat_tree(test[0], save_fig=True, fig_name=f'{args.name}_stream_flattree.pdf')
    st.plot_visualization_2D(test[0], fig_legend_ncol=9, save_fig=True, fig_name=f'{args.name}_stream_umap_branch.pdf', color_by='branch')
    ## st.subwaymap_plot(test[0], root='S1', fig_legend_ncol=6) 
    print('stream plot...')
    st.stream_plot(test[0], root='S1',fig_legend_ncol=9, save_fig=True, fig_name=f'{args.name}_stream.pdf')
    print('done')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--name',    default='test', help='output file name')
    parser.add_argument('-s', '--seurat',  default='test', help='seurat object path')
    parser.add_argument('--species',       default='human', help='output file name')
    args = parser.parse_args()

    pandas2ri.activate()
    print(args)
    if args.name == 'test':
        parser.print_help()
        sys.exit(1)

    if os.path.exists(f'../results/stream/{args.name}_stream.pkl'):
        with open(f'../results/stream/{args.name}_stream.pkl', 'rb') as f:
            test = [pickle.loads(f.read())]
            load_stream(test, args, False)
    else:
        test = []
        test.append(sc.read_h5ad(f'../results/{args.name}_velocity.h5ad'))
        test.append(np.load(f'../results/{args.name}_velocity.npy'))
        clusters, reductions, pos = load_seurat_umap(args, test)
        load_stream(test, args, True)
