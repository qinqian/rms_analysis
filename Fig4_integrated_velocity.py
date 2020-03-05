# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import loompy
import scanpy as sc
from glob import glob 
import os
import sys
import os
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import matplotlib.pyplot as plt
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import seaborn as sns
loompy.__version__
seurat = importr("Seurat")
scv.settings.set_figure_params('scvelo')
scv.settings.verbosity = 3
sc.logging.print_versions()
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
pandas2ri.activate()
import matplotlib
print(matplotlib.__version__)
from scanpy.plotting import _utils
from scanpy.plotting._tools.paga import paga
from scanpy.plotting._tools.scatterplots import _get_data_points

from anndata import AnnData

## https://romanhaa.github.io/blog/paga_to_r/
def paga_compare(
    adata: AnnData,
    ax=None,
    basis=None,
    edges=False,
    color=None,
    alpha=None,
    groups=None,
    components=None,
    projection='2d',
    legend_loc='on data',
    legend_fontsize=None,
    legend_fontweight='bold',
    legend_fontoutline=None,
    color_map=None,
    palette=None,
    frameon=False,
    size=None,
    title=None,
    right_margin=None,
    left_margin=0.05,
    show=None,
    save=None,
    title_graph=None,
    groups_graph=None,
    **paga_graph_params,
):
    """\
    Scatter and PAGA graph side-by-side.
    Consists in a scatter plot and the abstracted graph. See
    :func:`~scanpy.pl.paga` for all related parameters.
    See :func:`~scanpy.pl.paga_path` for visualizing gene changes along paths
    through the abstracted graph.
    Additional parameters are as follows.
    Parameters
    ----------
    adata
        Annotated data matrix.
    kwds_scatter
        Keywords for :func:`~scanpy.pl.scatter`.
    kwds_paga
        Keywords for :func:`~scanpy.pl.paga`.
    Returns
    -------
    A list of :class:`~matplotlib.axes.Axes` if `show` is `False`.
    """
    if color is None:
        color = adata.uns['paga']['groups']
    suptitle = None  # common title for entire figure
    if title_graph is None:
        suptitle = color if title is None else title
        title, title_graph = '', ''
    if basis is None:
        if 'X_draw_graph_fa' in adata.obsm.keys():
            basis = 'draw_graph_fa'
        elif 'X_umap' in adata.obsm.keys():
            basis = 'umap'
        elif 'X_tsne' in adata.obsm.keys():
            basis = 'tsne'
        elif 'X_draw_graph_fr' in adata.obsm.keys():
            basis = 'draw_graph_fr'
        else:
            basis = 'umap'

    if 'labels' in paga_graph_params:
        labels = paga_graph_params.pop('labels')
    else:
        labels = groups_graph
    if legend_fontsize is not None:
        paga_graph_params['fontsize'] = legend_fontsize
    if legend_fontweight is not None:
        paga_graph_params['fontweight'] = legend_fontweight
    if legend_fontoutline is not None:
        paga_graph_params['fontoutline'] = legend_fontoutline
    sc.pl.umap(adata, ax=ax, show=False)
#     if 'pos' not in paga_graph_params:
#         if color == adata.uns['paga']['groups']:
#             paga_graph_params['pos'] = _utils._tmp_cluster_pos
#         else:
#             paga_graph_params['pos'] = adata.uns['paga']['pos'] 
    categories = list(adata.obs['clusters'].cat.categories)
    all_pos = np.zeros((len(categories), 2))
    datapoints, components = _get_data_points(adata, 'umap', '2d', 'all')
    for ilabel, label in enumerate(categories):
        _scatter = datapoints[0][adata.obs['clusters'] == label, :]
        x_pos, y_pos = np.median(_scatter, axis=0)
        all_pos[ilabel] = [x_pos, y_pos]
#         if mast111.obs.index[0].split('_')[1] == 'MSK82489':
#             if label == 'G2M':
#                 all_pos[ilabel] = [all_pos[ilabel][0]*4.25, all_pos[ilabel][1]*1.2]    
    paga_graph_params['pos'] = all_pos
    
    paga(
        adata,
        ax=ax,
        show=False,
        save=False,
        title=title_graph,
        labels=labels,
        colors=color,
        frameon=frameon,
        node_size_scale=3, 
        **paga_graph_params,
    )
    _utils.savefig_or_show('paga_compare', show=show, save=save)
    if show == False: return ax

# +
label_stat = {}
fig, axes = plt.subplots(3, 5)
fig.set_size_inches(38, 24)
col = 0
row = 0
degrees = []
adatas = []

# for g in ['../results/MAST35_velocity.h5ad', 
#           '../results/MAST85_velocity.h5ad', 
#           '../results/MAST139_velocity.h5ad',
#           '../results/MAST95_velocity.h5ad', 
#           '../results/MSK82489_velocity.h5ad']:
labels = []
for g in ['../results/velocity_dynamical/MAST35_velocity.h5ad', 
          '../results/velocity_dynamical/MAST85_velocity.h5ad', 
          '../results/velocity_dynamical/MAST139_velocity.h5ad',
          '../results/velocity_dynamical/MAST95_velocity.h5ad', 
          '../results/velocity_dynamical/MSK82489_velocity.h5ad']:
    mast111 = sc.read(g)
    print(g)
    label = os.path.basename(g).split('_')[0]
    labels.append(label)
    seurat = glob(f'../results/seurat_intersect_velocity/*{label}*_seu.rds')
    test_seu = r('readRDS')(seurat[0])
    with localconverter(ro.default_converter + pandas2ri.converter):
        meta = ro.conversion.rpy2py(test_seu.slots['meta.data'])    
#     meta = pandas2ri.ri2py(test_seu.slots['meta.data'])

#    selection = np.load(f'../results/{label}_velocity.npy')

    selection = np.load(f'../results/velocity_dynamical/{label}_velocity.npy')
    clusters = meta.loc[selection, 'RNA_snn_res.0.8']
    metalabels = np.array(["GROUND", "Hypoxia", "EMT", "G1S", "UNASSIGNED", "G2M", 
                           "MUSCLE", "INTERFERON", "PROLIF", "Histones"])
    metacolors = np.array(["#8D510B", "#F19545", "#672366", "#3465FC", "#F2F2F2",
                           "#3465FC", "#E93F33", "#418107", "#3465FC", "#F769A1"])
    
    colortab = pd.read_table('color_table.xls')    
    color_dict = {}
    color_dict2 = {}
    metalabels_dict = {}
    for i, m in enumerate(metalabels):
        metalabels_dict[m] = i
    colors = np.array([ metalabels_dict[m] for m in colortab.loc[:, label].dropna() ])
    for i, j, c in zip(range(len(colortab.loc[:, label].dropna())), 
                       colortab.loc[:, label].dropna(), 
                       metacolors[colors]):
        print(i, j, c)
        color_dict[str(i)] = j
        color_dict2[j] = c
    mast111.obs['clusters'] = clusters.values.map(color_dict)
    mast111.obs['label'] = clusters.values.map(color_dict)
    sc.tl.paga(mast111)
    sc.tl.paga(mast111, groups='clusters', use_rna_velocity=False) ## use_rna_velocity=True is buggy and unstable...!!!!!

    mast111.uns['clusters_colors'] = [color_dict2[c] for c in mast111.obs['clusters'].cat.categories]
    
    print(mast111.obs['clusters'].cat.categories)
    if col == 5:
        row = 1
        col = 0
#         break

    top_genes = mast111.var_names[mast111.var.fit_likelihood.argsort()[::-1]][:300]
    scv.pl.heatmap(mast111, var_names=top_genes, tkey='latent_time', n_convolve=100,
                   col_color='label')
    scv.pl.scatter(mast111, basis=top_genes[:10], legend_loc='none',
               size=80, frameon=False, ncols=5, fontsize=20)
    scv.pl.scatter(mast111, x='latent_time', y=top_genes[:10], 
                   fontsize=16, size=100,
                   n_convolve=100, frameon=False, legend_loc='none')

    scv.pl.velocity_embedding_stream(mast111, basis='umap', color=['clusters'], title=label,
                                     legend_fontsize=25, alpha=.6, show=False,
                                     legend_loc='none',
                                     linewidth=2,
                                     ax=axes[row, col])
    scv.pl.velocity_embedding_stream(mast111, basis='umap', color=['clusters'], title=label,
                                     legend_fontsize=25, alpha=.6, show=True,
                                     linewidth=2)    
    scv.pl.scatter(mast111, basis='umap', color=['latent_time'], title=label,
                   size=100, legend_fontsize=40, fontsize=40,        
                   vmin=0, vmax=1.0, 
                   show=False, ax=axes[row+1, col])
    paga_compare(mast111, title='', legend_loc='none', threshold=0.1, 
                 right_margin=0.2, size=20, edge_width_scale=1, basis='umap',
                 legend_fontsize=25, frameon=False, show=False, ax=axes[row+2, col])
    col += 1
    
    adatas.append(mast111)
    
for a in axes:
    for ax in a:
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                      ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(25)
# -

fig

fig.savefig("Fig4_latent_time_examples.png")

for i in range(len(adatas)):
    plt.figure()
    sns.violinplot(x=adatas[i].obs.loc[:, 'clusters'], 
                   y=adatas[i].obs.loc[:, 'latent_time'])
    plt.title(label=labels[i])

# +
label_stat = {}
fig, axes = plt.subplots(3, 5)
fig.set_size_inches(38, 24)
col = 0
row = 0
degrees = []
adatas = []
labels = []

# for g in ['../results/velocity_dynamical/MAST111_velocity.h5ad', 
#           '../results/velocity_dynamical/RH74_velocity.h5ad']:
# for g in ['../results/MAST35_velocity.h5ad', 
#           '../results/MAST85_velocity.h5ad', 
#           '../results/MAST139_velocity.h5ad',
#           '../results/MAST95_velocity.h5ad', 
#           '../results/MSK82489_velocity.h5ad']:
for g in ['../results/velocity_dynamical/MAST35_velocity.h5ad', 
          '../results/velocity_dynamical/MAST85_velocity.h5ad', 
          '../results/velocity_dynamical/MAST139_velocity.h5ad']:
    mast111 = sc.read(g)
    print(g)
    label = os.path.basename(g).split('_')[0]
    labels.append(label)
    seurat = glob(f'../results/seurat_intersect_velocity/*{label}*_seu.rds')
    test_seu = r('readRDS')(seurat[0])
    with localconverter(ro.default_converter + pandas2ri.converter):
        meta = ro.conversion.rpy2py(test_seu.slots['meta.data'])    
#     meta = pandas2ri.ri2py(test_seu.slots['meta.data'])

#    selection = np.load(f'../results/{label}_velocity.npy')

    selection = np.load(f'../results/velocity_dynamical/{label}_velocity.npy')
    clusters = meta.loc[selection, 'RNA_snn_res.0.8']
    metalabels = np.array(["GROUND", "Hypoxia", "EMT", "G1S", "UNASSIGNED", "G2M", 
                           "MUSCLE", "INTERFERON", "PROLIF", "Histones"])
    metacolors = np.array(["#8D510B", "#F19545", "#672366", "#3465FC", "#F2F2F2",
                           "#3465FC", "#E93F33", "#418107", "#3465FC", "#F769A1"])
    
    colortab = pd.read_table('color_table.xls')    
    color_dict = {}
    color_dict2 = {}
    metalabels_dict = {}
    for i, m in enumerate(metalabels):
        metalabels_dict[m] = i
    colors = np.array([ metalabels_dict[m] for m in colortab.loc[:, label].dropna() ])
    for i, j, c in zip(range(len(colortab.loc[:, label].dropna())), 
                       colortab.loc[:, label].dropna(), 
                       metacolors[colors]):
        print(i, j, c)
        color_dict[str(i)] = j
        color_dict2[j] = c
    mast111.obs['clusters'] = clusters.values.map(color_dict)
    mast111.obs['label'] = clusters.values.map(color_dict)
    sc.tl.paga(mast111)
    sc.tl.paga(mast111, groups='clusters', use_rna_velocity=False) ## use_rna_velocity=True is buggy and unstable...!!!!!

    mast111.uns['clusters_colors'] = [color_dict2[c] for c in mast111.obs['clusters'].cat.categories]
    
    print(mast111.obs['clusters'].cat.categories)
    if col == 5:
        row = 1
        col = 0
#         break
    scv.pl.velocity_embedding_stream(mast111, basis='umap', color=['clusters'], title=label,
                                     legend_fontsize=25, alpha=.6, show=False,
                                     legend_loc='none',
                                     linewidth=2,
                                     ax=axes[row, col])
    scv.pl.velocity_embedding_stream(mast111, basis='umap', color=['clusters'], title=label,
                                     legend_fontsize=25, alpha=.6, show=True,
                                     linewidth=2)    
    scv.pl.scatter(mast111, basis='umap', color=['latent_time'], title=label,
                   size=100, legend_fontsize=40, fontsize=40,        
                   vmin=0, vmax=1.0, 
                   show=False, ax=axes[row+1, col])
    paga_compare(mast111, title='', legend_loc='none', threshold=0.1, 
                 right_margin=0.2, size=20, edge_width_scale=1, basis='umap',
                 legend_fontsize=25, frameon=False, show=False, ax=axes[row+2, col])
    col += 1
    
    adatas.append(mast111)
    
for a in axes:
    for ax in a:
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                      ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(25)
# -

fig

adatas[0].obs.head()

scv.pl.velocity_embedding(adatas[0], color=['root_cells', 'end_points'])

for i in range(len(adatas)):
    plt.figure()
    sns.violinplot(x=adatas[i].obs.loc[:, 'clusters'], 
                   y=adatas[i].obs.loc[:, 'latent_time'], title=labels[i])
    plt.title(labels[i])

# # possibly mesenchymal cell types

for a in adatas:
    scv.pl.velocity_embedding_stream(a, color=['VIM', 'latent_time'])
    scv.pl.velocity_embedding_stream(adatas[2], color=['ZEB1', 'latent_time'])

scv.pl.velocity_embedding_stream(adatas[2], color=['TGFB1', 'latent_time'])

adatas[2].var.index[adatas[2].var.index.str.startswith('SN')]

from glob import glob
looms = glob('../results/seurat_intersect_velocity/*loom')

looms

# +
# %reload_ext autoreload
# %autoreload 2

from resvel import Projector, SCIntegration
rms_int = SCIntegration(adatas, method='conos', strategy=1, seed=99)
rms_int = rms_int.get_integration()
# -

rms_int.obs.loc[:, 'samples'] = rms_int.obs.index.map(lambda x: x.split('_')[1])

sc.pl.umap(rms_int, color='samples')

rms_int_copy = rms_int.copy()

# +
# rms_int_copy = rms_int_copy[rms_int_copy.obs.percent_mt < 10, ]
# -

rms_int_copy.shape

# +
from resvel import run_velocity

rms_int_run = run_velocity(rms_int_copy, mode='dynamical')
# -

rms_int_run.obs.columns

# +
fig, ax = plt.subplots(1, 2)
fig.set_size_inches(18, 9)

scv.pl.velocity_embedding_stream(rms_int_run, ax=ax[0])
scv.pl.velocity_embedding_stream(rms_int_run, color='latent_time', ax=ax[1])
plt.savefig('Integrated_ERMS.png')
# -

scv.pl.velocity_embedding_grid(rms_int_run)

print(rms_int_run.var['velocity_genes'].sum(), rms_int_run.n_vars)
top_genes = rms_int_run.var_names[rms_int_run.var.fit_likelihood.argsort()[::-1]]

order = rms_int_run.obs.groupby(by=["clusters"])["latent_time"].median()
order.index.values[np.argsort(order.values)]
sns.violinplot(x=rms_int_run.obs.loc[:, 'clusters'], 
               y=rms_int_run.obs.loc[:, 'latent_time'], title='Integrated ERMS',
               order=order.index.values[np.argsort(order.values)])

rms_int_run.var.columns = ['gene'] + list(rms_int_run.var.columns[1:])
top_genes = rms_int_run.var_names[rms_int_run.var.fit_likelihood.argsort()[::-1]][:300]
scv.pl.heatmap(rms_int_run, var_names=top_genes, tkey='latent_time', n_convolve=100,
               col_color='label')
scv.pl.scatter(rms_int_run, basis=top_genes[:10], legend_loc='none',
           size=80, frameon=False, ncols=5, fontsize=20)
scv.pl.scatter(rms_int_run, x='latent_time', y=top_genes[:10], 
               fontsize=16, size=100,
               n_convolve=100, frameon=False, legend_loc='none')

scv.pl.scatter(rms_int_run, x='latent_time', y=top_genes[:10], 
               fontsize=16, size=100, color='samples',
               n_convolve=100, frameon=False, legend_loc='none')

# +
# %reload_ext autoreload
# %autoreload 2

label_stat = {}
fig, axes = plt.subplots(3, 5)
fig.set_size_inches(38, 24)
col = 0
row = 0
degrees = []

adatas = []
labels = []

# for g in ['../results/velocity_dynamical/MAST111_velocity.h5ad', 
#           '../results/velocity_dynamical/RH74_velocity.h5ad']:
# for g in ['../results/MAST35_velocity.h5ad', 
#           '../results/MAST85_velocity.h5ad', 
#           '../results/MAST139_velocity.h5ad',
#           '../results/MAST95_velocity.h5ad', 
#           '../results/MSK82489_velocity.h5ad']:
for g in ['../results/velocity_dynamical/MAST85_velocity.h5ad', 
          '../results/velocity_dynamical/MSK82489_velocity.h5ad']:
    mast111 = sc.read(g)
    print(g)
    label = os.path.basename(g).split('_')[0]
    labels.append(label)
    seurat = glob(f'../results/seurat_intersect_velocity/*{label}*_seu.rds')
    test_seu = r('readRDS')(seurat[0])
    with localconverter(ro.default_converter + pandas2ri.converter):
        meta = ro.conversion.rpy2py(test_seu.slots['meta.data'])    
#     meta = pandas2ri.ri2py(test_seu.slots['meta.data'])

#    selection = np.load(f'../results/{label}_velocity.npy')

    selection = np.load(f'../results/velocity_dynamical/{label}_velocity.npy')
    clusters = meta.loc[selection, 'RNA_snn_res.0.8']
    metalabels = np.array(["GROUND", "Hypoxia", "EMT", "G1S", "UNASSIGNED", "G2M", 
                           "MUSCLE", "INTERFERON", "PROLIF", "Histones"])
    metacolors = np.array(["#8D510B", "#F19545", "#672366", "#3465FC", "#F2F2F2",
                           "#3465FC", "#E93F33", "#418107", "#3465FC", "#F769A1"])
    
    colortab = pd.read_table('color_table.xls')    
    color_dict = {}
    color_dict2 = {}
    metalabels_dict = {}
    for i, m in enumerate(metalabels):
        metalabels_dict[m] = i
    colors = np.array([ metalabels_dict[m] for m in colortab.loc[:, label].dropna() ])
    for i, j, c in zip(range(len(colortab.loc[:, label].dropna())), 
                       colortab.loc[:, label].dropna(), 
                       metacolors[colors]):
        print(i, j, c)
        color_dict[str(i)] = j
        color_dict2[j] = c
    mast111.obs['clusters'] = clusters.values.map(color_dict)
    mast111.obs['label'] = clusters.values.map(color_dict)
    sc.tl.paga(mast111)
    sc.tl.paga(mast111, groups='clusters', use_rna_velocity=False) ## use_rna_velocity=True is buggy and unstable...!!!!!

    mast111.uns['clusters_colors'] = [color_dict2[c] for c in mast111.obs['clusters'].cat.categories]
    
    print(mast111.obs['clusters'].cat.categories)
    if col == 5:
        row = 1
        col = 0
#         break
    scv.pl.velocity_embedding_stream(mast111, basis='umap', color=['clusters'], title=label,
                                     legend_fontsize=25, alpha=.6, show=False,
                                     legend_loc='none',
                                     linewidth=2,
                                     ax=axes[row, col])
    scv.pl.velocity_embedding_stream(mast111, basis='umap', color=['clusters'], title=label,
                                     legend_fontsize=25, alpha=.6, show=True,
                                     linewidth=2)    
    scv.pl.scatter(mast111, basis='umap', color=['latent_time'], title=label,
                   size=100, legend_fontsize=40, fontsize=40,        
                   vmin=0, vmax=1.0, 
                   show=False, ax=axes[row+1, col])
    paga_compare(mast111, title='', legend_loc='none', threshold=0.1, 
                 right_margin=0.2, size=20, edge_width_scale=1, basis='umap',
                 legend_fontsize=25, frameon=False, show=False, ax=axes[row+2, col])
    col += 1
    
    adatas.append(mast111)
    
for a in axes:
    for ax in a:
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                      ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(25)
            
from resvel import Projector, SCIntegration
arms_int = SCIntegration(adatas, method='conos', strategy=1, seed=99)
arms_int = arms_int.get_integration()

# +
from resvel import run_velocity

arms_int_run = run_velocity(arms_int, mode='dynamical')

# +
fig, ax = plt.subplots(1, 2)
fig.set_size_inches(18, 9)

scv.pl.velocity_embedding_stream(arms_int_run, ax=ax[0])
scv.pl.velocity_embedding_stream(arms_int_run, color='latent_time', ax=ax[1])
plt.savefig('Integrated_ARMS.png')
# -
order = arms_int_run.obs.groupby(by=["clusters"])["latent_time"].median()
order.index.values[np.argsort(order.values)]
sns.violinplot(x=arms_int_run.obs.loc[:, 'clusters'], 
               y=arms_int_run.obs.loc[:, 'latent_time'], title='Integrated ARMS',
               order=order.index.values[np.argsort(order.values)])


arms_int_run.var.columns = ['gene'] + list(arms_int_run.var.columns[1:])
atop_genes = arms_int_run.var_names[arms_int_run.var.fit_likelihood.argsort()[::-1]][:300]
scv.pl.heatmap(arms_int_run, var_names=atop_genes, tkey='latent_time', n_convolve=100,
               col_color='label')
scv.pl.scatter(arms_int_run, basis=atop_genes[:10], legend_loc='none',
           size=80, frameon=False, ncols=5, fontsize=20)
scv.pl.scatter(arms_int_run, x='latent_time', y=atop_genes[:10], 
               fontsize=16, size=100,
               n_convolve=100, frameon=False, legend_loc='none')


