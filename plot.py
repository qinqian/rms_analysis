import seaborn as sns
import matplotlib.pyplot as plt


def embedding_with_density():
    """ generating a 2x2 grid figure
    """
    fig, axes = plt.subplots(2, 2)
    fig.set_size_inches(15, 12)
    scv.pl.umap(adata, ax=axes[0][0], show=False)
    sns.kdeplot(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1],ax=axes[0][1])
    scv.pl.umap(adata_den, ax=axes[1][0], show=False)
    sns.kdeplot(adata_den.obsm['X_umap'][:, 0], adata_den.obsm['X_umap'][:, 1], ax=axes[1][1])
    return
