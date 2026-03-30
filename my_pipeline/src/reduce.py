
# src/reduce.py

import scanpy as sc

def run_umap(adata, n_neighbors=15, resolution=0.8, use_rep=None):
    print("[UMAP]")


    if use_rep:
        sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=40)

    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)

    return adata

