# src/preprocess.py

import scanpy as sc

def run_preprocess(adata, n_hvg=2000):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata = adata[:, adata.var.highly_variable].copy()

    return adata