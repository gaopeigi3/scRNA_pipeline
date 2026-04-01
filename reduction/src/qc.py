
# src/qc.py

import scanpy as sc

def run_qc(adata, min_genes=200, max_mt=0.05):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    adata = adata[adata.obs.n_genes_by_counts > min_genes, :].copy()
    adata = adata[adata.obs.pct_counts_mt < max_mt * 100, :].copy()

    return adata