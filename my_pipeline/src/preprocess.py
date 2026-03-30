# src/preprocess.py

import scanpy as sc

def run_preprocess(adata, n_hvg=2000):

    # ✅ 保留 raw counts（极其重要）
    adata.layers["counts"] = adata.X.copy()

    # normalize + log
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ✅ 标记 HVG（但不删基因）
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)

    # ✅ 保存一份 raw（用于 DEG）
    adata.raw = adata

    return adata



