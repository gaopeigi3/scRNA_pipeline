# src/preprocess.py

import scanpy as sc
def clean_genes(
    adata,
    remove_mt=True,
    remove_ribo=True,
    remove_ensg=False,
    extra_remove=None,
    copy=True,
    verbose=True
):
    import numpy as np

    if copy:
        adata = adata.copy()

    # gene name
    if "symbol" in adata.var.columns:
        gene_names = adata.var["symbol"].astype(str)
    else:
        gene_names = adata.var_names.astype(str)

    gene_upper = gene_names.str.upper()

    remove = np.zeros(len(gene_names), dtype=bool)

    if remove_mt:
        remove |= gene_upper.str.startswith("MT-")

    if remove_ribo:
        remove |= (
            gene_upper.str.startswith("RPS") |
            gene_upper.str.startswith("RPL")
        )

    if remove_ensg:
        remove |= gene_upper.str.startswith("ENSG")

    if extra_remove:
        remove |= gene_upper.isin([g.upper() for g in extra_remove])

    keep = ~remove

    if verbose:
        print(f"[FILTER] Removed genes: {np.sum(remove)}")
        print(f"[FILTER] Remaining genes: {np.sum(keep)}")

    return adata[:, keep].copy()

def run_preprocess(adata, n_hvg=2000):

    # 🔥 1️⃣ 保存 counts（工程关键）
    adata.layers["counts"] = adata.X.copy()

    # 🔥 2️⃣ normalize + log
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # 🔥 3️⃣ 设置 raw（唯一一次）
    adata.raw = adata.copy()

    # HVG
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)

    return adata