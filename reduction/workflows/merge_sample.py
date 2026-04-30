import scanpy as sc
import scanpy.external as sce
import anndata as ad
import harmonypy as hm
import numpy as np
from src.qc import run_qc
from src.preprocess import clean_genes, run_preprocess
from src.reduce import run_umap
from src.robustness import check_10x
# import os
# os.environ["OMP_NUM_THREADS"] = "1"
# os.environ["MKL_NUM_THREADS"] = "1"
# os.environ["OPENBLAS_NUM_THREADS"] = "1"


def run_multi_sample(samples, params):
    adatas = []

    for sample, path in samples.items():
        print(f"[LOAD] {sample}")
        check_10x(path)
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)

        adata.obs["sample"] = sample

        # QC（per sample）
        adata = run_qc(
            adata,
            min_genes=params["qc"]["min_genes"],
            max_mt=params["qc"]["max_mt"]
        )
        adata.obs_names_make_unique()

        adatas.append(adata)

    # merge
    print("[MERGE]")
    # adata = adatas[0].concatenate(*adatas[1:], batch_key="sample")

    adata = ad.concat(
        adatas,
        label="sample",
        keys=list(samples.keys())
    )

    adata = clean_genes(adata, remove_mt=True, remove_ribo=True)

    adata.obs["sample"] = adata.obs["sample"].astype("category")
    # preprocess
    adata = run_preprocess(
        adata,
        n_hvg=params["preprocess"]["n_hvg"]
    )

    # 🔥 Integration（核心）
    method = params["integrate"]["method"]
    key = params["integrate"]["key"]

    use_rep = None

    if method == "harmony":
        print("[INTEGRATION] Harmony")

        # 🔥 必须加这三步
        # adata.raw = adata.copy()
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(
            adata,
            svd_solver="arpack",
            n_comps=50,
            use_highly_variable=True
        )

        # X = adata.obsm["X_pca"]
        # print("PCA shape:", X.shape)
        # print("NaN:", np.isnan(X).sum())
        # print("Inf:", np.isinf(X).sum())
        # print("max:", X.max(), "min:", X.min())
        sce.pp.harmony_integrate(adata, key=key)

        use_rep = "X_pca_harmony"

    elif method == "none":
        print("[INTEGRATION] None")

    else:
        raise ValueError(f"Unknown integration method: {method}")

    # reduce
    adata = run_umap(
        adata,
        n_neighbors=params["reduce"]["n_neighbors"],
        resolution=params["reduce"]["resolution"],
        use_rep=use_rep
    )

    return adata