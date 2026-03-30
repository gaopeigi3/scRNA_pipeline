import scanpy as sc
import scanpy.external as sce
import anndata as ad

from src.qc import run_qc
from src.preprocess import run_preprocess
from src.reduce import run_umap
from src.robustness import check_10x,check_marker_coverage
from src.markers import hierarchical_markers 

def run_multi_sample(samples, params):
    adatas = []

    for sample, path in samples.items():
        print(f"[LOAD] {sample}")
        check_10x(path)
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)

        # 🔥 metadata（必须）
        adata.obs["sample"] = sample

        # QC（per sample）
        adata = run_qc(
            adata,
            min_genes=params["qc"]["min_genes"],
            max_mt=params["qc"]["max_mt"]
        )
        # 🔥 保证 clean index
        adata.obs_names_make_unique()

        adatas.append(adata)

    # merge
    print("[MERGE]")
    # adata = adatas[0].concatenate(*adatas[1:], batch_key="sample")

    adata = ad.concat(
        adatas,
        label="sample",
        keys=list(samples.keys())   # 🔥关键
    )
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    # preprocess
    adata = run_preprocess(
        adata,
        n_hvg=params["preprocess"]["n_hvg"]
    )
    missing = check_marker_coverage(adata, hierarchical_markers)
    if missing:
        print(f"⚠️ Adding {len(missing)} missing markers into HVG")
        print(missing[:10], "...")

        adata.var["highly_variable"] |= adata.var_names.isin(missing)
    # 🔥 Integration（核心）
    method = params["integrate"]["method"]
    key = params["integrate"]["key"]

    use_rep = None

    if method == "harmony":
        print("[INTEGRATION] Harmony")

        sc.pp.scale(adata, max_value=10)
        # sc.tl.pca(adata, svd_solver="arpack", n_comps=50) # pca用全部基因
        sc.tl.pca(adata, use_highly_variable=True) # pca用HVG基因
        # print("X_pca shape:", adata.obsm["X_pca"].shape)
        # print("X_pca:", type(adata.obsm["X_pca"]))
        # print(adata.obs["sample"].value_counts())
        # print(adata.obs["sample"].unique())
        hvg = adata.var["highly_variable"]
        ery_genes = ["HBA1", "HBA2", "HBD", "ALAS2"]
        print(adata.var.loc[ery_genes, "highly_variable"])
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