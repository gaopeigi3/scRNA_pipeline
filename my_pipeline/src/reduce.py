import scanpy as sc


def run_umap(adata, n_neighbors=10, resolution=0.5):
    print("[UMAP]")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=40)
    sc.tl.umap(adata)

    sc.tl.leiden(adata, resolution=resolution)

    return adata