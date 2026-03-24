import scanpy as sc

def load_data(path):
    print(f"[LOAD] {path}")
    adata = sc.read_10x_mtx(
        path,
        var_names="gene_symbols",
        cache=True
    )
    return adata
