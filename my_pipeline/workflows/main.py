# import yaml
# from src.load import load_data
# from src.qc import run_qc
# from src.preprocess import preprocess
# from src.reduce import run_umap
# from src.visualize import plot_umap

# def main(config_path="config/config.yaml"):
#     config = yaml.safe_load(open(config_path))

#     for sample, path in config["samples"].items():
#         print(f"=== Processing {sample} ===")

#         adata = load_data(path)
#         adata = run_qc(adata, **config["params"])
#         adata = preprocess(adata)
#         adata = run_umap(adata)
#         adata.write(f"data/processed/{sample}.h5ad")
#         plot_umap(adata)
        
# if __name__ == "__main__":
#     main()
import scanpy as sc
import yaml
import os

def run_umap_pipeline(sample, path, outdir="data/processed"):
    print(f"=== {sample} ===")

    # 1️⃣ load
    adata = sc.read_10x_mtx(
        path,
        var_names="gene_symbols",
        cache=True
    )

    # 2️⃣ QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # 3️⃣ filter（保守一点，别翻车）
    adata = adata[adata.obs.n_genes_by_counts > 200, :]
    adata = adata[adata.obs.pct_counts_mt < 10, :]

    # 4️⃣ normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 5️⃣ HVG
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable]

    # 6️⃣ scale + PCA
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")

    # 7️⃣ neighbors + UMAP
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)

    # 8️⃣ clustering（加分项）
    sc.tl.leiden(adata, resolution=0.5)

    # 9️⃣ 保存
    os.makedirs(outdir, exist_ok=True)

    sc.pl.umap(
        adata,
        color=["leiden"],
        save=f"_{sample}.png",   # 会保存到 figures/
        show=False
    )

    adata.write(f"{outdir}/{sample}.h5ad")

    print(f"✅ Done {sample}")


def main():
    config = yaml.safe_load(open("config/config.yaml"))

    for sample, path in config["samples"].items():
        run_umap_pipeline(sample, path)


if __name__ == "__main__":
    main()