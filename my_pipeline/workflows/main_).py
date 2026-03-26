import scanpy as sc
import yaml
import os

from src.reduce import run_umap
from src.visualize import plot_umap
from src.markers import hierarchical_markers
from src.annotation import annotate_by_marker_voting, apply_celltype_colors
from src.annotation import celltype_colors_dict


def run_umap_pipeline(sample, path, params, outdir="data/processed"):
    print(f"=== {sample} ===")

    # load
    adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)

    # QC
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    adata = adata[adata.obs.n_genes_by_counts > 200, :].copy()
    adata = adata[adata.obs.pct_counts_mt < 10, :].copy()

    # preprocess
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable]

    adata = run_umap(
        adata,
        n_neighbors=params["n_neighbors"],
        resolution=params["resolution"]
    )

    # annotation
    adata, cluster_map, summary = annotate_by_marker_voting(
        adata,
        hierarchical_markers,
        threshold_main=params["threshold_main"],
        threshold_sub=params["threshold_sub"]
    )

    os.makedirs("logs", exist_ok=True)
    summary.to_csv(f"logs/{sample}_annotation_summary.csv")

    # color
    adata = apply_celltype_colors(adata, celltype_colors_dict)

    # plot
    plot_umap(adata, sample)

    # save
    os.makedirs(outdir, exist_ok=True)
    adata.write(f"{outdir}/{sample}.h5ad")

    print(f"✅ Done {sample}")


def main():
    config = yaml.safe_load(open("config/config.yaml"))
    params = config["params"]

    for sample, path in config["samples"].items():
        run_umap_pipeline(sample, path, params)


if __name__ == "__main__":
    main()