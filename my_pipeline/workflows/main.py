# workflows/main.py

import scanpy as sc
import yaml
import os
import argparse
from src.qc import run_qc
from src.preprocess import run_preprocess
from src.reduce import run_umap
# from src.visualize import plot_umap
from src.markers import hierarchical_markers
from src.annotation import annotate_by_marker_voting, apply_celltype_colors
from src.annotation import celltype_colors_dict
import datetime
from src.visualize import plot_umap_all, plot_umap_per_sample


# def run_pipeline(sample, path, params, outdir="data/processed"):
#     print(f"=== {sample} ===")

#     # load
#     adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)

#     # QC
#     adata = run_qc(
#         adata,
#         min_genes=params["qc"]["min_genes"],
#         max_mt=params["qc"]["max_mt"]
#     )

#     # preprocess
#     adata = run_preprocess(
#         adata,
#         n_hvg=params["preprocess"]["n_hvg"]
#     )


#     # reduce
#     adata = run_umap(
#         adata,
#         n_neighbors=params["reduce"]["n_neighbors"],
#         resolution=params["reduce"]["resolution"]
#     )

#     # annotation
#     adata, cluster_map, summary = annotate_by_marker_voting(
#         adata,
#         hierarchical_markers,
#         threshold_main=params["annotation"]["threshold_main"],
#         threshold_sub=params["annotation"]["threshold_sub"]
#     )

#     os.makedirs("logs", exist_ok=True)
#     summary.to_csv(f"logs/{sample}_annotation_summary.csv")

#     # color
#     adata = apply_celltype_colors(adata, celltype_colors_dict)

#     # plot
#     plot_umap(adata, sample)

#     # save
#     os.makedirs(outdir, exist_ok=True)
#     adata.write(f"{outdir}/{sample}.h5ad")

#     print(f"✅ Done {sample}")

from workflows.single_sample import run_single_sample
from workflows.merge_sample import run_multi_sample
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--config",
        type=str,
        default="config/config.yaml",
        help="Path to config file"
    )

    parser.add_argument(
        "--mode",
        type=str,
        choices=["single", "merge"],
        help="Override mode in config"
    )

    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    # 🔥 CLI优先级 > config
    mode = args.mode if args.mode else config.get("mode", "single")
    samples = config["samples"]
    params = config["params"]
    print(f"[MODE] {mode}")
    if mode == "single":
        run_single_sample(samples, params)

    elif mode == "merge":
        adata = run_multi_sample(samples, params)

        adata, cluster_map, summary, score_df = annotate_by_marker_voting(
            adata,
            hierarchical_markers,
            threshold_main=params["annotation"]["threshold_main"],
            threshold_sub=params["annotation"]["threshold_sub"]
        )
        os.makedirs("logs", exist_ok=True)
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M")
        summary.to_csv(f"logs/merged_annotation_{timestamp}.csv")

        adata = apply_celltype_colors(adata, celltype_colors_dict)
        print("Number of clusters:", adata.obs["leiden"].nunique())
        print("HBA1" in adata.var_names)
        print("HBA1" in adata.raw.var_names)
        score_df.to_csv(f"logs/cluster_scores_{timestamp}.csv", index=False)
        os.makedirs("data/processed", exist_ok=True)
        adata.write("data/processed/merged.h5ad")
        # 🔥 merged
        plot_umap_all(adata)

        # 🔥 per sample
        plot_umap_per_sample(adata)
        sc.pl.umap(adata, color="sample", save="_sample.png")
    else:
        raise ValueError("mode must be 'single' or 'merge'")


if __name__ == "__main__":
    main()
