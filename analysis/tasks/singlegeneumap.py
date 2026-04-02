import scanpy as sc
import yaml
import argparse
import os
import matplotlib.pyplot as plt


def apply_alias(genes, alias_dict):
    new_genes = []
    for g in genes:
        if g in alias_dict:
            print(f"[INFO] Alias {g} → {alias_dict[g]}")
            new_genes.extend(alias_dict[g])
        else:
            new_genes.append(g)
    return new_genes


def expand_gene_list(genes, var_names, enable_family=True):
    expanded = []
    missing = []

    for g in genes:
        if enable_family:
            matched = [v for v in var_names if v.startswith(g)]
            if len(matched) > 1:
                print(f"[INFO] Expand family {g} → {matched}")
                expanded.extend(matched)
                continue

        if g in var_names:
            expanded.append(g)
        else:
            print(f"[WARN] Missing gene: {g}")
            missing.append(g)

    return list(dict.fromkeys(expanded)), missing

def run(cfg):
    print("[INFO] Loading adata...")
    adata = sc.read_h5ad(cfg["adata"])

    sample_key = cfg["sample_key"]
    outdir = cfg["output_dir"]

    os.makedirs(outdir, exist_ok=True)

    # categorical优化
    adata.obs[sample_key] = adata.obs[sample_key].astype("category")
    samples = adata.obs[sample_key].cat.categories

    # 检查UMAP
    if "X_umap" not in adata.obsm:
        raise ValueError("UMAP not found")

    # 👉 Step1：alias
    genes = apply_alias(cfg["genes"], cfg.get("gene_alias", {}))

    # 👉 Step2：family expand（只做一次！）
    genes, missing_genes = expand_gene_list(
        genes,
        adata.var_names,
        enable_family=cfg.get("gene_family", True)
    )

    print(f"[INFO] Final genes: {genes}")

    # 👉 Step3：记录 missing（只写一次）
    if missing_genes:
        with open(os.path.join(outdir, "missing_genes.txt"), "w") as f:
            for g in missing_genes:
                f.write(g + "\n")

    # 👉 Step4：plot
    for gene in genes:
        gene_dir = os.path.join(outdir, gene)
        os.makedirs(gene_dir, exist_ok=True)

        for s in samples:
            print(f"[INFO] Plot {gene} | {s}")

            adata_sub = adata[adata.obs[sample_key] == s]

            sc.pl.umap(
                adata_sub,
                color=gene,
                show=False,
                cmap=cfg["plot"]["cmap"],
                vmin=cfg["plot"]["vmin"],
                vmax=cfg["plot"]["vmax"],
                title=f"{gene}_{s}" ,
            )

            outfile = os.path.join(gene_dir, f"{gene}_{s}.png")
            plt.savefig(outfile, dpi=cfg["plot"]["dpi"])
            plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    run(cfg)


if __name__ == "__main__":
    main()