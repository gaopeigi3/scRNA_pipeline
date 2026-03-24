from pathlib import Path
import argparse

# =========================
# 📦 模板定义
# =========================

TEMPLATES = {
    "config/config.yaml": """\
project: scRNA_pipeline

samples:
  sample1: data/raw/sample1/filtered_feature_bc_matrix

paths:
  raw: data/raw
  processed: data/processed

params:
  min_genes: 200
  max_mt: 0.1
""",

    "src/load.py": """\
import scanpy as sc

def load_data(path):
    print(f"[LOAD] {path}")
    adata = sc.read_10x_mtx(
        path,
        var_names="gene_symbols",
        cache=True
    )
    return adata
""",

    "src/qc.py": """\
def run_qc(adata, min_genes=200, max_mt=0.1):
    print("[QC]")
    return adata
""",

    "src/preprocess.py": """\
def preprocess(adata):
    print("[PREPROCESS]")
    return adata
""",

    "src/reduce.py": """\
def run_umap(adata):
    print("[UMAP]")
    return adata
""",

    "src/visualize.py": """\
def plot_umap(adata):
    print("[PLOT]")
""",

    "src/utils.py": """\
def check_10x(path):
    import os
    required = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]
    for f in required:
        if not os.path.exists(os.path.join(path, f)):
            raise FileNotFoundError(f"{f} missing in {path}")
""",

    "workflows/main.py": """\
import yaml
from src.load import load_data
from src.qc import run_qc
from src.preprocess import preprocess
from src.reduce import run_umap
from src.visualize import plot_umap

def main(config_path="config/config.yaml"):
    config = yaml.safe_load(open(config_path))

    for sample, path in config["samples"].items():
        print(f"=== Processing {sample} ===")

        adata = load_data(path)
        adata = run_qc(adata, **config["params"])
        adata = preprocess(adata)
        adata = run_umap(adata)
        plot_umap(adata)

if __name__ == "__main__":
    main()
""",

    "Snakefile": """\
rule all:
    input:
        expand("data/processed/{sample}.h5ad", sample=["sample1"])

rule process:
    input:
        "data/raw/{sample}/filtered_feature_bc_matrix"
    output:
        "data/processed/{sample}.h5ad"
    shell:
        "python workflows/main.py"
"""
}

# =========================
# 📁 目录结构
# =========================

DIRS = [
    "config",
    "data/raw",
    "data/processed",
    "src",
    "workflows",
    "logs"
]


# =========================
# 🚀 创建函数
# =========================

def create_project(base="sc_pipeline", force=False):
    base_path = Path(base)

    if base_path.exists() and not force:
        raise FileExistsError(f"{base} already exists. Use --force to overwrite.")

    print(f"🚀 Creating project: {base}")

    # 创建目录
    for d in DIRS:
        (base_path / d).mkdir(parents=True, exist_ok=True)

    # 写入模板
    for file, content in TEMPLATES.items():
        file_path = base_path / file
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_path.write_text(content)

    # package init
    (base_path / "src/__init__.py").touch()

    print(f"✅ Done: {base_path.resolve()}")


# =========================
# CLI入口
# =========================

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", default="sc_pipeline")
    parser.add_argument("--force", action="store_true")

    args = parser.parse_args()

    create_project(args.name, args.force)