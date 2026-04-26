Task: per-sample gene UMAP visualization

Input:
  - h5ad (merged.h5ad)
  - gene list
  - sample_key

Output:
  - results/umap_gene/{gene}/{sample}.png


~/projects/scRNA_pipeline/analysis$ PYTHONPATH=. python ./tasks/singlegeneumap.py \
    --config ../config/config_singlegeneumap.yaml