from src.annotation import apply_celltype_colors
# def plot_umap(adata):
#     if "celltype" in adata.obs:
#         add_celltype_colors(adata, key="celltype")

#         sc.pl.umap(
#             adata,
#             color="celltype",
#             legend_loc="on data",
#             frameon=False,
#             show=False,
#             save="_celltype.png"
#         )
#     else:
#         sc.pl.umap(
#             adata,
#             color="leiden",
#             show=False,
#             save="_leiden.png"
#         )

# def plot_umap(adata):
#     print("[PLOT]")
import scanpy as sc

def plot_umap(adata, sample):
    sc.pl.umap(
        adata,
        color="celltype",
        size=6,
        frameon=False,
        legend_loc="on data",   # 👈 关键
        title=sample,
        show=False,
        save=f"_{sample}_leiden.png"
    )