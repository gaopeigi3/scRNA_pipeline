# src/annotation.py

# =========================
# 🎨 Cell type colors
# =========================
celltype_colors_dict = {
    # ---- CD4 T ----
    "CD4": "#8dd3c7",
    "CD4 Memory T cells": "#66c2a5",
    "CD4 Naive T cells": "#99d8c9",
    "CD4 CTL T cells": "#41ae76",
    "CD4 Exhausted T cells": "#238b45",
    "CD4 Th1 T cells": "#005824",
    "CD4 Th2 T cells": "#b2e2e2",
    "CD4 Th17 T cells": "#7bccc4",
    "CD4 Tfh T cells": "#2ca25f",
    "Treg": "#006d2c",

    # ---- CD8 T ----
    "CD8": "#fb8072",
    "CD8 Naive T cells": "#fdbb84",
    "CD8 Effector Memory T cells": "#e34a33",
    "CD8 Exhausted T cells": "#b30000",
    "CD8 CTL T cells": "#f16913",
    "MAIT": "#fb6a4a",

    # ---- B cells ----
    "B": "#80b1d3",
    "B intermediate": "#4eb3d3",
    "B memory": "#2b8cbe",
    "B naive": "#7fcdbb",
    "Plasmablast": "#1c9099",

    # ---- pre B ----
    "pre B": "#08306b",

    # ---- Monocytes ----
    "Mono": "#bebada",
    "CD14 Mono": "#9e9ac8",
    "CD16 Mono": "#756bb1",

    # ---- NK ----
    "NK": "#fccde5",
    "NK CD56-dim": "#fa9fb5",
    "NK Proliferating": "#dd3497",
    "NK CD56-bright": "#980043",

    # ---- EMPs ----
    # ---- EMPs ----
    "EMPs": "#fdb462", 
    "Megakaryocyte": "#f47c3c",
    "Erythroid progenitor": "#fed98e",

    # ---- DC ----
    "DC": "#b3de69",
    "ASDC": "#66c2a5",
    "cDC1": "#31a354",
    "cDC2": "#006d2c",
    "pDC": "#b2df8a",

    # ---- Progenitors / Stem cells ----
    "HSC": "#bc80bd",
    "MonocyticLineage": "#9e5fa8",   # 更深紫（偏分化/activated）
    "Antigen": "#d9a8db",            # 更浅紫（更primitive/接近HSC）
    "CLP": "#8c6bb1",
    "GMP": "#88419d",

    # ---- Macrophage ----
    "Macrophage": "#ffed6f",

    # ---- Erythroid ----
    "Erythroid": "#d9d9d9",
    "Early Erythroid": "#bdbdbd",
    "Late Erythroid": "#969696",

    # "Erythroid": "#f0f0f0",                 # general / pan-erythroid (very light gray)

    # "BFU_E": "#e0ecf4",                      # very early progenitor (pale blue-gray)
    # "CFU_E": "#bfd3e6",                      # committed progenitor (light blue)

    # "Proerythroblast_ProE": "#9ebcda",       # start terminal differentiation
    # "Basophilic_Erythroblast": "#8c96c6",    # early terminal
    # "Polychromatic_Erythroblast": "#8c6bb1", # mid terminal (more purple)
    # "Orthochromatic_Erythroblast": "#88419d",# late terminal (deep purple)

    # "Reticulocyte": "#810f7c",               # enucleated but immature (dark purple)
    # "RBC_Mature_Erythrocyte": "#4d004b",     # mature RBC (very dark purple)

    # ---- Platelets ----
    "Platelets": "#a6cee3",

    # ---- Default ----
    "Unknown": "#999999"
}

# =========================
# 🧬 Marker system
# =========================

import numpy as np
import pandas as pd
import scanpy as sc
from typing import Dict, Tuple

def compute_cluster_means(adata) -> pd.DataFrame:
    """
    高效计算每个 cluster 的平均表达（避免 to_df 内存爆炸）
    返回: DataFrame (clusters x genes)
    """
    clusters = adata.obs["leiden"].unique()

    means = {}
    for cluster in clusters:
        idx = adata.obs["leiden"] == cluster
        # mean over cells → (1, genes)
        means[cluster] = np.asarray(adata[idx].X.mean(axis=0)).ravel()

    cluster_means = pd.DataFrame(means, index=adata.var_names).T

    # z-score normalization（带稳定项）
    cluster_means = (cluster_means - cluster_means.mean()) / (cluster_means.std() + 1e-6)
    cluster_means = cluster_means.fillna(0)

    return cluster_means

def annotate_by_marker_voting(
    adata,
    hierarchical_markers: Dict,
    threshold_main: float = 0.3,
    threshold_sub: float = 0.4
) -> Tuple:

    if "leiden" not in adata.obs:
        raise ValueError("Leiden clustering not found. Run clustering first.")

    """
    基于 hierarchical markers 的 cluster annotation

    Returns:
        adata
        cluster_annotations (dict)
        score_summary (DataFrame)
    """

    cluster_means = compute_cluster_means(adata)

    cluster_annotations = {}
    cluster_scores = {}

    for cluster in cluster_means.index:
        lineage_scores = {}

        # =========================
        # 1️⃣ lineage scoring
        # =========================
        for lineage, info in hierarchical_markers.items():
            genes = [g for g in info["general"] if g in cluster_means.columns]
            if not genes:
                continue

            # 用 median 更稳（比 mean 抗噪）
            lineage_scores[lineage] = cluster_means.loc[cluster, genes].median()

        if not lineage_scores:
            cluster_annotations[cluster] = "Unknown"
            continue

        best_lineage, best_score = max(lineage_scores.items(), key=lambda x: x[1])

        # =========================
        # 2️⃣ main threshold
        # =========================
        if best_score < threshold_main:
            cluster_annotations[cluster] = "Unknown"
            continue

        # =========================
        # 3️⃣ subtype scoring
        # =========================
        subtypes = hierarchical_markers[best_lineage]["subtypes"]
        subtype_scores = {}

        for subtype, genes in subtypes.items():
            genes = [g for g in genes if g in cluster_means.columns]
            if not genes:
                continue

            subtype_scores[subtype] = cluster_means.loc[cluster, genes].median()

        # =========================
        # 4️⃣ final decision
        # =========================
        if subtype_scores:
            best_sub, sub_score = max(subtype_scores.items(), key=lambda x: x[1])

            if sub_score >= threshold_sub:
                final_label = best_sub
            else:
                final_label = best_lineage
        else:
            final_label = best_lineage

        cluster_annotations[cluster] = final_label

        # =========================
        # 5️⃣ 保存score信息（用于debug/论文）
        # =========================
        cluster_scores[cluster] = {
            "best_lineage": best_lineage,
            "best_lineage_score": best_score,
            "lineage_scores": lineage_scores,
            "top_subtypes": sorted(subtype_scores.items(), key=lambda x: x[1], reverse=True)[:3]
        }

    # =========================
    # 6️⃣ 写入 adata
    # =========================
    adata.obs["celltype"] = adata.obs["leiden"].map(cluster_annotations)

    # 防止category报错
    adata.obs["celltype"] = adata.obs["celltype"].astype("category")

    # =========================
    # 7️⃣ summary table（可选但强烈推荐）
    # =========================
    score_summary = pd.DataFrame({
        cluster: {
            **{
                f"main_type_{i+1}": sorted(v["lineage_scores"].items(), key=lambda x: x[1], reverse=True)[i][0]
                if len(v["lineage_scores"]) > i else None
                for i in range(3)
            },
            **{
                f"main_score_{i+1}": sorted(v["lineage_scores"].items(), key=lambda x: x[1], reverse=True)[i][1]
                if len(v["lineage_scores"]) > i else None
                for i in range(3)
            }
        }
        for cluster, v in cluster_scores.items()
    }).T

    return adata, cluster_annotations, score_summary

def apply_celltype_colors(adata, celltype_colors_dict):
    """
    根据 celltype 设置 scanpy 颜色
    """
    if "celltype" not in adata.obs:
        raise ValueError("adata.obs['celltype'] not found")

    adata.obs["celltype"] = adata.obs["celltype"].astype("category")

    adata.uns["celltype_colors"] = [
        celltype_colors_dict.get(ct, "#999999")
        for ct in adata.obs["celltype"].cat.categories
    ]

    return adata

    