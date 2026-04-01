# src/annotation.py

# =========================
# 🎨 Cell type colors
# =========================

celltype_colors_dict = {

    # =========================
    # 🔴 T lineage（红系）
    # =========================
    "CD4": "#fb6a4a",
    "CD4 Memory T cells": "#ef3b2c",
    "CD4 Naive T cells": "#fcae91",
    "CD4 CTL T cells": "#cb181d",
    "CD4 Exhausted T cells": "#99000d",
    "CD4 Th1 T cells": "#a50f15",
    "CD4 Th2 T cells": "#fb6a4a",
    "CD4 Th17 T cells": "#de2d26",
    "CD4 Tfh T cells": "#fc9272",
    "Treg": "#67000d",

    "CD8": "#f16913",
    "CD8 Naive T cells": "#fdae6b",
    "CD8 Effector Memory T cells": "#e6550d",
    "CD8 Exhausted T cells": "#a63603",
    "CD8 CTL T cells": "#d94801",
    "MAIT": "#fd8d3c",

    # =========================
    # 🔵 B lineage（蓝系）
    # =========================
    "B": "#3182bd",
    "B intermediate": "#6baed6",
    "B memory": "#08519c",
    "B naive": "#9ecae1",
    "Plasmablast": "#08306b",

    "pre B": "#4292c6",

    # =========================
    # 🌸 NK lineage（粉系）
    # =========================
    "NK": "#f768a1",
    "NK CD56-dim": "#dd3497",
    "NK Proliferating": "#ae017e",
    "NK CD56-bright": "#7a0177",

    # =========================
    # 🟢 Myeloid lineage（绿系）
    # =========================
    "Mono": "#41ae76",
    "CD14 Mono": "#66c2a4",
    "CD16 Mono": "#238b45",

    "GMP": "#74c476",          # 🔥 改成绿（关键）
    "Macrophage": "#006d2c",

    # =========================
    # 🟣 Progenitor / Stem（紫系）
    # =========================
    "HSC": "#9e9ac8",
    "CLP": "#807dba",          # lymphoid progenitor（保留紫）
    "MonocyticLineage": "#6a51a3",
    "Antigen": "#bcbddc",

    # =========================
    # 🟢 DC（橄榄绿，独立于Mono）
    # =========================
    "DC": "#8c9e3f",
    "ASDC": "#b5cf6b",
    "cDC1": "#637939",
    "cDC2": "#4d5d2c",
    "pDC": "#9c9ede",

    # =========================
    # ⚪ Erythroid（灰系）
    # =========================
    "Erythroid": "#bdbdbd",
    "Early Erythroid": "#d9d9d9",
    "Late Erythroid": "#969696",

    # =====（你之前注释掉的完整体系，我也帮你整理好了）=====
    # "BFU_E": "#e0e0e0",
    # "CFU_E": "#cccccc",
    # "Proerythroblast_ProE": "#bdbdbd",
    # "Basophilic_Erythroblast": "#a6a6a6",
    # "Polychromatic_Erythroblast": "#8c8c8c",
    # "Orthochromatic_Erythroblast": "#737373",
    # "Reticulocyte": "#595959",
    # "RBC_Mature_Erythrocyte": "#404040",

    # =========================
    # 🟠 EMP / Megakaryocyte（橙系）
    # =========================
    "EMPs": "#fdae6b",
    "Megakaryocyte": "#e6550d",
    "Erythroid progenitor": "#fdd0a2",

    # =========================
    # 🔷 Platelets
    # =========================
    "Platelets": "#9ecae1",

    # =========================
    # ⚫ Default
    # =========================
    "Unknown": "#969696"
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
# def annotate_by_marker_voting(
#     adata,
#     hierarchical_markers,
#     threshold_main: float = 0.3,
#     threshold_sub: float = 0.4
# ):

#     if "leiden" not in adata.obs:
#         raise ValueError("Leiden clustering not found. Run clustering first.")

#     import numpy as np
#     import pandas as pd

#     # =========================
#     # 🔥 使用 raw（避免HVG丢marker）
#     # =========================
#     use_raw = adata.raw is not None
#     X = adata.raw.X if use_raw else adata.X
#     genes = adata.raw.var_names if use_raw else adata.var_names

#     clusters = adata.obs["leiden"].unique()

#     mean_expr = {}
#     pct_expr = {}

#     for c in clusters:
#         idx = np.where(adata.obs["leiden"].values == c)[0]
#         sub = X[idx]
#         mean_expr[c] = np.asarray(sub.mean(axis=0)).ravel()
#         pct_expr[c] = np.asarray((sub > 0).mean(axis=0)).ravel()

#     mean_df = pd.DataFrame(mean_expr, index=genes).T
#     pct_df = pd.DataFrame(pct_expr, index=genes).T

#     # z-score（按cluster）
#     z_df = (mean_df - mean_df.mean()) / (mean_df.std() + 1e-6)
#     z_df = z_df.fillna(0)

#     # =========================
#     # 🔥 annotation
#     # =========================
#     cluster_annotations = {}
#     cluster_scores = {}

#     for cluster in mean_df.index:

#         lineage_scores = {}

#         # =========================
#         # 1️⃣ lineage scoring（升级版）
#         # =========================
#         for lineage, info in hierarchical_markers.items():
#             genes_list = [g for g in info["general"] if g in mean_df.columns]
#             if not genes_list:
#                 continue

#             mean = mean_df.loc[cluster, genes_list].mean()
#             pct = pct_df.loc[cluster, genes_list].mean()
#             z = z_df.loc[cluster, genes_list].mean()

#             score = 0.5 * mean + 0.3 * pct + 0.2 * z
#             lineage_scores[lineage] = score

#         if not lineage_scores:
#             cluster_annotations[cluster] = "Unknown"
#             continue

#         # 排序
#         sorted_lineages = sorted(lineage_scores.items(), key=lambda x: x[1], reverse=True)
#         best_lineage, best_score = sorted_lineages[0]

#         # =========================
#         # 2️⃣ main threshold
#         # =========================
#         if best_score < threshold_main:
#             cluster_annotations[cluster] = "Unknown"
#             continue

#         # =========================
#         # 3️⃣ subtype scoring（同样升级）
#         # =========================
#         subtypes = hierarchical_markers[best_lineage]["subtypes"]
#         subtype_scores = {}

#         for subtype, genes_list in subtypes.items():
#             genes_list = [g for g in genes_list if g in mean_df.columns]
#             if not genes_list:
#                 continue

#             mean = mean_df.loc[cluster, genes_list].mean()
#             pct = pct_df.loc[cluster, genes_list].mean()
#             z = z_df.loc[cluster, genes_list].mean()

#             subtype_scores[subtype] = 0.5 * mean + 0.3 * pct + 0.2 * z

#         # =========================
#         # 4️⃣ final decision
#         # =========================
#         if subtype_scores:
#             sorted_sub = sorted(subtype_scores.items(), key=lambda x: x[1], reverse=True)
#             best_sub, sub_score = sorted_sub[0]

#             if sub_score >= threshold_sub:
#                 final_label = best_sub
#             else:
#                 final_label = best_lineage
#         else:
#             final_label = best_lineage

#         cluster_annotations[cluster] = final_label

#         # =========================
#         # 🔥 新增：confidence / entropy
#         # =========================
#         scores = np.array([v for _, v in sorted_lineages])
#         probs = np.exp(scores) / np.exp(scores).sum()

#         top1 = scores[0]
#         top2 = scores[1] if len(scores) > 1 else 0

#         confidence = top1 - top2
#         entropy = -np.sum(probs * np.log(probs + 1e-9))

#         # =========================
#         # 5️⃣ 保存score信息（兼容旧结构 + 增强）
#         # =========================
#         cluster_scores[cluster] = {
#             "best_lineage": best_lineage,
#             "best_lineage_score": best_score,
#             "lineage_scores": lineage_scores,
#             "top_subtypes": sorted(subtype_scores.items(), key=lambda x: x[1], reverse=True)[:3],
#             # 🔥新增
#             "confidence": float(confidence),
#             "entropy": float(entropy)
#         }

#     # =========================
#     # 6️⃣ 写入 adata
#     # =========================
#     adata.obs["celltype"] = adata.obs["leiden"].map(cluster_annotations)
#     adata.obs["celltype"] = adata.obs["celltype"].astype("category")

#     # =========================
#     # 7️⃣ summary（保持你原来的格式）
#     # =========================
#     score_summary = pd.DataFrame({
#         cluster: {
#             **{
#                 f"main_type_{i+1}": sorted(v["lineage_scores"].items(), key=lambda x: x[1], reverse=True)[i][0]
#                 if len(v["lineage_scores"]) > i else None
#                 for i in range(3)
#             },
#             **{
#                 f"main_score_{i+1}": sorted(v["lineage_scores"].items(), key=lambda x: x[1], reverse=True)[i][1]
#                 if len(v["lineage_scores"]) > i else None
#                 for i in range(3)
#             },
#             # 🔥额外信息（不影响原逻辑）
#             "confidence": v.get("confidence"),
#             "entropy": v.get("entropy")
#         }
#         for cluster, v in cluster_scores.items()
#     }).T

#     return adata, cluster_annotations, score_summary
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

    


# def compute_cluster_stats(adata, use_raw=True, threshold=0):
#     X = adata.raw.X if (use_raw and adata.raw is not None) else adata.X
#     genes = adata.raw.var_names if (use_raw and adata.raw is not None) else adata.var_names

#     clusters = adata.obs["leiden"].unique()

#     mean_expr = {}
#     pct_expr = {}

#     for c in clusters:
#         idx = adata.obs["leiden"] == c
#         sub = X[idx]

#         mean_expr[c] = np.asarray(sub.mean(axis=0)).ravel()
#         pct_expr[c] = np.asarray((sub > threshold).mean(axis=0)).ravel()

#     mean_df = pd.DataFrame(mean_expr, index=genes).T
#     pct_df = pd.DataFrame(pct_expr, index=genes).T

#     # z-score（按cluster）
#     z_df = (mean_df - mean_df.mean()) / (mean_df.std() + 1e-6)

#     return mean_df, pct_df, z_df