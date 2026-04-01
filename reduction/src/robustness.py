def check_10x(path):
    import os
    required = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]
    for f in required:
        if not os.path.exists(os.path.join(path, f)):
            raise FileNotFoundError(f"{f} missing in {path}")


def check_marker_coverage(adata, marker_dict):
    missing = []
    for lineage, info in marker_dict.items():
        genes = info["general"]
        for g in genes:
            if g in adata.var_names:
                if not adata.var.loc[g, "highly_variable"]:
                    missing.append(g)
    return missing