import os

def check_10x(path):
    required = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]
    for f in required:
        if not os.path.exists(os.path.join(path, f)):
            raise FileNotFoundError(f"{f} missing in {path}")