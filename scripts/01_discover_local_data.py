#!/usr/bin/env python
import argparse
from pathlib import Path

import pandas as pd

from project_paths import project_path


EXCLUDE_DIR_NAMES = {
    "h5ad_per_sample",
    "h5ad_merged",
    "metadata",
    "results",
    "logs",
    "__pycache__",
}

DEFAULT_DATASETS = [
    "CNP0005824",
    "GSE166489",
    "GSE167029",
    "GSE180045",
    "GSE183716",
    # STT0000127 是 GEM 空间数据，默认不进第二层 PBMC V0。
]


def should_skip(path: Path) -> bool:
    return bool(set(path.parts) & EXCLUDE_DIR_NAMES)


def classify_file(path: Path) -> str:
    name = path.name.lower()

    if name.endswith(".h5ad"):
        return "h5ad"
    if name.endswith(".h5") or name.endswith(".hdf5"):
        return "h5_or_hdf5"
    if name.endswith(".rds"):
        return "seurat_rds"
    if name in {"matrix.mtx", "matrix.mtx.gz"}:
        return "10x_mtx_matrix"
    if name in {"barcodes.tsv", "barcodes.tsv.gz"}:
        return "10x_mtx_barcodes"
    if name in {"features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"}:
        return "10x_mtx_features"
    if name.endswith(".tar") or name.endswith(".tar.gz") or name.endswith(".tgz"):
        return "archive_tar"
    if name.endswith(".zip"):
        return "archive_zip"
    if name.endswith(".soft") or name.endswith(".soft.gz"):
        return "geo_soft"
    if "series_matrix" in name and (name.endswith(".txt") or name.endswith(".txt.gz")):
        return "geo_series_matrix"
    return "other"


def is_10x_mtx_dir(path: Path) -> bool:
    if not path.is_dir():
        return False

    names = {x.name.lower() for x in path.iterdir() if x.is_file()}
    has_matrix = "matrix.mtx" in names or "matrix.mtx.gz" in names
    has_barcodes = "barcodes.tsv" in names or "barcodes.tsv.gz" in names
    has_features = (
        "features.tsv" in names or "features.tsv.gz" in names
        or "genes.tsv" in names or "genes.tsv.gz" in names
    )
    return has_matrix and has_barcodes and has_features


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_root", default=str(project_path("data")))
    ap.add_argument("--out", default=str(project_path("data", "metadata", "local_data_inventory.tsv")))
    ap.add_argument("--datasets", nargs="*", default=DEFAULT_DATASETS)
    args = ap.parse_args()

    data_root = Path(args.data_root)
    rows = []

    for dataset in args.datasets:
        ds_root = data_root / dataset
        if not ds_root.exists():
            print(f"[warn] missing dataset dir: {ds_root}")
            continue

        for f in sorted(ds_root.rglob("*")):
            if should_skip(f):
                continue
            if not f.is_file():
                continue

            ftype = classify_file(f)
            if ftype == "other":
                continue

            rows.append({
                "dataset": dataset,
                "record_type": "file",
                "path": str(f),
                "parent": str(f.parent),
                "name": f.name,
                "file_type": ftype,
                "size_bytes": f.stat().st_size,
                "size_gb": round(f.stat().st_size / 1024**3, 5),
            })

        for d in sorted(ds_root.rglob("*")):
            if should_skip(d):
                continue
            if d.is_dir() and is_10x_mtx_dir(d):
                rows.append({
                    "dataset": dataset,
                    "record_type": "directory",
                    "path": str(d),
                    "parent": str(d.parent),
                    "name": d.name,
                    "file_type": "10x_mtx_dir",
                    "size_bytes": 0,
                    "size_gb": 0.0,
                })

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"[done] inventory written: {args.out}")

    if df.empty:
        print("[warn] no candidate files found")
        return

    print("\n[file type summary]")
    print(
        df.groupby(["dataset", "file_type"])
        .size()
        .reset_index(name="n")
        .sort_values(["dataset", "file_type"])
        .to_string(index=False)
    )

    print("\n[high priority candidates]")
    hp = df[df["file_type"].isin(["h5ad", "h5_or_hdf5", "seurat_rds", "10x_mtx_dir"])]
    if len(hp):
        print(hp[["dataset", "file_type", "size_gb", "path"]].to_string(index=False))
    else:
        print("[warn] no h5ad / h5 / rds / 10x_mtx_dir found")


if __name__ == "__main__":
    main()
