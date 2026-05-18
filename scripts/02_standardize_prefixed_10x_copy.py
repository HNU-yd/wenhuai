#!/usr/bin/env python
from pathlib import Path
import argparse
import shutil

from project_paths import project_path


EXCLUDE_PARTS = {
    "h5ad_per_sample",
    "h5ad_merged",
    "tenx_standard",
    "metadata",
    "__pycache__",
}


def should_skip(path: Path) -> bool:
    return bool(set(path.parts) & EXCLUDE_PARTS)


def copy_one(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists() and dst.stat().st_size == src.stat().st_size:
        print(f"[skip] exists: {dst}")
        return

    tmp = dst.with_name(dst.name + ".tmp")
    if tmp.exists():
        tmp.unlink()

    print(f"[copy] {src} -> {dst}")
    shutil.copy2(src, tmp)
    tmp.rename(dst)


def standardize_dataset(data_root: Path, dataset: str):
    ds_root = data_root / dataset
    if not ds_root.exists():
        print(f"[skip] missing {ds_root}")
        return 0

    out_root = ds_root / "tenx_standard"
    matrix_files = []

    for mtx in ds_root.rglob("*_matrix.mtx.gz"):
        if should_skip(mtx):
            continue
        matrix_files.append(mtx)

    print(f"==== {dataset}: {len(matrix_files)} prefixed matrix files ====")

    n_ok = 0
    for mtx in sorted(matrix_files):
        prefix = mtx.name.replace("_matrix.mtx.gz", "")

        barcodes = mtx.parent / f"{prefix}_barcodes.tsv.gz"
        features = mtx.parent / f"{prefix}_features.tsv.gz"
        genes = mtx.parent / f"{prefix}_genes.tsv.gz"

        if not barcodes.exists():
            print(f"[skip] missing barcodes: {mtx}")
            continue

        if features.exists():
            feature_src = features
        elif genes.exists():
            feature_src = genes
        else:
            print(f"[skip] missing features/genes: {mtx}")
            continue

        sample = prefix
        out = out_root / sample

        copy_one(mtx, out / "matrix.mtx.gz")
        copy_one(barcodes, out / "barcodes.tsv.gz")
        copy_one(feature_src, out / "features.tsv.gz")

        print(f"[ok] {dataset}/{sample}")
        n_ok += 1

    return n_ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_root", default=str(project_path("data")))
    ap.add_argument("--datasets", nargs="*", default=[
        "CNP0005824",
        "GSE166489",
        "GSE167029",
        "GSE180045",
    ])
    args = ap.parse_args()

    data_root = Path(args.data_root)
    total = 0

    for dataset in args.datasets:
        total += standardize_dataset(data_root, dataset)

    print(f"[done] standardized samples: {total}")


if __name__ == "__main__":
    main()
