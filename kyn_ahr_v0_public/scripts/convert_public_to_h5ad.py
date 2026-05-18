#!/usr/bin/env python3
import argparse
import gzip
import re
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


DEFAULT_ROOT = Path(__file__).resolve().parents[1]


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(DEFAULT_ROOT))
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--datasets", nargs="*", default=["GSE167029", "GSE166489", "GSE183716"])
    return ap.parse_args()


def safe_extract_tar(tar_path: Path, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path) as tar:
        for member in tar.getmembers():
            target = outdir / member.name
            if not str(target.resolve()).startswith(str(outdir.resolve())):
                raise RuntimeError(f"Unsafe tar member: {member.name}")
        tar.extractall(outdir)


def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t")
    meta["gsm"] = meta["gsm"].astype(str)
    meta["sample_id"] = meta["sample_id"].astype(str)
    return meta


def infer_meta_for_file(path: Path, meta: pd.DataFrame, dataset: str):
    name = path.name
    subset = meta[(meta.dataset == dataset) & (meta.include_v0 == 1)].copy()
    # Prefer GSM accession if present in filename.
    m = re.search(r"GSM\d+", name)
    if m:
        hit = subset[subset.gsm == m.group(0)]
        if len(hit):
            return hit.iloc[0].to_dict()
    # Then sample labels.
    for _, row in subset.iterrows():
        sid = str(row.sample_id)
        candidates = [sid, sid.replace(".", "-"), sid.replace(".", "_"), sid.split("_")[-1]]
        for c in candidates:
            if c and c in name:
                return row.to_dict()
    # GSE183716 h5 files use Sample1..Sample4.
    m = re.search(r"Sample(\d+)", name)
    if m and dataset == "GSE183716":
        hit = subset[subset.sample_id.str.startswith(f"Sample{m.group(1)}_")]
        if len(hit):
            return hit.iloc[0].to_dict()
    return None


def read_10x_h5(path: Path):
    try:
        x = sc.read_10x_h5(path, gex_only=True)
    except Exception:
        x = sc.read_10x_h5(path)
    x.var_names_make_unique()
    return x


def read_csv_counts(path: Path):
    # Assumes genes/features in rows and cells in columns. This is slower but robust for GEO CSV matrices.
    print(f"[read_csv_counts] {path}")
    df = pd.read_csv(path, index_col=0)
    # Drop non-expression feature columns if present.
    non_expr_cols = {"gene", "genes", "feature", "features", "gene_id", "gene_name"}
    keep_cols = [c for c in df.columns if str(c).lower() not in non_expr_cols]
    df = df[keep_cols]
    x = sparse.csr_matrix(df.T.values)
    obs = pd.DataFrame(index=df.columns.astype(str))
    var = pd.DataFrame(index=df.index.astype(str))
    a = ad.AnnData(X=x, obs=obs, var=var)
    a.var_names_make_unique()
    return a


def read_any_matrix(path: Path):
    if path.suffix == ".h5":
        return read_10x_h5(path)
    if path.suffix == ".csv" or path.name.endswith(".csv.gz"):
        return read_csv_counts(path)
    return None


def collect_matrix_files(dataset_dir: Path):
    files = []
    for pat in ["**/*.h5", "**/*.csv", "**/*.csv.gz"]:
        files.extend(dataset_dir.glob(pat))
    # Avoid duplicate files in macOS hidden dirs or metadata csvs.
    files = [f for f in sorted(set(files)) if not f.name.startswith("._")]
    return files


def main():
    args = parse_args()
    root = Path(args.root)
    raw = root / "raw"
    processed = root / "processed"
    processed.mkdir(parents=True, exist_ok=True)
    meta_path = Path(args.metadata) if args.metadata else root / "metadata" / "metadata_public_v0.tsv"
    meta = load_metadata(meta_path)
    out = Path(args.out) if args.out else processed / "public_mis_c_combined_raw.h5ad"

    adatas = []
    used = []
    for dataset in args.datasets:
        ddir = raw / dataset
        # Extract tar files if not already extracted.
        for tar_path in ddir.glob("*.tar"):
            extracted = ddir / "extracted"
            if not extracted.exists() or not any(extracted.iterdir()):
                print(f"[extract] {tar_path} -> {extracted}")
                safe_extract_tar(tar_path, extracted)
        files = collect_matrix_files(ddir)
        print(f"[{dataset}] candidate matrix files: {len(files)}")
        for f in files:
            info = infer_meta_for_file(f, meta, dataset)
            if info is None:
                print(f"[skip: no metadata match] {f}")
                continue
            try:
                a = read_any_matrix(f)
            except Exception as e:
                print(f"[skip: read failed] {f}: {e}")
                continue
            if a is None:
                continue
            for k, v in info.items():
                a.obs[k] = v
            a.obs["source_file"] = str(f)
            a.obs_names = [f"{info['dataset']}::{info['sample_id']}::{bc}" for bc in a.obs_names]
            adatas.append(a)
            used.append((dataset, info["gsm"], info["sample_id"], str(f), a.n_obs, a.n_vars))
            print(f"[loaded] {dataset} {info['sample_id']} cells={a.n_obs} genes={a.n_vars}")
    if not adatas:
        raise SystemExit("No matrices loaded. Check extraction and metadata filename matching.")
    combined = ad.concat(adatas, join="outer", label="concat_batch", fill_value=0, index_unique=None)
    # Ensure sparse matrix.
    if not sparse.issparse(combined.X):
        combined.X = sparse.csr_matrix(combined.X)
    combined.var_names_make_unique()
    combined.write_h5ad(out, compression="gzip")
    pd.DataFrame(used, columns=["dataset", "gsm", "sample_id", "source_file", "n_cells", "n_genes"]).to_csv(
        processed / "loaded_samples.tsv", sep="\t", index=False
    )
    print(f"[write] {out}")
    print(f"[write] {processed / 'loaded_samples.tsv'}")


if __name__ == "__main__":
    main()
