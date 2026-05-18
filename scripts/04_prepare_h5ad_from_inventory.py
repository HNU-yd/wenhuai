#!/usr/bin/env python
import os
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import scanpy as sc
from tqdm import tqdm

import sys
from project_paths import add_src_to_path, project_path

add_src_to_path()

from kyn_ahr_v0.metadata_rules import (
    infer_meta_from_path,
    load_geo_sample_table,
    safe_name,
)


SUPPORTED_TYPES = {"h5ad", "h5_or_hdf5", "10x_mtx_dir"}


def set_worker_threads():
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"


def read_candidate(path: str, file_type: str):
    if file_type == "h5ad":
        return sc.read_h5ad(path)

    if file_type == "h5_or_hdf5":
        try:
            return sc.read_10x_h5(path, gex_only=True)
        except TypeError:
            return sc.read_10x_h5(path)
        except Exception as e:
            raise RuntimeError(f"Failed to read as 10x_h5: {path}; error={repr(e)}")

    if file_type == "10x_mtx_dir":
        try:
            return sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
        except Exception:
            return sc.read_10x_mtx(path, var_names="gene_ids", cache=False)

    raise ValueError(f"unsupported file_type: {file_type}")


def standardize_adata(adata, meta):
    adata.var_names_make_unique()

    sample_id = safe_name(meta["sample_id"])
    adata.obs_names = [f"{sample_id}::{x}" for x in adata.obs_names]

    for k, v in meta.items():
        if k != "source_file":
            adata.obs[k] = v

    adata.uns["source_file"] = meta["source_file"]
    return adata


def convert_one(job):
    set_worker_threads()

    dataset, path, file_type, out_root, geo_table_path = job
    geo_df = load_geo_sample_table(geo_table_path)
    meta = infer_meta_from_path(dataset, path, geo_df=geo_df)

    sample_id = safe_name(meta["sample_id"])
    out_dir = Path(out_root) / dataset
    out_dir.mkdir(parents=True, exist_ok=True)

    out_h5ad = out_dir / f"{sample_id}.h5ad"
    out_meta = out_dir / f"{sample_id}.metadata.tsv"

    if out_h5ad.exists() and out_h5ad.stat().st_size > 1024:
        return {
            "status": "skip",
            "dataset": dataset,
            "file_type": file_type,
            "path": path,
            "out": str(out_h5ad),
            "n_obs": "",
            "n_vars": "",
        }

    try:
        adata = read_candidate(path, file_type)
        adata = standardize_adata(adata, meta)

        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        adata.write_h5ad(out_h5ad, compression="gzip")
        pd.DataFrame([meta]).to_csv(out_meta, sep="\t", index=False)

        return {
            "status": "ok",
            "dataset": dataset,
            "file_type": file_type,
            "path": path,
            "out": str(out_h5ad),
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
        }

    except Exception as e:
        return {
            "status": "error",
            "dataset": dataset,
            "file_type": file_type,
            "path": path,
            "out": "",
            "n_obs": "",
            "n_vars": "",
            "error": repr(e),
        }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inventory", default=str(project_path("data", "metadata", "local_data_inventory.tsv")))
    ap.add_argument("--geo_table", default=str(project_path("data", "metadata", "geo_sample_table.tsv")))
    ap.add_argument("--out_root", default=str(project_path("data", "h5ad_per_sample")))
    ap.add_argument("--report", default=str(project_path("data", "metadata", "prepare_h5ad_from_inventory_report.tsv")))
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--datasets", nargs="*", default=None)
    args = ap.parse_args()

    inv = pd.read_csv(args.inventory, sep="\t")

    cand = inv[inv["file_type"].isin(SUPPORTED_TYPES)].copy()

    if args.datasets:
        cand = cand[cand["dataset"].isin(args.datasets)].copy()

    # 避免错误地把 tenx_standard 根目录本身当样本。
    cand = cand[~cand["path"].str.endswith("/tenx_standard")].copy()

    cand = cand.drop_duplicates(subset=["dataset", "path", "file_type"])
    cand = cand.sort_values(["dataset", "file_type", "path"])

    print(f"[candidates] {len(cand)}")
    if len(cand):
        print(cand[["dataset", "file_type", "size_gb", "path"]].to_string(index=False))

    jobs = [
        (r["dataset"], r["path"], r["file_type"], args.out_root, args.geo_table)
        for _, r in cand.iterrows()
    ]

    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(convert_one, j) for j in jobs]
        for fut in tqdm(as_completed(futs), total=len(futs)):
            res = fut.result()
            print(res)
            results.append(res)

    df = pd.DataFrame(results)
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.report, sep="\t", index=False)

    ok = df[df["status"].isin(["ok", "skip"])] if len(df) else df
    geo_df = load_geo_sample_table(args.geo_table)

    metas = []
    for _, r in ok.iterrows():
        metas.append(infer_meta_from_path(r["dataset"], r["path"], geo_df=geo_df))

    if metas:
        pd.DataFrame(metas).drop_duplicates().to_csv(
            project_path("data", "metadata", "sample_sheet.auto.tsv"),
            sep="\t",
            index=False,
        )

    print(f"[done] report: {args.report}")


if __name__ == "__main__":
    main()
