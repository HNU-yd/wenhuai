#!/usr/bin/env python
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

import sys
from project_paths import add_src_to_path, project_path

add_src_to_path()

from kyn_ahr_v0.gene_sets import (
    get_gene_sets_for_dataset,
    COMPOSITE_SCORE_CONFIGS,
)


def present_genes(adata, genes):
    var_map = {str(g).upper(): str(g) for g in adata.var_names}
    out = []
    for g in genes:
        if g in adata.var_names:
            out.append(g)
        elif g.upper() in var_map:
            out.append(var_map[g.upper()])
    return sorted(set(out))


def to_dense_small(x):
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def score_gene_set(adata, genes):
    idx = [adata.var_names.get_loc(g) for g in genes]
    X = to_dense_small(adata.X[:, idx]).astype(np.float32)

    expr = X.mean(axis=1)

    mu = X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    z = ((X - mu) / sd).mean(axis=1)

    return expr, z


def add_composite_scores(adata, composite_name, parts):
    expr_cols = [f"{x}_expr" for x in parts if f"{x}_expr" in adata.obs.columns]
    z_cols = [f"{x}_z" for x in parts if f"{x}_z" in adata.obs.columns]

    if expr_cols:
        mat = adata.obs[expr_cols].to_numpy(dtype=float)
        # 对 module expr 在样本内部标准化后再平均，避免不同模块尺度差异。
        mat_z = (mat - np.nanmean(mat, axis=0)) / np.nanstd(mat, axis=0)
        mat_z = np.nan_to_num(mat_z, nan=0.0, posinf=0.0, neginf=0.0)
        adata.obs[f"{composite_name}_expr"] = mat_z.mean(axis=1)

    if z_cols:
        adata.obs[f"{composite_name}_z"] = adata.obs[z_cols].mean(axis=1)


def process_one(in_h5ad, out_h5ad, dataset="STT0000127", target_sum=1e4, overwrite=False):
    out_h5ad = Path(out_h5ad)
    if out_h5ad.exists() and not overwrite:
        print(f"[skip] exists: {out_h5ad}")
        return

    print(f"[read] {in_h5ad}")
    adata = sc.read_h5ad(in_h5ad)

    # 保留原始 counts 的 QC 信息，不在 scored 文件里复制 counts layer，节省空间。
    print("[normalize/log1p]")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    gene_sets = get_gene_sets_for_dataset(dataset)
    presence_rows = []

    print("[score modules]")
    for score_name, genes in gene_sets.items():
        use = present_genes(adata, genes)

        presence_rows.append({
            "sample_id": adata.obs["sample_id"].iloc[0],
            "score": score_name,
            "n_requested": len(genes),
            "n_present": len(use),
            "genes_present": ",".join(use),
        })

        if len(use) < 2:
            adata.obs[f"{score_name}_expr"] = 0.0
            adata.obs[f"{score_name}_z"] = 0.0
            continue

        expr, z = score_gene_set(adata, use)
        adata.obs[f"{score_name}_expr"] = expr
        adata.obs[f"{score_name}_z"] = z

    print("[score composites]")
    for comp, parts in COMPOSITE_SCORE_CONFIGS.items():
        add_composite_scores(adata, comp, parts)

    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    print(f"[write] {out_h5ad}")
    adata.write_h5ad(out_h5ad, compression="gzip")

    presence = pd.DataFrame(presence_rows)
    presence_path = out_h5ad.with_suffix(".gene_presence.tsv")
    presence.to_csv(presence_path, sep="\t", index=False)

    print(f"[done] {out_h5ad}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_dir", default=str(project_path("data", "spatial_h5ad", "STT0000127", "bin50_tissuecut")))
    ap.add_argument("--out_dir", default=str(project_path("data", "spatial_h5ad", "STT0000127", "bin50_tissuecut_scored")))
    ap.add_argument("--dataset", default="STT0000127")
    ap.add_argument("--target_sum", type=float, default=1e4)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--samples", nargs="*", default=None)
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(in_dir.glob("*.bin50.tissuecut.h5ad"))
    if args.samples:
        wanted = set(args.samples)
        files = [f for f in files if f.name.split(".")[0] in wanted]

    reports = []
    for f in files:
        sample_id = f.name.split(".")[0]
        out = out_dir / f"{sample_id}.bin50.tissuecut.scored.h5ad"

        process_one(
            in_h5ad=f,
            out_h5ad=out,
            dataset=args.dataset,
            target_sum=args.target_sum,
            overwrite=args.overwrite,
        )

        reports.append({
            "sample_id": sample_id,
            "input_h5ad": str(f),
            "output_h5ad": str(out),
        })

    pd.DataFrame(reports).to_csv(
        out_dir / "STT0000127_spatial_scoring_report.tsv",
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
