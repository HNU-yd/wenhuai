#!/usr/bin/env python
import argparse
import gzip
import gc
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

from project_paths import project_path


def parse_gem_stream(path, bin_size=50, use_exon=False, max_records=0):
    """
    读取 GEM 文件并按 bin_size 聚合。
    支持 header:
      geneID x y MIDCount ExonCount
    """
    counts = defaultdict(float)
    bin_xy = {}
    n_records = 0
    n_bad = 0
    header_seen = False

    with gzip.open(path, "rt", errors="ignore") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) < 4:
                n_bad += 1
                continue

            # header
            if not header_seen and parts[0].lower() in {"geneid", "gene"}:
                header_seen = True
                continue

            try:
                gene = parts[0]
                x = int(float(parts[1]))
                y = int(float(parts[2]))
                mid = float(parts[3])
                exon = float(parts[4]) if use_exon and len(parts) > 4 else mid
                value = exon if use_exon else mid
            except Exception:
                n_bad += 1
                continue

            xb = x // bin_size
            yb = y // bin_size
            spot = f"{xb}_{yb}"

            counts[(spot, gene)] += value

            if spot not in bin_xy:
                bin_xy[spot] = {
                    "x_bin": xb,
                    "y_bin": yb,
                    "x_center": xb * bin_size + bin_size / 2.0,
                    "y_center": yb * bin_size + bin_size / 2.0,
                }

            n_records += 1

            if max_records > 0 and n_records >= max_records:
                break

    return counts, bin_xy, n_records, n_bad


def build_anndata_from_counts(counts, bin_xy, sample_meta, min_counts_per_bin=5, min_genes_per_bin=3):
    spots = sorted({k[0] for k in counts.keys()})
    genes = sorted({k[1] for k in counts.keys()})

    spot_to_i = {s: i for i, s in enumerate(spots)}
    gene_to_j = {g: j for j, g in enumerate(genes)}

    row = []
    col = []
    data = []

    for (spot, gene), val in counts.items():
        if val <= 0:
            continue
        row.append(spot_to_i[spot])
        col.append(gene_to_j[gene])
        data.append(val)

    X = sparse.csr_matrix(
        (data, (row, col)),
        shape=(len(spots), len(genes)),
        dtype=np.float32,
    )

    obs = pd.DataFrame(index=[f"{sample_meta['sample_id']}::{s}" for s in spots])
    obs["spot"] = spots

    obs["sample_id"] = sample_meta["sample_id"]
    obs["group"] = sample_meta["group"]
    obs["stage"] = sample_meta["stage"]
    obs["stage_order"] = sample_meta["stage_order"]
    obs["file_role"] = sample_meta["preferred_role"]

    obs["x_bin"] = [bin_xy[s]["x_bin"] for s in spots]
    obs["y_bin"] = [bin_xy[s]["y_bin"] for s in spots]
    obs["x_center"] = [bin_xy[s]["x_center"] for s in spots]
    obs["y_center"] = [bin_xy[s]["y_center"] for s in spots]

    var = pd.DataFrame(index=genes)
    var["gene_symbols"] = genes

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # QC
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()

    keep = (adata.obs["n_counts"] >= min_counts_per_bin) & (adata.obs["n_genes"] >= min_genes_per_bin)
    adata = adata[keep].copy()

    adata.uns["spatial_bin_info"] = {
        "bin_size": sample_meta["bin_size"],
        "source_path": sample_meta["preferred_path"],
        "use_exon": sample_meta["use_exon"],
    }

    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample_design", default=str(project_path("results", "spatial_inventory", "STT0000127_spatial_sample_design.tsv")))
    ap.add_argument("--out_dir", default=str(project_path("data", "spatial_h5ad", "STT0000127", "bin50_tissuecut")))
    ap.add_argument("--bin_size", type=int, default=50)
    ap.add_argument("--use_exon", action="store_true")
    ap.add_argument("--max_records", type=int, default=0, help="0 means full file")
    ap.add_argument("--min_counts_per_bin", type=float, default=5)
    ap.add_argument("--min_genes_per_bin", type=int, default=3)
    ap.add_argument("--samples", nargs="*", default=None)
    args = ap.parse_args()

    sample_design = pd.read_csv(args.sample_design, sep="\t")

    if args.samples:
        sample_design = sample_design[sample_design["sample_id"].isin(args.samples)].copy()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    report = []

    for _, r in sample_design.iterrows():
        sample_id = r["sample_id"]
        path = r["preferred_path"]

        if not isinstance(path, str) or not path:
            print(f"[skip] empty path for {sample_id}")
            continue

        path = Path(path)
        if not path.exists():
            print(f"[skip] missing {path}")
            continue

        out_h5ad = out_dir / f"{sample_id}.bin{args.bin_size}.tissuecut.h5ad"

        print("=" * 80)
        print(f"[sample] {sample_id}")
        print(f"[path] {path}")
        print(f"[out] {out_h5ad}")

        counts, bin_xy, n_records, n_bad = parse_gem_stream(
            path,
            bin_size=args.bin_size,
            use_exon=args.use_exon,
            max_records=args.max_records,
        )

        sample_meta = {
            "sample_id": sample_id,
            "group": r["group"],
            "stage": r["stage"],
            "stage_order": r["stage_order"],
            "preferred_role": r["preferred_role"],
            "preferred_path": str(path),
            "bin_size": args.bin_size,
            "use_exon": bool(args.use_exon),
        }

        adata = build_anndata_from_counts(
            counts,
            bin_xy,
            sample_meta=sample_meta,
            min_counts_per_bin=args.min_counts_per_bin,
            min_genes_per_bin=args.min_genes_per_bin,
        )

        adata.write_h5ad(out_h5ad, compression="gzip")

        report.append({
            "sample_id": sample_id,
            "group": r["group"],
            "stage": r["stage"],
            "path": str(path),
            "out_h5ad": str(out_h5ad),
            "n_records": n_records,
            "n_bad": n_bad,
            "n_bins": adata.n_obs,
            "n_genes": adata.n_vars,
            "total_counts": float(adata.obs["n_counts"].sum()),
            "mean_counts_per_bin": float(adata.obs["n_counts"].mean()) if adata.n_obs else 0,
            "median_counts_per_bin": float(adata.obs["n_counts"].median()) if adata.n_obs else 0,
        })

        del counts
        del bin_xy
        del adata
        gc.collect()

    report_df = pd.DataFrame(report)
    report_path = out_dir / "STT0000127_bin_h5ad_build_report.tsv"
    report_df.to_csv(report_path, sep="\t", index=False)

    print(f"[done] report: {report_path}")
    print(report_df.to_string(index=False))


if __name__ == "__main__":
    main()
