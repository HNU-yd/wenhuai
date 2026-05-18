#!/usr/bin/env python
import os
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import zscore

import sys
from project_paths import add_src_to_path, project_path

add_src_to_path()

from kyn_ahr_v0.gene_sets import (
    get_gene_sets_for_dataset,
    CELLTYPE_MARKERS,
    MYELOID_CELLTYPES,
    COMPOSITE_SCORE_CONFIGS,
    MODULE_SCORE_COLUMNS,
    COMPOSITE_SCORE_COLUMNS,
    ALL_SCORE_COLUMNS,
)


def setup_threads(n: int):
    os.environ["OMP_NUM_THREADS"] = str(n)
    os.environ["OPENBLAS_NUM_THREADS"] = str(n)
    os.environ["MKL_NUM_THREADS"] = str(n)
    os.environ["NUMEXPR_NUM_THREADS"] = str(n)


def present_genes(adata, genes):
    # 大小写不敏感匹配，兼容 human 全大写和 mouse 首字母大写。
    var_map = {str(g).upper(): str(g) for g in adata.var_names}
    out = []

    for g in genes:
        if g in adata.var_names:
            out.append(g)
        elif g.upper() in var_map:
            out.append(var_map[g.upper()])

    return sorted(set(out))


def qc_filter(adata, min_genes=200, max_mito=20.0, min_cells=3):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    before_cells = adata.n_obs
    before_genes = adata.n_vars

    adata = adata[adata.obs["n_genes_by_counts"] >= min_genes].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= max_mito].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)

    print(
        f"[qc] cells {before_cells} -> {adata.n_obs}; "
        f"genes {before_genes} -> {adata.n_vars}"
    )

    return adata


def add_module_scores(adata, outdir: Path, dataset: str):
    rows = []
    gene_sets = get_gene_sets_for_dataset(dataset)

    for score_name, genes in gene_sets.items():
        use = present_genes(adata, genes)
        rows.append({
            "dataset": dataset,
            "score": score_name,
            "n_requested": len(genes),
            "n_present": len(use),
            "genes_present": ",".join(use),
        })

        if len(use) >= 2:
            sc.tl.score_genes(
                adata,
                gene_list=use,
                score_name=score_name,
                use_raw=False,
            )
        else:
            adata.obs[score_name] = 0.0

    # marker scores
    for ct, genes in CELLTYPE_MARKERS.items():
        use = present_genes(adata, genes)
        rows.append({
            "dataset": dataset,
            "score": f"marker_{ct}",
            "n_requested": len(genes),
            "n_present": len(use),
            "genes_present": ",".join(use),
        })

        if len(use) >= 2:
            sc.tl.score_genes(
                adata,
                gene_list=use,
                score_name=f"marker_{ct}",
                use_raw=False,
            )
        else:
            adata.obs[f"marker_{ct}"] = 0.0

    marker_cols = [f"marker_{x}" for x in CELLTYPE_MARKERS.keys()]
    marker_mat = adata.obs[marker_cols].to_numpy(dtype=float)
    labels = list(CELLTYPE_MARKERS.keys())
    best = np.argmax(marker_mat, axis=1)
    adata.obs["cell_type_v0"] = [labels[i] for i in best]

    # composite scores
    for composite_name, parts in COMPOSITE_SCORE_CONFIGS.items():
        z_list = []
        for c in parts:
            arr = adata.obs[c].to_numpy(dtype=float)
            zz = zscore(arr, nan_policy="omit")
            zz = np.nan_to_num(zz, nan=0.0, posinf=0.0, neginf=0.0)
            z_list.append(zz)

        adata.obs[composite_name] = np.vstack(z_list).mean(axis=0)

    report = pd.DataFrame(rows)
    report.to_csv(outdir / f"{dataset}.score_gene_presence.tsv", sep="\t", index=False)


def run_embedding(adata, dataset: str, n_hvg: int):
    print("[HVG]")
    batch_key = "sample_id" if "sample_id" in adata.obs.columns else None

    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_hvg,
            flavor="seurat",
            batch_key=batch_key,
            subset=True,
        )
    except Exception as e:
        print(f"[warn] seurat HVG failed: {repr(e)}")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_hvg,
            flavor="cell_ranger",
            batch_key=batch_key,
            subset=True,
        )

    print("[scale/PCA]")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

    use_rep = "X_pca"

    try:
        import scanpy.external as sce
        if batch_key and adata.obs[batch_key].nunique() > 1:
            print("[Harmony]")
            sce.pp.harmony_integrate(
                adata,
                key=batch_key,
                basis="X_pca",
                adjusted_basis="X_pca_harmony",
            )
            use_rep = "X_pca_harmony"

    except Exception as e:
        print(f"[warn] Harmony skipped: {repr(e)}")

    print("[neighbors/leiden/umap]")
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=40, use_rep=use_rep)
    sc.tl.leiden(adata, resolution=0.6, key_added="leiden_v0")
    sc.tl.umap(adata)


def make_plots(adata, outdir: Path, dataset: str):
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(figdir)

    colors = [
        "group",
        "severity",
        "stage",
        "sample_id",
        "patient_id",
        "cell_type_v0",
        "leiden_v0",
    ] + ALL_SCORE_COLUMNS

    for color in colors:
        if color in adata.obs.columns:
            try:
                sc.pl.umap(
                    adata,
                    color=color,
                    show=False,
                    save=f"_{dataset}_{color}.png",
                )
            except Exception as e:
                print(f"[warn] plot failed {color}: {repr(e)}")


def export_summary_tables(adata, outdir: Path, dataset: str):
    obs = adata.obs.copy()

    base_cols = ["dataset", "sample_id", "patient_id", "group", "severity", "stage"]
    for c in base_cols:
        if c not in obs.columns:
            obs[c] = "unknown"

    score_cols = [c for c in ALL_SCORE_COLUMNS if c in obs.columns]

    print("[export] cell type proportions")
    ct = (
        obs.groupby(base_cols + ["cell_type_v0"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )

    totals = (
        ct.groupby(["dataset", "sample_id"], observed=True)["n_cells"]
        .sum()
        .reset_index(name="sample_total_cells")
    )

    ct = ct.merge(totals, on=["dataset", "sample_id"], how="left")
    ct["fraction"] = ct["n_cells"] / ct["sample_total_cells"]
    ct.to_csv(outdir / f"{dataset}.celltype_proportions.tsv", sep="\t", index=False)

    print("[export] module + composite score pseudobulk")
    pseudo = (
        obs.groupby(base_cols + ["cell_type_v0"], observed=True)[score_cols]
        .mean()
        .reset_index()
    )
    pseudo.to_csv(
        outdir / f"{dataset}.score_pseudobulk_by_sample_celltype.tsv",
        sep="\t",
        index=False,
    )

    myeloid = obs[obs["cell_type_v0"].isin(MYELOID_CELLTYPES)].copy()
    if len(myeloid):
        myeloid_summary = (
            myeloid.groupby(base_cols + ["cell_type_v0"], observed=True)[score_cols]
            .mean()
            .reset_index()
        )
        myeloid_summary.to_csv(
            outdir / f"{dataset}.myeloid_score_pseudobulk.tsv",
            sep="\t",
            index=False,
        )

    print("[export] top Kyn-AHR-high cells")
    q = obs["Kyn_AHR_myeloid_score"].quantile(0.90)
    high = obs.loc[obs["Kyn_AHR_myeloid_score"] >= q].copy()

    keep_cols = base_cols + [
        "cell_type_v0",
        "leiden_v0",
    ] + score_cols

    high[keep_cols].to_csv(
        outdir / f"{dataset}.top10pct_Kyn_AHR_high_cells.tsv",
        sep="\t",
    )

    high_summary = (
        high.groupby(["dataset", "group", "severity", "stage", "cell_type_v0"], observed=True)
        .size()
        .reset_index(name="n_top10pct_cells")
        .sort_values("n_top10pct_cells", ascending=False)
    )
    high_summary.to_csv(
        outdir / f"{dataset}.top10pct_Kyn_AHR_high_summary.tsv",
        sep="\t",
        index=False,
    )


def run_quick_de(adata, outdir: Path, dataset: str):
    print("[DE] quick wilcoxon disease vs control")

    if "severity" not in adata.obs.columns:
        return

    sev = adata.obs["severity"].astype(str)
    adata.obs["v0_binary"] = np.where(
        sev.str.contains("control", case=False),
        "control",
        "disease",
    )

    if adata.obs["v0_binary"].nunique() >= 2:
        try:
            sc.tl.rank_genes_groups(adata, groupby="v0_binary", method="wilcoxon")
            df = sc.get.rank_genes_groups_df(adata, group="disease")
            df.to_csv(
                outdir / f"{dataset}.DE_disease_vs_control.wilcoxon.tsv",
                sep="\t",
                index=False,
            )
        except Exception as e:
            print(f"[warn] global DE failed: {repr(e)}")

    myeloid_mask = adata.obs["cell_type_v0"].isin(MYELOID_CELLTYPES)
    if int(myeloid_mask.sum()) > 50:
        sub = adata[myeloid_mask].copy()
        if sub.obs["v0_binary"].nunique() >= 2:
            try:
                sc.tl.rank_genes_groups(sub, groupby="v0_binary", method="wilcoxon")
                df = sc.get.rank_genes_groups_df(sub, group="disease")
                df.to_csv(
                    outdir / f"{dataset}.DE_myeloid_disease_vs_control.wilcoxon.tsv",
                    sep="\t",
                    index=False,
                )
            except Exception as e:
                print(f"[warn] myeloid DE failed: {repr(e)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--out_root", default=str(project_path("results", "v0")))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min_genes", type=int, default=200)
    ap.add_argument("--max_mito", type=float, default=20.0)
    ap.add_argument("--n_hvg", type=int, default=3000)
    args = ap.parse_args()

    setup_threads(args.threads)

    outdir = Path(args.out_root) / args.dataset
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[read] {args.input}")
    adata = sc.read_h5ad(args.input)

    for c in ["dataset", "sample_id", "patient_id", "group", "severity", "stage", "batch"]:
        if c not in adata.obs.columns:
            adata.obs[c] = "unknown"

    print("[qc]")
    adata = qc_filter(
        adata,
        min_genes=args.min_genes,
        max_mito=args.max_mito,
        min_cells=3,
    )
    adata.write_h5ad(outdir / f"{args.dataset}.01_qc.h5ad", compression="gzip")

    print("[normalize/log1p]")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    print("[scores]")
    add_module_scores(adata, outdir, args.dataset)
    adata.obs.to_csv(outdir / f"{args.dataset}.obs_after_scores.tsv", sep="\t")

    run_embedding(adata, args.dataset, n_hvg=args.n_hvg)
    adata.write_h5ad(outdir / f"{args.dataset}.02_v0_processed.h5ad", compression="gzip")

    make_plots(adata, outdir, args.dataset)
    export_summary_tables(adata, outdir, args.dataset)
    run_quick_de(adata, outdir, args.dataset)

    adata.write_h5ad(outdir / f"{args.dataset}.final_v0.h5ad", compression="gzip")

    print(f"[done] {args.dataset}")


if __name__ == "__main__":
    main()
