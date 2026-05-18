#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import zscore

try:
    import yaml
except Exception:
    yaml = None

GENESETS = {
    "KYN_metabolism_score": ["IDO1", "IDO2", "TDO2", "AFMID", "KMO", "KYNU", "HAAO", "QPRT", "AADAT"],
    "AHR_response_score": ["AHR", "ARNT", "CYP1A1", "CYP1B1", "AHRR", "TIPARP", "NQO1", "ALDH1A1", "PTGS2"],
    "myeloid_inflammation_score": ["S100A8", "S100A9", "IL1B", "TNF", "CXCL8", "CCL2", "CCL3", "CCL4", "NFKBIA", "STAT1", "IRF1"],
    "chemotaxis_score": ["CCR1", "CCR2", "CCR5", "CXCR4", "CXCL8", "CCL2", "CCL3", "CCL4", "CCL5", "ITGAM", "ITGB2"],
    "antigen_presentation_score": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CD74", "CIITA"],
}

MARKERS = {
    "T_cell": ["CD3D", "CD3E", "TRAC"],
    "CD4_T": ["CD3D", "CD4", "IL7R", "CCR7"],
    "CD8_T": ["CD8A", "CD8B", "GZMK", "GZMB"],
    "NK": ["NKG7", "GNLY", "KLRD1", "NCAM1"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "CD74"],
    "plasmablast": ["MZB1", "JCHAIN", "XBP1", "PRDM1"],
    "CD14_mono": ["LYZ", "LST1", "S100A8", "S100A9", "FCN1", "CD14"],
    "FCGR3A_mono": ["LYZ", "LST1", "FCGR3A", "MS4A7", "CTSS"],
    "DC": ["FCER1A", "CST3", "CLEC10A", "CD1C"],
    "pDC": ["LILRA4", "TCF4", "IRF7", "GZMB"],
    "neutrophil_like": ["S100A8", "S100A9", "CSF3R", "FCGR3B", "CXCR2"],
    "platelet": ["PPBP", "PF4", "GP9"],
}

MYELOID_LABELS = ["CD14_mono", "FCGR3A_mono", "DC", "pDC", "neutrophil_like"]


DEFAULT_ROOT = Path(__file__).resolve().parents[1]


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_h5ad", default=str(DEFAULT_ROOT / "processed" / "public_mis_c_combined_raw.h5ad"))
    ap.add_argument("--out_h5ad", default=str(DEFAULT_ROOT / "processed" / "public_mis_c_v0_scanpy.h5ad"))
    ap.add_argument("--figdir", default=str(DEFAULT_ROOT / "figures"))
    ap.add_argument("--tables", default=str(DEFAULT_ROOT / "tables"))
    ap.add_argument("--min_genes", type=int, default=200)
    ap.add_argument("--max_genes", type=int, default=7000)
    ap.add_argument("--max_pct_mt", type=float, default=20)
    ap.add_argument("--min_counts", type=int, default=500)
    ap.add_argument("--n_hvg", type=int, default=3000)
    ap.add_argument("--n_pcs", type=int, default=50)
    ap.add_argument("--leiden_resolution", type=float, default=0.8)
    ap.add_argument("--batch_key", default="dataset")
    ap.add_argument("--group_key", default="group_v0")
    return ap.parse_args()


def present_genes(adata, genes):
    return [g for g in genes if g in adata.var_names]


def score_set(adata, name, genes):
    genes = present_genes(adata, genes)
    if len(genes) == 0:
        adata.obs[name] = np.nan
        print(f"[score missing] {name}: no genes present")
    else:
        sc.tl.score_genes(adata, gene_list=genes, score_name=name, use_raw=False)
        print(f"[score] {name}: {len(genes)} genes")


def z_col(s):
    x = pd.to_numeric(s, errors="coerce").astype(float)
    sd = np.nanstd(x)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(x)), index=s.index)
    return (x - np.nanmean(x)) / sd


def main():
    args = parse_args()
    figdir = Path(args.figdir); figdir.mkdir(parents=True, exist_ok=True)
    tables = Path(args.tables); tables.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(figdir)
    sc.settings.verbosity = 2

    adata = sc.read_h5ad(args.in_h5ad)
    adata.var_names_make_unique()
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    adata.layers["counts"] = adata.X.copy()

    # QC
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    qc = adata.obs[["dataset", "sample_id", "group_v0", "n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
    qc.to_csv(tables / "cell_qc_before_filter.tsv", sep="\t")
    keep = (
        (adata.obs["n_genes_by_counts"] >= args.min_genes) &
        (adata.obs["n_genes_by_counts"] <= args.max_genes) &
        (adata.obs["total_counts"] >= args.min_counts) &
        (adata.obs["pct_counts_mt"] <= args.max_pct_mt)
    )
    print(f"[QC] cells before={adata.n_obs}, after={int(keep.sum())}")
    adata = adata[keep].copy()

    # Normalize/log/HVG/PCA
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, batch_key=args.batch_key, flavor="seurat")
    except Exception as e:
        print(f"[HVG fallback] {e}")
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, flavor="cell_ranger")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=args.n_pcs, svd_solver="arpack")

    # Harmony if installed; otherwise use PCA.
    rep = "X_pca"
    try:
        import scanpy.external as sce
        sce.pp.harmony_integrate(adata, key=args.batch_key, basis="X_pca")
        rep = "X_pca_harmony"
        print("[batch] using Harmony representation: X_pca_harmony")
    except Exception as e:
        print(f"[batch] Harmony unavailable or failed; using X_pca. Reason: {e}")

    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=args.n_pcs, use_rep=rep)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.leiden_resolution, key_added="leiden")

    # Basic marker-based cell-type labels.
    for label, genes in MARKERS.items():
        score_set(adata, f"marker_{label}", genes)
    marker_cols = [f"marker_{k}" for k in MARKERS]
    marker_mat = adata.obs[marker_cols].fillna(-999)
    inv = {f"marker_{k}": k for k in MARKERS}
    adata.obs["cell_type_v0"] = marker_mat.idxmax(axis=1).map(inv).astype("category")

    # Kyn-AHR/pathway scores on all cells.
    for name, genes in GENESETS.items():
        score_set(adata, name, genes)
    # Placeholder proxy until pySCENIC/DoRothEA regulon is run.
    adata.obs["AHR_regulon_score_proxy"] = adata.obs["AHR_response_score"]
    components = ["KYN_metabolism_score", "AHR_response_score", "AHR_regulon_score_proxy", "myeloid_inflammation_score"]
    z = pd.concat([z_col(adata.obs[c]) for c in components], axis=1)
    z.columns = components
    adata.obs["Kyn_AHR_myeloid_score"] = z.mean(axis=1)

    # Tables.
    adata.obs[["dataset", "sample_id", "patient_id", "group_v0", "severity", "stage", "cell_type_v0"] + list(GENESETS.keys()) + ["AHR_regulon_score_proxy", "Kyn_AHR_myeloid_score"]].to_csv(
        tables / "cell_metadata_with_v0_scores.tsv", sep="\t"
    )
    prop = (
        adata.obs.groupby(["dataset", "sample_id", "group_v0", "cell_type_v0"], observed=True)
        .size().rename("n_cells").reset_index()
    )
    prop["sample_total"] = prop.groupby(["dataset", "sample_id"], observed=True)["n_cells"].transform("sum")
    prop["fraction"] = prop["n_cells"] / prop["sample_total"]
    prop.to_csv(tables / "celltype_fraction_by_sample.tsv", sep="\t", index=False)

    score_summary = (
        adata.obs.groupby(["dataset", "group_v0", "cell_type_v0"], observed=True)[list(GENESETS.keys()) + ["AHR_regulon_score_proxy", "Kyn_AHR_myeloid_score"]]
        .agg(["mean", "median", "std", "count"])
    )
    score_summary.to_csv(tables / "score_summary_by_group_celltype.tsv", sep="\t")

    # Marker genes by Leiden cluster.
    try:
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
        sc.get.rank_genes_groups_df(adata, group=None).to_csv(tables / "rank_genes_leiden.tsv", sep="\t", index=False)
    except Exception as e:
        print(f"[rank_genes] skipped: {e}")

    # Figures.
    sc.pl.umap(adata, color=["dataset", args.group_key, "cell_type_v0", "leiden"], save="_v0_overview.png", show=False)
    sc.pl.umap(adata, color=["KYN_metabolism_score", "AHR_response_score", "myeloid_inflammation_score", "Kyn_AHR_myeloid_score"], save="_v0_scores.png", show=False)
    sc.pl.violin(adata, keys=["KYN_metabolism_score", "AHR_response_score", "myeloid_inflammation_score", "Kyn_AHR_myeloid_score"], groupby=args.group_key, rotation=90, save="_v0_scores_by_group.png", show=False)
    sc.pl.dotplot(adata, var_names={k: present_genes(adata.raw.to_adata(), v) for k, v in GENESETS.items()}, groupby="cell_type_v0", save="_v0_genesets_by_celltype.png", show=False)

    # Myeloid-focused subset.
    myeloid = adata[adata.obs["cell_type_v0"].isin(MYELOID_LABELS)].copy()
    if myeloid.n_obs > 100:
        sc.pp.neighbors(myeloid, n_neighbors=20, n_pcs=min(args.n_pcs, 30), use_rep=rep if rep in myeloid.obsm else "X_pca")
        sc.tl.umap(myeloid)
        sc.tl.leiden(myeloid, resolution=0.6, key_added="myeloid_leiden")
        sc.pl.umap(myeloid, color=[args.group_key, "cell_type_v0", "myeloid_leiden", "Kyn_AHR_myeloid_score"], save="_v0_myeloid.png", show=False)
        myeloid.write_h5ad(Path(args.out_h5ad).with_name("public_mis_c_v0_myeloid.h5ad"), compression="gzip")
        myeloid.obs.to_csv(tables / "myeloid_cells_with_scores.tsv", sep="\t")

    adata.write_h5ad(args.out_h5ad, compression="gzip")
    print(f"[write] {args.out_h5ad}")
    print(f"[write tables] {tables}")
    print(f"[write figures] {figdir}")


if __name__ == "__main__":
    main()
