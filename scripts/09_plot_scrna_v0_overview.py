#!/usr/bin/env python
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from project_paths import project_path, project_root


MYELOID_TYPES = ["CD16_Monocyte", "DC", "Monocyte", "pDC"]

SCORE_COLS = [
    "KYN_metabolism_score",
    "AHR_response_score",
    "AHR_regulon_proxy_score",
    "myeloid_inflammation_score",
    "chemotaxis_score",
    "antigen_presentation_score",
    "Kyn_AHR_axis_score",
    "Kyn_AHR_myeloid_score",
    "myeloid_effector_score",
]

COMPARISONS = [
    {
        "dataset": "CNP0005824",
        "label": "CVB3 vs Control",
        "positive": "CVB3_myocarditis",
        "negative": "control",
    },
    {
        "dataset": "CNP0005824",
        "label": "CVB3 vs IVIG",
        "positive": "CVB3_myocarditis",
        "negative": "IVIG_treated",
    },
    {
        "dataset": "GSE166489",
        "label": "MIS-C vs pediatric healthy",
        "positive": "MIS-C",
        "negative": "pediatric_healthy",
    },
    {
        "dataset": "GSE167029",
        "label": "MIS-C_MYO vs control",
        "positive": "MIS-C_MYO",
        "negative": "control",
    },
]


def load_pseudobulk(root: Path, datasets):
    rows = []
    for ds in datasets:
        f = root / "results/v0" / ds / f"{ds}.myeloid_score_pseudobulk.tsv"
        if not f.exists():
            print(f"[skip] missing {f}")
            continue
        df = pd.read_csv(f, sep="\t")
        df["dataset"] = ds
        rows.append(df)

    if not rows:
        raise SystemExit("[error] no pseudobulk files found")

    out = pd.concat(rows, ignore_index=True)

    for s in SCORE_COLS:
        if s not in out.columns:
            out[s] = np.nan

    return out


def compute_effect_table(df):
    rows = []

    for comp in COMPARISONS:
        ds = comp["dataset"]
        pos = comp["positive"]
        neg = comp["negative"]

        sub = df[df["dataset"] == ds].copy()

        for ct in MYELOID_TYPES:
            one = sub[sub["cell_type_v0"] == ct].copy()

            pos_df = one[one["group"] == pos]
            neg_df = one[one["group"] == neg]

            if pos_df.empty or neg_df.empty:
                continue

            for score in SCORE_COLS:
                x = pos_df[score].dropna().to_numpy(dtype=float)
                y = neg_df[score].dropna().to_numpy(dtype=float)

                if len(x) == 0 or len(y) == 0:
                    continue

                rows.append({
                    "dataset": ds,
                    "comparison": comp["label"],
                    "positive_group": pos,
                    "negative_group": neg,
                    "cell_type_v0": ct,
                    "score": score,
                    "positive_mean": float(np.mean(x)),
                    "negative_mean": float(np.mean(y)),
                    "diff_positive_minus_negative": float(np.mean(x) - np.mean(y)),
                    "n_positive": len(x),
                    "n_negative": len(y),
                })

    return pd.DataFrame(rows)


def short_score_name(score):
    mapping = {
        "KYN_metabolism_score": "KYN",
        "AHR_response_score": "AHR resp.",
        "AHR_regulon_proxy_score": "AHR proxy",
        "myeloid_inflammation_score": "Myeloid infl.",
        "chemotaxis_score": "Chemotaxis",
        "antigen_presentation_score": "Antigen pres.",
        "Kyn_AHR_axis_score": "KYN/AHR axis",
        "Kyn_AHR_myeloid_score": "KYN/AHR/myeloid",
        "myeloid_effector_score": "Myeloid eff.",
    }
    return mapping.get(score, score)


def dataset_row_order(dataset):
    comp_labels = [x["label"] for x in COMPARISONS if x["dataset"] == dataset]
    rows = []
    for comp in comp_labels:
        for ct in MYELOID_TYPES:
            rows.append(comp + " | " + ct)
    return rows


def plot_boundary_heatmap_per_dataset(effect, dataset, out_png):
    """
    每个数据集单独画一张 heatmap，
    并且使用该数据集自己的色条范围。
    """
    sub = effect[effect["dataset"] == dataset].copy()
    if sub.empty:
        print(f"[skip] no effect rows for {dataset}")
        return

    sub["row"] = sub["comparison"] + " | " + sub["cell_type_v0"]
    sub["score_short"] = sub["score"].map(short_score_name)

    row_order = dataset_row_order(dataset)
    col_order = [short_score_name(x) for x in SCORE_COLS]

    mat = (
        sub.pivot_table(
            index="row",
            columns="score_short",
            values="diff_positive_minus_negative",
            aggfunc="mean",
        )
        .reindex(index=row_order, columns=col_order)
    )

    values = mat.to_numpy(dtype=float)
    finite = values[np.isfinite(values)]

    if len(finite) == 0:
        print(f"[skip] no finite values for {dataset}")
        return

    # 每个数据集单独定色条范围
    vmax = np.nanquantile(np.abs(finite), 0.95)
    if vmax == 0 or not np.isfinite(vmax):
        vmax = np.nanmax(np.abs(finite))
    if vmax == 0 or not np.isfinite(vmax):
        vmax = 1.0

    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig_h = max(4.5, 0.55 * len(mat.index))
    fig_w = max(9, 0.95 * len(mat.columns))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(values, aspect="auto", cmap="RdBu_r", norm=norm)

    ax.set_xticks(np.arange(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right")

    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels(mat.index)

    ax.set_title(
        f"{dataset} boundary effect heatmap\n(target group - reference group, dataset-specific scale)",
        fontsize=14
    )

    # 网格
    ax.set_xticks(np.arange(-0.5, len(mat.columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(mat.index), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    # 数字标注
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            v = values[i, j]
            if np.isfinite(v):
                ax.text(
                    j,
                    i,
                    f"{v:.2f}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="black" if abs(v) < 0.6 * vmax else "white",
                )

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Mean score difference", rotation=270, labelpad=16)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=260)
    plt.close(fig)


def plot_cnp_group_celltype(df, out_png, score="Kyn_AHR_myeloid_score"):
    sub = df[(df["dataset"] == "CNP0005824") & (df["cell_type_v0"].isin(MYELOID_TYPES))].copy()

    group_order = ["control", "IVIG_treated", "CVB3_myocarditis"]
    sub["group"] = pd.Categorical(sub["group"], categories=group_order, ordered=True)
    sub = sub.sort_values(["cell_type_v0", "group"])

    fig, axes = plt.subplots(1, len(MYELOID_TYPES), figsize=(4 * len(MYELOID_TYPES), 4), sharey=True)

    if len(MYELOID_TYPES) == 1:
        axes = [axes]

    for ax, ct in zip(axes, MYELOID_TYPES):
        one = sub[sub["cell_type_v0"] == ct].copy()

        data = [one[one["group"] == g][score].dropna().to_numpy(dtype=float) for g in group_order]
        means = [np.mean(x) if len(x) else np.nan for x in data]

        ax.boxplot(data, labels=group_order, showfliers=False)
        for i, vals in enumerate(data, start=1):
            if len(vals):
                jitter = np.random.default_rng(0).normal(0, 0.04, size=len(vals))
                ax.scatter(np.full(len(vals), i) + jitter, vals, s=18, alpha=0.7)

        ax.plot(np.arange(1, len(group_order) + 1), means, marker="o", linewidth=2)
        ax.set_title(ct)
        ax.tick_params(axis="x", rotation=35)
        ax.axhline(0, linestyle="--", linewidth=0.8)

    axes[0].set_ylabel(score)
    fig.suptitle(f"CNP0005824 myeloid pseudobulk score: {score}", fontsize=14)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=260)
    plt.close(fig)


def plot_human_boundary_box(df, out_png, score="Kyn_AHR_myeloid_score"):
    datasets = ["GSE166489", "GSE167029"]
    fig, axes = plt.subplots(len(datasets), len(MYELOID_TYPES), figsize=(4 * len(MYELOID_TYPES), 7), sharey=False)

    for r, ds in enumerate(datasets):
        sub = df[(df["dataset"] == ds) & (df["cell_type_v0"].isin(MYELOID_TYPES))].copy()

        if ds == "GSE166489":
            group_order = ["pediatric_healthy", "MIS-C", "adult_healthy"]
        else:
            group_order = ["control", "acute_infection", "MIS-C", "MIS-C_MYO", "KD"]

        sub["group"] = pd.Categorical(sub["group"], categories=group_order, ordered=True)
        sub = sub.sort_values(["cell_type_v0", "group"])

        for c, ct in enumerate(MYELOID_TYPES):
            ax = axes[r, c]
            one = sub[sub["cell_type_v0"] == ct].copy()

            data = [one[one["group"] == g][score].dropna().to_numpy(dtype=float) for g in group_order]
            means = [np.mean(x) if len(x) else np.nan for x in data]

            ax.boxplot(data, labels=group_order, showfliers=False)
            for i, vals in enumerate(data, start=1):
                if len(vals):
                    jitter = np.random.default_rng(1).normal(0, 0.04, size=len(vals))
                    ax.scatter(np.full(len(vals), i) + jitter, vals, s=14, alpha=0.7)

            ax.plot(np.arange(1, len(group_order) + 1), means, marker="o", linewidth=1.5)
            ax.axhline(0, linestyle="--", linewidth=0.8)
            ax.set_title(f"{ds}\n{ct}")
            ax.tick_params(axis="x", rotation=45)

            if c == 0:
                ax.set_ylabel(score)

    fig.suptitle(f"Human PBMC boundary analysis: {score}", fontsize=14)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=260)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(project_root()))
    ap.add_argument("--datasets", nargs="*", default=["CNP0005824", "GSE166489", "GSE167029"])
    ap.add_argument("--out_dir", default=str(project_path("results", "v0_figures")))
    args = ap.parse_args()

    root = Path(args.root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_pseudobulk(root, args.datasets)
    df.to_csv(out_dir / "scrna_myeloid_pseudobulk_merged.tsv", sep="\t", index=False)

    effect = compute_effect_table(df)
    effect.to_csv(out_dir / "scrna_boundary_effect_table.tsv", sep="\t", index=False)

    # 改成每个数据集单独一张 heatmap
    for ds in ["CNP0005824", "GSE166489", "GSE167029"]:
        plot_boundary_heatmap_per_dataset(
            effect,
            ds,
            out_dir / f"Fig_scRNA_1_{ds}_boundary_heatmap.png",
        )

    plot_cnp_group_celltype(
        df,
        out_dir / "Fig_scRNA_2_CNP0005824_Kyn_AHR_myeloid_score.png",
        score="Kyn_AHR_myeloid_score",
    )

    plot_cnp_group_celltype(
        df,
        out_dir / "Fig_scRNA_2b_CNP0005824_myeloid_inflammation_score.png",
        score="myeloid_inflammation_score",
    )

    plot_human_boundary_box(
        df,
        out_dir / "Fig_scRNA_3_human_boundary_Kyn_AHR_myeloid_score.png",
        score="Kyn_AHR_myeloid_score",
    )

    plot_human_boundary_box(
        df,
        out_dir / "Fig_scRNA_3b_human_boundary_myeloid_inflammation_score.png",
        score="myeloid_inflammation_score",
    )

    print(f"[done] figures -> {out_dir}")
    print("[files]")
    for f in sorted(out_dir.glob("*.png")):
        print(f)


if __name__ == "__main__":
    main()
