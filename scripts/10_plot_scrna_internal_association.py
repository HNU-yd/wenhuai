#!/usr/bin/env python
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu

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

SHORT_NAMES = {
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

# 这些是非自证型关键关系：
# 尽量避免用 composite 去和组成它的同一模块直接作唯一证据。
KEY_PAIRS = [
    ("KYN_metabolism_score", "AHR_response_score"),
    ("KYN_metabolism_score", "AHR_regulon_proxy_score"),
    ("AHR_response_score", "myeloid_inflammation_score"),
    ("AHR_regulon_proxy_score", "myeloid_inflammation_score"),
    ("Kyn_AHR_axis_score", "myeloid_inflammation_score"),
    ("Kyn_AHR_axis_score", "chemotaxis_score"),
    ("Kyn_AHR_axis_score", "antigen_presentation_score"),
]


def short(x):
    return SHORT_NAMES.get(x, x)


def load_pseudobulk(root: Path, dataset: str):
    f = root / "results/v0" / dataset / f"{dataset}.myeloid_score_pseudobulk.tsv"
    if not f.exists():
        raise FileNotFoundError(f"Missing pseudobulk file: {f}")
    df = pd.read_csv(f, sep="\t")
    df["dataset"] = dataset

    for c in SCORE_COLS:
        if c not in df.columns:
            df[c] = np.nan

    return df


def safe_spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]

    if len(x) < 4:
        return np.nan, np.nan, len(x)

    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan, np.nan, len(x)

    r, p = spearmanr(x, y)
    return float(r), float(p), int(len(x))


def compute_corr_long(df):
    rows = []

    for (dataset, ct), sub in df.groupby(["dataset", "cell_type_v0"], observed=True):
        for xscore in SCORE_COLS:
            for yscore in SCORE_COLS:
                r, p, n = safe_spearman(sub[xscore], sub[yscore])
                rows.append({
                    "dataset": dataset,
                    "cell_type_v0": ct,
                    "x_score": xscore,
                    "y_score": yscore,
                    "spearman_r": r,
                    "p_value": p,
                    "n": n,
                })

    return pd.DataFrame(rows)


def plot_corr_heatmap(df, dataset, cell_type, out_png):
    sub = df[(df["dataset"] == dataset) & (df["cell_type_v0"] == cell_type)].copy()
    if sub.empty:
        print(f"[skip] no data for {dataset} {cell_type}")
        return

    mat = sub[SCORE_COLS].corr(method="spearman")
    mat = mat.reindex(index=SCORE_COLS, columns=SCORE_COLS)

    values = mat.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    im = ax.imshow(values, cmap="RdBu_r", vmin=-1, vmax=1)

    ax.set_xticks(np.arange(len(SCORE_COLS)))
    ax.set_xticklabels([short(x) for x in SCORE_COLS], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(SCORE_COLS)))
    ax.set_yticklabels([short(x) for x in SCORE_COLS])

    ax.set_title(f"{dataset} | {cell_type}\nModule internal Spearman correlation", fontsize=13)

    ax.set_xticks(np.arange(-0.5, len(SCORE_COLS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(SCORE_COLS), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            v = values[i, j]
            if np.isfinite(v):
                ax.text(
                    j, i, f"{v:.2f}",
                    ha="center", va="center",
                    fontsize=7,
                    color="white" if abs(v) > 0.6 else "black",
                )

    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Spearman r", rotation=270, labelpad=14)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=260)
    plt.close(fig)


def plot_key_scatter(df, dataset, cell_type, xscore, yscore, out_png):
    sub = df[(df["dataset"] == dataset) & (df["cell_type_v0"] == cell_type)].copy()
    if sub.empty:
        return

    x = sub[xscore].to_numpy(dtype=float)
    y = sub[yscore].to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    sub = sub.iloc[np.where(ok)[0]].copy()

    if len(sub) < 4:
        return

    r, p, n = safe_spearman(x, y)

    groups = sorted(sub["group"].dropna().unique())

    fig, ax = plt.subplots(figsize=(5, 4.5))

    for g in groups:
        one = sub[sub["group"] == g]
        ax.scatter(
            one[xscore],
            one[yscore],
            s=45,
            alpha=0.8,
            label=g,
        )

    # 简单趋势线
    if len(x) >= 3 and np.nanstd(x) > 0:
        coef = np.polyfit(x, y, 1)
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
        ys = coef[0] * xs + coef[1]
        ax.plot(xs, ys, linewidth=1.5)

    ax.axhline(0, linestyle="--", linewidth=0.8)
    ax.axvline(0, linestyle="--", linewidth=0.8)

    ax.set_xlabel(short(xscore))
    ax.set_ylabel(short(yscore))
    ax.set_title(
        f"{dataset} | {cell_type}\n{short(xscore)} vs {short(yscore)}\nSpearman r={r:.2f}, p={p:.3g}, n={n}",
        fontsize=11,
    )
    ax.legend(fontsize=8, frameon=False)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=260)
    plt.close(fig)


def high_low_by_score(df, dataset, cell_type, split_score, target_scores):
    sub = df[(df["dataset"] == dataset) & (df["cell_type_v0"] == cell_type)].copy()
    if sub.empty:
        return pd.DataFrame()

    v = sub[split_score].to_numpy(dtype=float)
    if len(v) < 6 or np.nanstd(v) == 0:
        return pd.DataFrame()

    q50 = np.nanmedian(v)
    sub["split_group"] = np.where(sub[split_score] >= q50, "High", "Low")

    rows = []
    for target in target_scores:
        hi = sub[sub["split_group"] == "High"][target].dropna().to_numpy(dtype=float)
        lo = sub[sub["split_group"] == "Low"][target].dropna().to_numpy(dtype=float)

        if len(hi) == 0 or len(lo) == 0:
            continue

        try:
            stat, p = mannwhitneyu(hi, lo, alternative="two-sided")
        except Exception:
            p = np.nan

        rows.append({
            "dataset": dataset,
            "cell_type_v0": cell_type,
            "split_score": split_score,
            "target_score": target,
            "split_threshold_median": q50,
            "high_mean": float(np.mean(hi)),
            "low_mean": float(np.mean(lo)),
            "diff_high_minus_low": float(np.mean(hi) - np.mean(lo)),
            "n_high": len(hi),
            "n_low": len(lo),
            "p_value": p,
        })

    return pd.DataFrame(rows)


def plot_high_low_boxplot(df, dataset, cell_type, split_score, target_scores, out_png):
    sub = df[(df["dataset"] == dataset) & (df["cell_type_v0"] == cell_type)].copy()
    if sub.empty:
        return

    v = sub[split_score].to_numpy(dtype=float)
    if len(v) < 6 or np.nanstd(v) == 0:
        return

    q50 = np.nanmedian(v)
    sub["split_group"] = np.where(sub[split_score] >= q50, "High", "Low")

    n = len(target_scores)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, target in zip(axes, target_scores):
        data = [
            sub[sub["split_group"] == "Low"][target].dropna().to_numpy(dtype=float),
            sub[sub["split_group"] == "High"][target].dropna().to_numpy(dtype=float),
        ]

        ax.boxplot(data, labels=["Low", "High"], showfliers=False)
        for i, vals in enumerate(data, start=1):
            if len(vals):
                jitter = np.random.default_rng(2).normal(0, 0.04, size=len(vals))
                ax.scatter(np.full(len(vals), i) + jitter, vals, s=35, alpha=0.8)

        ax.axhline(0, linestyle="--", linewidth=0.8)
        ax.set_title(short(target))
        ax.set_xlabel(f"{short(split_score)} group")

    axes[0].set_ylabel("module score")
    fig.suptitle(f"{dataset} | {cell_type}\nHigh vs Low by {short(split_score)}", fontsize=12)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=260)
    plt.close(fig)


def summarize_key_pairs(corr_long):
    rows = []

    for _, r in corr_long.iterrows():
        x = r["x_score"]
        y = r["y_score"]
        if (x, y) in KEY_PAIRS:
            rows.append(r.to_dict())

    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(project_root()))
    ap.add_argument("--datasets", nargs="*", default=["GSE166489_child"])
    ap.add_argument("--out_dir", default=str(project_path("results", "v0_assoc")))
    args = ap.parse_args()

    root = Path(args.root)
    out_dir = Path(args.out_dir)
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    all_df = []
    for ds in args.datasets:
        try:
            df = load_pseudobulk(root, ds)
            all_df.append(df)
        except FileNotFoundError as e:
            print(f"[skip] {e}")

    if not all_df:
        raise SystemExit("[error] no datasets loaded")

    df = pd.concat(all_df, ignore_index=True)
    df.to_csv(out_dir / "internal_association_input_pseudobulk.tsv", sep="\t", index=False)

    corr_long = compute_corr_long(df)
    corr_long.to_csv(out_dir / "internal_association_spearman_all_pairs.tsv", sep="\t", index=False)

    key_corr = summarize_key_pairs(corr_long)
    key_corr.to_csv(out_dir / "internal_association_key_pairs.tsv", sep="\t", index=False)

    # 1. correlation heatmaps
    for ds in args.datasets:
        for ct in MYELOID_TYPES:
            plot_corr_heatmap(
                df,
                dataset=ds,
                cell_type=ct,
                out_png=fig_dir / f"Fig_assoc_corr_{ds}_{ct}.png",
            )

    # 2. key scatter plots
    for ds in args.datasets:
        for ct in MYELOID_TYPES:
            for xscore, yscore in KEY_PAIRS:
                out_png = fig_dir / f"Fig_assoc_scatter_{ds}_{ct}_{xscore}_vs_{yscore}.png"
                plot_key_scatter(df, ds, ct, xscore, yscore, out_png)

    # 3. high-low boxplots
    highlow_rows = []
    target_scores = [
        "myeloid_inflammation_score",
        "chemotaxis_score",
        "antigen_presentation_score",
        "AHR_response_score",
        "AHR_regulon_proxy_score",
    ]

    for ds in args.datasets:
        for ct in MYELOID_TYPES:
            hl = high_low_by_score(
                df,
                dataset=ds,
                cell_type=ct,
                split_score="Kyn_AHR_axis_score",
                target_scores=target_scores,
            )
            if len(hl):
                highlow_rows.append(hl)

            plot_high_low_boxplot(
                df,
                dataset=ds,
                cell_type=ct,
                split_score="Kyn_AHR_axis_score",
                target_scores=target_scores,
                out_png=fig_dir / f"Fig_assoc_highlow_{ds}_{ct}_by_Kyn_AHR_axis.png",
            )

    if highlow_rows:
        highlow = pd.concat(highlow_rows, ignore_index=True)
    else:
        highlow = pd.DataFrame()

    highlow.to_csv(out_dir / "internal_association_highlow_by_Kyn_AHR_axis.tsv", sep="\t", index=False)

    print(f"[done] tables -> {out_dir}")
    print(f"[done] figures -> {fig_dir}")

    print("\n[key correlations]")
    if len(key_corr):
        show = key_corr[
            (key_corr["dataset"].isin(args.datasets))
            & (key_corr["x_score"].isin([x for x, y in KEY_PAIRS]))
        ].copy()
        print(show.to_string(index=False))

    print("\n[figure files]")
    for f in sorted(fig_dir.glob("*.png"))[:30]:
        print(f)


if __name__ == "__main__":
    main()
