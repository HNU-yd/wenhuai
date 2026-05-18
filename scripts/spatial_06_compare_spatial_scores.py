#!/usr/bin/env python
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import wilcoxon

from project_paths import project_path


SCORES = [
    "KYN_metabolism_score_expr",
    "AHR_response_score_expr",
    "AHR_regulon_proxy_score_expr",
    "myeloid_inflammation_score_expr",
    "Kyn_AHR_axis_score_expr",
    "Kyn_AHR_myeloid_score_expr",
    "Kyn_AHR_myeloid_score_z",
]


def load_obs(f):
    adata = sc.read_h5ad(f)
    obs = adata.obs.copy()
    return obs


def summarize_one(f):
    obs = load_obs(f)

    sample_id = obs["sample_id"].iloc[0]
    group = obs["group"].iloc[0]
    stage = obs["stage"].iloc[0]
    stage_order = int(obs["stage_order"].iloc[0])

    rows = []

    for score in SCORES:
        if score not in obs.columns:
            continue

        v = obs[score].astype(float).to_numpy()
        v = v[np.isfinite(v)]

        if len(v) == 0:
            continue

        q95 = float(np.quantile(v, 0.95))
        q99 = float(np.quantile(v, 0.99))
        top5 = v[v >= q95]
        top1 = v[v >= q99]

        rows.append({
            "sample_id": sample_id,
            "group": group,
            "stage": stage,
            "stage_order": stage_order,
            "score": score,
            "n_bins": len(v),
            "mean": float(np.mean(v)),
            "median": float(np.median(v)),
            "q75": float(np.quantile(v, 0.75)),
            "q90": float(np.quantile(v, 0.90)),
            "q95": q95,
            "q99": q99,
            "top5_mean": float(np.mean(top5)) if len(top5) else np.nan,
            "top1_mean": float(np.mean(top1)) if len(top1) else np.nan,
            "max": float(np.max(v)),
        })

    return rows


def build_control_thresholds(summary):
    control = summary[summary["group"] == "control"].copy()
    rows = []

    for score, sub in control.groupby("score"):
        # 当前只有一个 Control 样本，直接用 Control 的 q95/q99 作为阈值。
        r = sub.iloc[0]
        rows.append({
            "score": score,
            "control_q90_threshold": r["q90"],
            "control_q95_threshold": r["q95"],
            "control_q99_threshold": r["q99"],
        })

    return pd.DataFrame(rows)


def hotspot_fraction(files, thresholds):
    threshold_map = {
        r["score"]: {
            "q90": r["control_q90_threshold"],
            "q95": r["control_q95_threshold"],
            "q99": r["control_q99_threshold"],
        }
        for _, r in thresholds.iterrows()
    }

    rows = []

    for f in files:
        obs = load_obs(f)

        sample_id = obs["sample_id"].iloc[0]
        group = obs["group"].iloc[0]
        stage = obs["stage"].iloc[0]
        stage_order = int(obs["stage_order"].iloc[0])

        for score in SCORES:
            if score not in obs.columns or score not in threshold_map:
                continue

            v = obs[score].astype(float).to_numpy()
            v = v[np.isfinite(v)]
            if len(v) == 0:
                continue

            for level in ["q90", "q95", "q99"]:
                thr = threshold_map[score][level]
                rows.append({
                    "sample_id": sample_id,
                    "group": group,
                    "stage": stage,
                    "stage_order": stage_order,
                    "score": score,
                    "threshold_level": f"control_{level}",
                    "threshold": thr,
                    "n_bins": len(v),
                    "n_hotspot_bins": int((v > thr).sum()),
                    "hotspot_fraction": float((v > thr).mean()),
                    "hotspot_mean_score": float(np.mean(v[v > thr])) if (v > thr).sum() else np.nan,
                })

    return pd.DataFrame(rows)


def day_matched_delta(summary, value_col):
    rows = []

    wide = summary.pivot_table(
        index=["stage_order", "score"],
        columns="group",
        values=value_col,
        aggfunc="mean",
    ).reset_index()

    for _, r in wide.iterrows():
        day = int(r["stage_order"])
        score = r["score"]

        if day == 0:
            continue

        cvb3 = r.get("CVB3_myocarditis", np.nan)
        ivig = r.get("IVIG_treated", np.nan)

        if pd.isna(cvb3) or pd.isna(ivig):
            continue

        rows.append({
            "metric": value_col,
            "stage_order": day,
            "score": score,
            "CVB3": cvb3,
            "IVIG": ivig,
            "delta_CVB3_minus_IVIG": cvb3 - ivig,
        })

    return pd.DataFrame(rows)


def paired_wilcoxon(delta):
    rows = []

    for (metric, score), sub in delta.groupby(["metric", "score"]):
        vals = sub["delta_CVB3_minus_IVIG"].dropna().to_numpy()

        # 全 0 或无变异时，wilcoxon 没有意义。
        if len(vals) >= 2 and np.nanstd(vals) > 0 and not np.allclose(vals, 0):
            try:
                stat, p = wilcoxon(vals)
            except Exception:
                stat, p = np.nan, np.nan
        else:
            stat, p = np.nan, np.nan

        rows.append({
            "metric": metric,
            "score": score,
            "n_days": len(vals),
            "mean_delta_CVB3_minus_IVIG": float(np.mean(vals)) if len(vals) else np.nan,
            "median_delta_CVB3_minus_IVIG": float(np.median(vals)) if len(vals) else np.nan,
            "n_positive_days": int((vals > 0).sum()) if len(vals) else 0,
            "wilcoxon_p_value": p,
        })

    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_dir", default=str(project_path("data", "spatial_h5ad", "STT0000127", "bin50_tissuecut_scored")))
    ap.add_argument("--out_dir", default=str(project_path("results", "spatial_summary", "STT0000127")))
    args = ap.parse_args()

    files = sorted(Path(args.in_dir).glob("*.scored.h5ad"))

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for f in files:
        print(f"[summarize] {f}")
        all_rows.extend(summarize_one(f))

    summary = pd.DataFrame(all_rows)
    summary.to_csv(
        out_dir / "STT0000127_spatial_score_sample_summary.tsv",
        sep="\t",
        index=False,
    )

    thresholds = build_control_thresholds(summary)
    thresholds.to_csv(
        out_dir / "STT0000127_control_hotspot_thresholds.tsv",
        sep="\t",
        index=False,
    )

    hotspot = hotspot_fraction(files, thresholds)
    hotspot.to_csv(
        out_dir / "STT0000127_spatial_hotspot_fraction.tsv",
        sep="\t",
        index=False,
    )

    delta_tables = []

    for metric in ["mean", "median", "q90", "q95", "q99", "top5_mean", "top1_mean"]:
        delta_tables.append(day_matched_delta(summary, metric))

    # hotspot fraction delta
    for level in ["control_q90", "control_q95", "control_q99"]:
        sub = hotspot[hotspot["threshold_level"] == level].copy()
        d = day_matched_delta(sub.rename(columns={"hotspot_fraction": f"hotspot_fraction_{level}"}), f"hotspot_fraction_{level}")
        delta_tables.append(d)

    delta = pd.concat(delta_tables, ignore_index=True)
    delta.to_csv(
        out_dir / "STT0000127_CVB3_vs_IVIG_day_matched_delta.tsv",
        sep="\t",
        index=False,
    )

    paired = paired_wilcoxon(delta)
    paired.to_csv(
        out_dir / "STT0000127_CVB3_vs_IVIG_paired_wilcoxon.tsv",
        sep="\t",
        index=False,
    )

    print("[done]")
    print("\n[paired wilcoxon: selected]")
    selected = paired[
        paired["score"].isin([
            "Kyn_AHR_myeloid_score_expr",
            "Kyn_AHR_axis_score_expr",
            "myeloid_inflammation_score_expr",
            "AHR_response_score_expr",
            "KYN_metabolism_score_expr",
        ])
    ].copy()
    print(selected.to_string(index=False))


if __name__ == "__main__":
    main()
