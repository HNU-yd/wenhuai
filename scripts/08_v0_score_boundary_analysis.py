#!/usr/bin/env python
import argparse
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests

from project_paths import project_root


SCORE_COLS = [
    "KYN_metabolism_score",
    "AHR_response_score",
    "AHR_regulon_proxy_score",
    "myeloid_inflammation_score",
    "chemotaxis_score",
    "antigen_presentation_score",
    "Kyn_AHR_myeloid_score",
    "Kyn_AHR_axis_score",
    "myeloid_effector_score",
]

MYELOID_TYPES = ["Monocyte", "CD16_Monocyte", "DC", "pDC"]


def load_all_myeloid(root: Path, datasets):
    rows = []
    for ds in datasets:
        f = root / "results/v0" / ds / f"{ds}.myeloid_score_pseudobulk.tsv"
        if not f.exists():
            print(f"[skip] missing: {f}")
            continue
        df = pd.read_csv(f, sep="\t")
        df["dataset"] = ds
        rows.append(df)

    if not rows:
        raise SystemExit("[error] no myeloid_score_pseudobulk files found")

    out = pd.concat(rows, ignore_index=True)

    for c in SCORE_COLS:
        if c not in out.columns:
            out[c] = np.nan

    return out


def group_means(df):
    keep_scores = [c for c in SCORE_COLS if c in df.columns]
    return (
        df.groupby(["dataset", "group", "severity", "stage", "cell_type_v0"], observed=True)[keep_scores]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
    )


def pairwise_tests(df):
    rows = []

    for (ds, ct), sub in df.groupby(["dataset", "cell_type_v0"], observed=True):
        groups = sorted([x for x in sub["group"].dropna().unique() if x != "unknown"])

        for g1, g2 in combinations(groups, 2):
            a = sub[sub["group"] == g1]
            b = sub[sub["group"] == g2]

            for score in SCORE_COLS:
                if score not in sub.columns:
                    continue

                x = a[score].dropna().to_numpy()
                y = b[score].dropna().to_numpy()

                if len(x) < 2 or len(y) < 2:
                    continue

                try:
                    stat, p = mannwhitneyu(x, y, alternative="two-sided")
                except Exception:
                    stat, p = np.nan, np.nan

                rows.append({
                    "dataset": ds,
                    "cell_type_v0": ct,
                    "group_1": g1,
                    "group_2": g2,
                    "score": score,
                    "mean_1": float(np.mean(x)),
                    "mean_2": float(np.mean(y)),
                    "diff_1_minus_2": float(np.mean(x) - np.mean(y)),
                    "n_1": int(len(x)),
                    "n_2": int(len(y)),
                    "test": "mannwhitneyu",
                    "p_value": p,
                })

    out = pd.DataFrame(rows)
    if len(out):
        out["fdr"] = multipletests(out["p_value"].fillna(1.0), method="fdr_bh")[1]
    return out


def gse166489_paired_acute_recovery(df):
    sub = df[df["dataset"] == "GSE166489"].copy()
    if sub.empty:
        return pd.DataFrame()

    # patient_id 已经由 metadata 修成 P3/P4；如果老结果没修，也从 sample_id 兜底抽取。
    sub["paired_patient_id"] = sub["patient_id"].astype(str)
    extracted = sub["sample_id"].astype(str).str.extract(r"(P\d+)\.\d+")[0]
    sub.loc[extracted.notna(), "paired_patient_id"] = extracted[extracted.notna()]

    paired = sub[sub["paired_patient_id"].isin(["P3", "P4"])].copy()
    if paired.empty:
        return pd.DataFrame()

    rows = []
    for (patient, ct), one in paired.groupby(["paired_patient_id", "cell_type_v0"], observed=True):
        stages = set(one["stage"].astype(str))
        if not {"acute", "recovery"}.issubset(stages):
            continue

        acute = one[one["stage"] == "acute"].iloc[0]
        recovery = one[one["stage"] == "recovery"].iloc[0]

        for score in SCORE_COLS:
            if score not in one.columns:
                continue
            rows.append({
                "dataset": "GSE166489",
                "paired_patient_id": patient,
                "cell_type_v0": ct,
                "score": score,
                "acute_value": acute[score],
                "recovery_value": recovery[score],
                "delta_acute_minus_recovery": acute[score] - recovery[score],
            })

    detail = pd.DataFrame(rows)
    if detail.empty:
        return detail

    summary_rows = []
    for (ct, score), one in detail.groupby(["cell_type_v0", "score"], observed=True):
        vals = one["delta_acute_minus_recovery"].dropna().to_numpy()
        if len(vals) == 0:
            continue

        if len(vals) >= 2:
            try:
                stat, p = wilcoxon(vals)
            except Exception:
                stat, p = np.nan, np.nan
        else:
            stat, p = np.nan, np.nan

        summary_rows.append({
            "dataset": "GSE166489",
            "cell_type_v0": ct,
            "score": score,
            "n_pairs": len(vals),
            "mean_delta_acute_minus_recovery": float(np.mean(vals)),
            "median_delta_acute_minus_recovery": float(np.median(vals)),
            "wilcoxon_p_value": p,
        })

    summary = pd.DataFrame(summary_rows)
    if len(summary):
        summary["wilcoxon_fdr"] = multipletests(
            summary["wilcoxon_p_value"].fillna(1.0),
            method="fdr_bh",
        )[1]

    detail["table_type"] = "detail"
    summary["table_type"] = "summary"

    return detail, summary


def build_boundary_flags(means, pairwise):
    rows = []

    # 这里只给研究叙事需要的方向性判断。
    def get_mean(ds, group, ct, score):
        m = means.copy()
        # means 是 multi-index columns flatten 前的表，这里要求输入前先 flatten。
        hit = m[
            (m["dataset"] == ds)
            & (m["group"] == group)
            & (m["cell_type_v0"] == ct)
        ]
        col = f"{score}__mean"
        if hit.empty or col not in hit.columns:
            return np.nan
        return float(hit.iloc[0][col])

    checks = [
        {
            "dataset": "CNP0005824",
            "claim": "CVB3_myocarditis > control",
            "positive_group": "CVB3_myocarditis",
            "negative_group": "control",
            "expected_direction": "positive",
        },
        {
            "dataset": "CNP0005824",
            "claim": "CVB3_myocarditis > IVIG_treated",
            "positive_group": "CVB3_myocarditis",
            "negative_group": "IVIG_treated",
            "expected_direction": "positive",
        },
        {
            "dataset": "GSE166489",
            "claim": "MIS-C > pediatric_healthy",
            "positive_group": "MIS-C",
            "negative_group": "pediatric_healthy",
            "expected_direction": "positive",
        },
        {
            "dataset": "GSE167029",
            "claim": "MIS-C_MYO > control",
            "positive_group": "MIS-C_MYO",
            "negative_group": "control",
            "expected_direction": "positive",
        },
    ]

    for check in checks:
        for ct in MYELOID_TYPES:
            for score in [
                "Kyn_AHR_myeloid_score",
                "Kyn_AHR_axis_score",
                "myeloid_inflammation_score",
                "AHR_response_score",
                "KYN_metabolism_score",
            ]:
                a = get_mean(check["dataset"], check["positive_group"], ct, score)
                b = get_mean(check["dataset"], check["negative_group"], ct, score)
                diff = a - b if np.isfinite(a) and np.isfinite(b) else np.nan

                if not np.isfinite(diff):
                    support = "not_evaluable"
                elif diff > 0:
                    support = "supports_expected_direction"
                elif diff < 0:
                    support = "opposes_expected_direction"
                else:
                    support = "neutral"

                rows.append({
                    "dataset": check["dataset"],
                    "claim": check["claim"],
                    "cell_type_v0": ct,
                    "score": score,
                    "positive_group": check["positive_group"],
                    "negative_group": check["negative_group"],
                    "positive_mean": a,
                    "negative_mean": b,
                    "diff_positive_minus_negative": diff,
                    "direction_flag": support,
                })

    return pd.DataFrame(rows)


def flatten_columns(df):
    new_cols = []
    for c in df.columns:
        if isinstance(c, tuple):
            c2 = "__".join([str(x) for x in c if str(x) != ""])
            new_cols.append(c2)
        else:
            new_cols.append(c)
    df.columns = new_cols
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(project_root()))
    ap.add_argument("--datasets", nargs="*", default=[
        "CNP0005824",
        "GSE166489",
        "GSE167029",
        "GSE183716",
    ])
    args = ap.parse_args()

    root = Path(args.root)
    outdir = root / "results/v0_summary"
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_all_myeloid(root, args.datasets)
    df.to_csv(outdir / "v0_all_myeloid_pseudobulk_merged.tsv", sep="\t", index=False)

    means = group_means(df)
    means = flatten_columns(means)
    means.to_csv(outdir / "v0_myeloid_group_means.tsv", sep="\t", index=False)

    tests = pairwise_tests(df)
    tests.to_csv(outdir / "v0_pairwise_tests.tsv", sep="\t", index=False)

    flags = build_boundary_flags(means, tests)
    flags.to_csv(outdir / "v0_boundary_flags.tsv", sep="\t", index=False)

    paired = gse166489_paired_acute_recovery(df)
    if isinstance(paired, tuple):
        detail, summary = paired
        detail.to_csv(outdir / "GSE166489_paired_acute_recovery_detail.tsv", sep="\t", index=False)
        summary.to_csv(outdir / "GSE166489_paired_acute_recovery_summary.tsv", sep="\t", index=False)
    else:
        paired.to_csv(outdir / "GSE166489_paired_acute_recovery_detail.tsv", sep="\t", index=False)

    print(f"[done] summary outputs -> {outdir}")

    print("\n[boundary flags: Kyn_AHR_myeloid_score]")
    show = flags[flags["score"] == "Kyn_AHR_myeloid_score"].copy()
    print(show.to_string(index=False))


if __name__ == "__main__":
    main()
