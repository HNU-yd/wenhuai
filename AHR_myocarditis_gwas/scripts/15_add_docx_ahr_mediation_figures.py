#!/usr/bin/env python3
"""Build the additional AHR mediation-MR deliverables requested in add-project docx.

The local repository contains AHR cis-eQTL data and myocarditis locus data, so the
two requested AHR instruments can be evaluated against the myocarditis outcomes.
This script intentionally limits the add-project outputs to the two SNPs named in
the docx request: rs17643734 and rs59291726.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch


PROJECT = Path(__file__).resolve().parents[1]

TARGET_SNPS = {
    "rs17643734": {"chr": 7, "pos": 17162882, "ea": "G", "oa": "A"},
    "rs59291726": {"chr": 7, "pos": 17211078, "ea": "T", "oa": "C"},
}

OUTCOMES = [
    {
        "dataset": "Sakaue2021_BBJ_Myocarditis",
        "population": "East Asian / BBJ",
        "path": PROJECT
        / "data/locus/outcome/Sakaue2021_BBJ_Myocarditis.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz",
    },
    {
        "dataset": "Sakaue2021_EUR_Myocarditis",
        "population": "European",
        "path": PROJECT
        / "data/locus/outcome/Sakaue2021_EUR_Myocarditis.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz",
    },
    {
        "dataset": "FinnGen_R12_I9_MYOCARD",
        "population": "Finnish / European",
        "path": PROJECT
        / "data/locus/outcome/FinnGen_R12_I9_MYOCARD.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz",
    },
]

EXPOSURE_PATH = PROJECT / "data/formatted/exposure/eqtlgen_AHR_full_cis_eqtl.standardized.tsv.gz"
COLOC_MAIN_PATH = PROJECT / "results/coloc/AHR_eqtlgen_coloc_main_result.tsv"

ADD_RESULTS = PROJECT / "results/add_project"
FIG_DIR = PROJECT / "results/figures"
DOCS_DIR = PROJECT / "docs"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def as_num(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def normal_allele(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def complement(allele: str) -> str:
    table = str.maketrans("ACGT", "TGCA")
    return allele.translate(table)


def align_outcome_beta(
    beta_out: float,
    out_ea: Any,
    out_oa: Any,
    exp_ea: Any,
    exp_oa: Any,
) -> tuple[float, str]:
    """Align outcome beta to the exposure effect allele."""
    oe = normal_allele(out_ea)
    oo = normal_allele(out_oa)
    ee = normal_allele(exp_ea)
    eo = normal_allele(exp_oa)

    if oe == ee and oo == eo:
        return beta_out, "aligned"
    if oe == eo and oo == ee:
        return -beta_out, "swapped"

    coe = complement(oe)
    coo = complement(oo)
    if coe == ee and coo == eo:
        return beta_out, "complement_aligned"
    if coe == eo and coo == ee:
        return -beta_out, "complement_swapped"

    return math.nan, "allele_mismatch"


def two_sided_p_from_z(z: float) -> float:
    if not math.isfinite(z):
        return math.nan
    return math.erfc(abs(z) / math.sqrt(2.0))


def neglog10(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    numeric = numeric.mask(numeric <= 0, 1e-300)
    return -np.log10(numeric)


def find_outcome_row(
    outcome_df: pd.DataFrame,
    snp: str,
    exp_row: pd.Series,
) -> tuple[pd.Series | None, str]:
    """Find the requested variant by rsID first, then by chr:pos and alleles."""
    if "SNP" in outcome_df.columns:
        match = outcome_df[outcome_df["SNP"].astype(str) == snp]
        if not match.empty:
            return match.sort_values("pval").iloc[0], "rsID"

    exp_chr = int(as_num(exp_row["chr"]))
    exp_pos = int(as_num(exp_row["pos"]))
    chr_num = pd.to_numeric(outcome_df["chr"], errors="coerce")
    pos_num = pd.to_numeric(outcome_df["pos"], errors="coerce")
    match = outcome_df[(chr_num == exp_chr) & (pos_num == exp_pos)]

    if not match.empty:
        exp_alleles = {normal_allele(exp_row["effect_allele"]), normal_allele(exp_row["other_allele"])}
        allele_match = match[
            match.apply(
                lambda row: {
                    normal_allele(row.get("effect_allele")),
                    normal_allele(row.get("other_allele")),
                }
                == exp_alleles,
                axis=1,
            )
        ]
        if not allele_match.empty:
            return allele_match.sort_values("pval").iloc[0], "chr_pos_alleles"
        return match.sort_values("pval").iloc[0], "chr_pos"

    return None, "missing"


def build_two_snp_wald(exposure_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    exposure_targets = {
        snp: exposure_df[exposure_df["SNP"].astype(str) == snp].iloc[0]
        for snp in TARGET_SNPS
        if not exposure_df[exposure_df["SNP"].astype(str) == snp].empty
    }

    for outcome in OUTCOMES:
        outcome_df = read_tsv(outcome["path"])
        for snp in TARGET_SNPS:
            exp_row = exposure_targets.get(snp)
            base = {
                "exposure": "AHR_expression_eQTLGen_whole_blood",
                "exposure_gwas_id": "eqtl-a-ENSG00000106546",
                "outcome_dataset": outcome["dataset"],
                "outcome_population": outcome["population"],
                "target_snp": snp,
                "status": "ok",
            }

            if exp_row is None:
                rows.append({**base, "status": "missing_exposure"})
                continue

            beta_exp = as_num(exp_row["beta"])
            se_exp = as_num(exp_row["se"])
            out_row, match_method = find_outcome_row(outcome_df, snp, exp_row)

            row = {
                **base,
                "exp_chr": int(as_num(exp_row["chr"])),
                "exp_pos": int(as_num(exp_row["pos"])),
                "effect_allele": normal_allele(exp_row["effect_allele"]),
                "other_allele": normal_allele(exp_row["other_allele"]),
                "eaf_exposure": as_num(exp_row.get("eaf")),
                "beta_exposure": beta_exp,
                "se_exposure": se_exp,
                "pval_exposure": as_num(exp_row["pval"]),
                "f_stat": as_num(exp_row["f_stat"]),
                "match_method": match_method,
            }

            if out_row is None:
                rows.append({**row, "status": "missing_outcome"})
                continue

            beta_out = as_num(out_row["beta"])
            se_out = as_num(out_row["se"])
            beta_out_aligned, allele_match = align_outcome_beta(
                beta_out,
                out_row.get("effect_allele"),
                out_row.get("other_allele"),
                exp_row.get("effect_allele"),
                exp_row.get("other_allele"),
            )

            ratio = math.nan
            ratio_se = math.nan
            ratio_p = math.nan
            ratio_or = math.nan
            ratio_or_lci95 = math.nan
            ratio_or_uci95 = math.nan
            if math.isfinite(beta_out_aligned) and beta_exp != 0:
                ratio = beta_out_aligned / beta_exp
                ratio_se = math.sqrt((se_out / beta_exp) ** 2 + ((beta_out_aligned * se_exp) / (beta_exp**2)) ** 2)
                ratio_p = two_sided_p_from_z(ratio / ratio_se)
                ratio_or = math.exp(ratio)
                ratio_or_lci95 = math.exp(ratio - 1.96 * ratio_se)
                ratio_or_uci95 = math.exp(ratio + 1.96 * ratio_se)

            if allele_match == "allele_mismatch":
                row["status"] = "allele_mismatch"

            rows.append(
                {
                    **row,
                    "outcome_chr": as_num(out_row.get("chr")),
                    "outcome_pos": as_num(out_row.get("pos")),
                    "outcome_SNP": out_row.get("SNP"),
                    "outcome_variant_id": out_row.get("variant_id"),
                    "beta_outcome_raw": beta_out,
                    "beta_outcome_aligned": beta_out_aligned,
                    "se_outcome": se_out,
                    "pval_outcome": as_num(out_row.get("pval")),
                    "effect_allele_outcome": normal_allele(out_row.get("effect_allele")),
                    "other_allele_outcome": normal_allele(out_row.get("other_allele")),
                    "allele_match": allele_match,
                    "ratio": ratio,
                    "ratio_se": ratio_se,
                    "ratio_pval": ratio_p,
                    "ratio_or": ratio_or,
                    "ratio_or_lci95": ratio_or_lci95,
                    "ratio_or_uci95": ratio_or_uci95,
                }
            )

    return pd.DataFrame(rows)


def build_mediation_status(wald_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    ok_wald = wald_df[wald_df["status"] == "ok"].copy()
    rows.append(
        {
            "mediation_arm": "AHR expression -> Myocarditis",
            "status": "completed_two_snp_wald_where_available",
            "method": "Single-SNP Wald ratio",
            "nsnp": ok_wald["target_snp"].nunique() if not ok_wald.empty else 0,
            "beta": pd.NA,
            "se": pd.NA,
            "pval": pd.NA,
            "or": pd.NA,
            "or_lci95": pd.NA,
            "or_uci95": pd.NA,
            "note": (
                "Docx-specified SNPs rs17643734 and rs59291726 were evaluated against "
                "the three myocarditis GWAS datasets; see AHR_two_significant_snp_wald_ratios.tsv."
            ),
        }
    )

    rows.append(
        {
            "mediation_arm": "AHR eQTL and Myocarditis coloc",
            "status": "completed",
            "method": "coloc.abf",
            "nsnp": len(OUTCOMES),
            "beta": pd.NA,
            "se": pd.NA,
            "pval": pd.NA,
            "or": pd.NA,
            "or_lci95": pd.NA,
            "or_uci95": pd.NA,
            "note": "Coloc summary and SNP-level posterior context were generated for the requested AHR SNPs.",
        }
    )

    return pd.DataFrame(rows)


def build_two_snp_coloc_context() -> pd.DataFrame:
    coloc = read_tsv(COLOC_MAIN_PATH)
    rows: list[pd.DataFrame] = []

    for record in coloc.itertuples(index=False):
        posterior_file = Path(getattr(record, "snp_posterior_file"))
        if not posterior_file.exists():
            continue
        posterior = read_tsv(posterior_file)
        snp_col = "snp" if "snp" in posterior.columns else "exp_SNP"
        selected = posterior[posterior[snp_col].astype(str).isin(TARGET_SNPS)].copy()
        if selected.empty:
            continue
        selected["outcome_dataset"] = getattr(record, "outcome_dataset")
        rows.append(selected)

    if not rows:
        return pd.DataFrame()

    combined = pd.concat(rows, ignore_index=True, sort=False)
    keep = [
        "outcome_dataset",
        "snp",
        "exp_chr",
        "exp_pos",
        "match_method",
        "exp_effect_allele",
        "exp_other_allele",
        "beta_eqtl",
        "se_eqtl",
        "pval_eqtl",
        "beta_gwas_aligned",
        "se_gwas",
        "pval_gwas",
        "SNP.PP.H4",
        "outcome_variant_id",
    ]
    keep = [col for col in keep if col in combined.columns]
    return combined[keep].sort_values(["outcome_dataset", "snp"])


def plot_mediation_model(output: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.set_axis_off()

    nodes = {
        "Upstream exposure\n(existing MR)": (0.17, 0.62),
        "AHR expression\n(eQTLGen whole blood)": (0.50, 0.62),
        "Myocarditis\nGWAS": (0.83, 0.62),
    }

    colors = ["#dbeafe", "#dcfce7", "#fee2e2"]
    for (label, (x, y)), color in zip(nodes.items(), colors):
        ax.text(
            x,
            y,
            label,
            ha="center",
            va="center",
            fontsize=13,
            bbox=dict(boxstyle="round,pad=0.45", facecolor=color, edgecolor="#334155", linewidth=1.3),
        )

    def arrow(start: tuple[float, float], end: tuple[float, float], rad: float = 0.0) -> None:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=16,
                linewidth=1.8,
                color="#1f2937",
                connectionstyle=f"arc3,rad={rad}",
            )
        )

    arrow((0.29, 0.62), (0.39, 0.62))
    arrow((0.61, 0.62), (0.73, 0.62))
    arrow((0.24, 0.52), (0.76, 0.52), rad=0.14)

    ax.text(0.50, 0.76, "Mediator MR framework", ha="center", fontsize=16, weight="bold")
    ax.text(0.50, 0.43, "AHR -> myocarditis tested with rs17643734 and rs59291726", ha="center", fontsize=11)
    ax.text(0.50, 0.35, "Only the two docx-specified AHR SNPs are included in this add-project output", ha="center", fontsize=11)
    ax.text(0.50, 0.27, "Coloc and locus plots summarize AHR eQTL and myocarditis evidence", ha="center", fontsize=11)

    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_exposure_manhattan(exposure_df: pd.DataFrame, output: Path) -> None:
    df = exposure_df.copy()
    df["pos_num"] = pd.to_numeric(df["pos"], errors="coerce")
    df["neglog10p"] = neglog10(df["pval"])
    df = df.dropna(subset=["pos_num", "neglog10p"]).sort_values("pos_num")

    fig, ax = plt.subplots(figsize=(11, 5.8))
    ax.scatter(df["pos_num"] / 1e6, df["neglog10p"], s=10, color="#64748b", alpha=0.55, linewidths=0)
    ax.axhline(-math.log10(5e-8), color="#dc2626", linestyle="--", linewidth=1, label="P=5e-8")

    for snp, color in zip(TARGET_SNPS, ["#b91c1c", "#f97316"]):
        target = df[df["SNP"].astype(str) == snp]
        if target.empty:
            continue
        row = target.iloc[0]
        ax.scatter(row["pos_num"] / 1e6, row["neglog10p"], s=72, color=color, edgecolor="black", linewidth=0.7, zorder=5)
        ax.annotate(
            snp,
            (row["pos_num"] / 1e6, row["neglog10p"]),
            xytext=(6, 8),
            textcoords="offset points",
            fontsize=10,
            color="#111827",
        )

    gene_pos = pd.to_numeric(df.get("GenePos"), errors="coerce").dropna()
    if not gene_pos.empty:
        ax.axvline(gene_pos.iloc[0] / 1e6, color="#2563eb", linestyle=":", linewidth=1.2, label="AHR gene")

    ax.set_title("AHR cis-eQTL locus: eQTLGen whole blood")
    ax.set_xlabel("Chromosome 7 position (Mb)")
    ax.set_ylabel("-log10(P)")
    ax.legend(frameon=False, loc="upper right")
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_outcome_locus(exposure_df: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(len(OUTCOMES), 1, figsize=(11, 9.5), sharex=False)

    if len(OUTCOMES) == 1:
        axes = [axes]

    for ax, outcome in zip(axes, OUTCOMES):
        df = read_tsv(outcome["path"])
        df["pos_num"] = pd.to_numeric(df["pos"], errors="coerce")
        df["neglog10p"] = neglog10(df["pval"])
        df = df.dropna(subset=["pos_num", "neglog10p"]).sort_values("pos_num")

        ax.scatter(df["pos_num"] / 1e6, df["neglog10p"], s=8, color="#64748b", alpha=0.55, linewidths=0)
        ax.axhline(-math.log10(5e-8), color="#dc2626", linestyle="--", linewidth=0.8)

        for snp, color in zip(TARGET_SNPS, ["#b91c1c", "#f97316"]):
            exp = exposure_df[exposure_df["SNP"].astype(str) == snp]
            if exp.empty:
                continue
            target, _ = find_outcome_row(df, snp, exp.iloc[0])
            if target is None:
                continue
            x = as_num(target.get("pos")) / 1e6
            y = -math.log10(max(as_num(target.get("pval")), 1e-300))
            ax.scatter(x, y, s=58, color=color, edgecolor="black", linewidth=0.6, zorder=5)
            ax.annotate(snp, (x, y), xytext=(5, 6), textcoords="offset points", fontsize=9)

        ax.set_title(outcome["dataset"])
        ax.set_ylabel("-log10(P)")

    axes[-1].set_xlabel("Chromosome 7 position in source build (Mb)")
    fig.suptitle("Myocarditis AHR-locus Manhattan plots", y=0.995, fontsize=14, weight="bold")
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_coloc_summary(output: Path) -> None:
    coloc = read_tsv(COLOC_MAIN_PATH)
    x = np.arange(len(coloc))
    width = 0.34

    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    ax.bar(x - width / 2, pd.to_numeric(coloc["PP.H3.abf"], errors="coerce"), width, label="PP.H3", color="#60a5fa")
    ax.bar(x + width / 2, pd.to_numeric(coloc["PP.H4.abf"], errors="coerce"), width, label="PP.H4", color="#f97316")
    ax.axhline(0.8, color="#dc2626", linestyle="--", linewidth=1, label="0.8 reference")
    ax.set_xticks(x)
    ax.set_xticklabels(coloc["outcome_dataset"], rotation=18, ha="right")
    ax.set_ylabel("Posterior probability")
    ax.set_title("AHR eQTL and myocarditis coloc summary")
    ax.set_ylim(0, max(0.85, float(pd.to_numeric(coloc[["PP.H3.abf", "PP.H4.abf"]].stack(), errors="coerce").max()) * 1.2))
    ax.legend(frameon=False)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_wald_forest(wald_df: pd.DataFrame, output: Path) -> None:
    df = wald_df[wald_df["status"] == "ok"].copy()
    df = df.dropna(subset=["ratio_or", "ratio_or_lci95", "ratio_or_uci95"])
    df = df.sort_values(["outcome_dataset", "target_snp"])

    fig, ax = plt.subplots(figsize=(10.8, max(4.8, 0.55 * len(df) + 1.8)))

    if df.empty:
        ax.text(0.5, 0.5, "No available Wald ratio rows", ha="center", va="center")
        ax.set_axis_off()
    else:
        labels = [f"{row.outcome_dataset} | {row.target_snp}" for row in df.itertuples()]
        y = np.arange(len(df))
        x = pd.to_numeric(df["ratio_or"], errors="coerce").to_numpy()
        lo = pd.to_numeric(df["ratio_or_lci95"], errors="coerce").to_numpy()
        hi = pd.to_numeric(df["ratio_or_uci95"], errors="coerce").to_numpy()
        ax.errorbar(x, y, xerr=[x - lo, hi - x], fmt="o", color="#1f2937", ecolor="#64748b", capsize=3)
        ax.axvline(1.0, color="#dc2626", linestyle="--", linewidth=1)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xscale("log")
        ax.set_xlabel("Odds ratio per genetically predicted AHR expression increase")
        ax.set_title("Two-SNP AHR expression -> myocarditis Wald ratios")
        ax.grid(axis="x", linestyle=":", linewidth=0.7, alpha=0.6)

    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def fmt(value: Any, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    try:
        value_float = float(value)
    except (TypeError, ValueError):
        return str(value)
    if value_float == 0:
        return "0"
    if abs(value_float) < 0.001:
        return f"{value_float:.{digits}e}"
    if abs(value_float - round(value_float)) < 1e-12 and abs(value_float) < 1e9:
        return str(int(round(value_float)))
    if abs(value_float) >= 1000:
        return f"{value_float:.{digits}e}"
    return f"{value_float:.{digits}f}"


def markdown_table(df: pd.DataFrame, columns: list[str]) -> str:
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |"
    lines = [header, sep]
    for _, row in df.iterrows():
        vals = []
        for col in columns:
            vals.append(fmt(row.get(col)) if isinstance(row.get(col), (float, int, np.floating)) or pd.isna(row.get(col)) else str(row.get(col)))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def write_report(
    wald_df: pd.DataFrame,
    status_df: pd.DataFrame,
    coloc_snp_df: pd.DataFrame,
) -> None:
    report_path = ADD_RESULTS / "AHR_add_project_report.md"
    coloc = read_tsv(COLOC_MAIN_PATH)

    wald_view = wald_df[
        [
            "outcome_dataset",
            "target_snp",
            "status",
            "match_method",
            "beta_outcome_aligned",
            "se_outcome",
            "pval_outcome",
            "ratio",
            "ratio_se",
            "ratio_pval",
            "ratio_or",
            "ratio_or_lci95",
            "ratio_or_uci95",
        ]
    ].copy()

    coloc_view = coloc[
        [
            "outcome_dataset",
            "n_overlap",
            "PP.H3.abf",
            "PP.H4.abf",
            "coloc_interpretation",
        ]
    ].copy()

    status_view = status_df[["mediation_arm", "status", "method", "nsnp", "beta", "se", "pval", "or", "note"]].copy()
    coloc_snp_view = coloc_snp_df[
        [
            "outcome_dataset",
            "snp",
            "match_method",
            "beta_gwas_aligned",
            "se_gwas",
            "pval_gwas",
            "SNP.PP.H4",
            "outcome_variant_id",
        ]
    ].copy()

    text = f"""# 加项目.docx 补充分析报告

## 完成内容

根据 `加项目.docx`，本次围绕 AHR 表达 GWAS `eqtl-a-ENSG00000106546` 和两个指定显著 SNP（`rs17643734`、`rs59291726`）补充了中介 MR 相关材料、共定位结果汇总和 Manhattan/locus 图。

已完成的本地可计算部分：

- 使用本地 eQTLGen AHR cis-eQTL 表提取两个 SNP 的 AHR 表达效应。
- 在三套心肌炎 GWAS 的 AHR locus 数据中匹配两个 SNP，并计算 AHR expression -> myocarditis 的单 SNP Wald ratio。
- 生成中介模式图、AHR eQTL locus Manhattan 图、心肌炎 AHR locus Manhattan 图、coloc 汇总图和两 SNP forest 图。

说明：本报告只纳入 `加项目.docx` 指定的 `rs17643734` 和 `rs59291726`。其他未在文档中指定的 SNP 排查不属于本次主交付。

## 两 SNP Wald ratio 结果

{markdown_table(wald_view, list(wald_view.columns))}

## 中介 MR 状态

{markdown_table(status_view, list(status_view.columns))}

## Coloc 汇总

{markdown_table(coloc_view, list(coloc_view.columns))}

三套心肌炎数据的 PP.H4 均低于常用强共定位参考阈值，既有结果解释为 `no_strong_colocalization`。

## 指定 SNP 的 coloc posterior

{markdown_table(coloc_snp_view, list(coloc_snp_view.columns))}

注：`SNP.PP.H4` 来自既有 SNP-level posterior 文件，用于显示在 H4 假设下指定 SNP 的相对定位信息；整体共定位证据仍以主结果表的 `PP.H4.abf` 为准。

## 输出文件

- `results/add_project/AHR_two_significant_snp_wald_ratios.tsv`
- `results/add_project/AHR_two_snp_mediation_status.tsv`
- `results/add_project/AHR_two_snp_coloc_posterior_context.tsv`
- `results/figures/AHR_two_snp_mediation_model.png`
- `results/figures/AHR_eqtlgen_AHR_locus_manhattan_two_snp.png`
- `results/figures/AHR_myocarditis_locus_manhattan_two_snp.png`
- `results/figures/AHR_coloc_summary_with_two_snp_context.png`
- `results/figures/AHR_two_snp_wald_ratio_forest.png`

## 复现命令

```bash
conda run -n wenhuai python AHR_myocarditis_gwas/scripts/15_add_docx_ahr_mediation_figures.py
```
"""

    report_path.write_text(text, encoding="utf-8")


def write_docx_task_note() -> None:
    note = """# 加项目.docx 任务摘录

原 docx 要求的核心事项：

- 目的：跑中介 MR。
- AHR GWAS ID：`eqtl-a-ENSG00000106546`。
- NCBI Gene 参考：AHR gene ID `196`，Ensembl ID `ENSG00000106546`。
- OpenGWAS 数据集检索：`https://opengwas.io/datasets/`。
- 指定 SNP：`rs17643734`、`rs59291726`；只跑两个显著 SNP。
- 补充：中介模式图、coloc、Manhattan plot 图。

本地执行说明：

- AHR expression -> Myocarditis 已用本地 AHR cis-eQTL 和心肌炎 locus 数据完成两 SNP Wald ratio。
- 本次 add-project 主交付只纳入文档指定的 `rs17643734` 和 `rs59291726`；其他未在文档中指定的 SNP 排查不写入主报告。
"""

    (DOCS_DIR / "add_project_docx_tasks.md").write_text(note, encoding="utf-8")


def main() -> None:
    ADD_RESULTS.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    DOCS_DIR.mkdir(parents=True, exist_ok=True)

    exposure_df = read_tsv(EXPOSURE_PATH)
    wald_df = build_two_snp_wald(exposure_df)
    status_df = build_mediation_status(wald_df)
    coloc_snp_df = build_two_snp_coloc_context()

    wald_df.to_csv(ADD_RESULTS / "AHR_two_significant_snp_wald_ratios.tsv", sep="\t", index=False)
    status_df.to_csv(ADD_RESULTS / "AHR_two_snp_mediation_status.tsv", sep="\t", index=False)
    coloc_snp_df.to_csv(ADD_RESULTS / "AHR_two_snp_coloc_posterior_context.tsv", sep="\t", index=False)

    plot_mediation_model(FIG_DIR / "AHR_two_snp_mediation_model.png")
    plot_exposure_manhattan(exposure_df, FIG_DIR / "AHR_eqtlgen_AHR_locus_manhattan_two_snp.png")
    plot_outcome_locus(exposure_df, FIG_DIR / "AHR_myocarditis_locus_manhattan_two_snp.png")
    plot_coloc_summary(FIG_DIR / "AHR_coloc_summary_with_two_snp_context.png")
    plot_wald_forest(wald_df, FIG_DIR / "AHR_two_snp_wald_ratio_forest.png")

    write_report(wald_df, status_df, coloc_snp_df)
    write_docx_task_note()

    print("Additional AHR mediation deliverables written to:")
    print(f"  {ADD_RESULTS}")
    print(f"  {FIG_DIR}")


if __name__ == "__main__":
    main()
