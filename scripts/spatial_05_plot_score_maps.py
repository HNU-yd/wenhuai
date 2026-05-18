#!/usr/bin/env python
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from project_paths import project_path


DEFAULT_SCORES = [
    "KYN_metabolism_score_expr",
    "AHR_response_score_expr",
    "AHR_regulon_proxy_score_expr",
    "myeloid_inflammation_score_expr",
    "Kyn_AHR_axis_score_expr",
    "Kyn_AHR_myeloid_score_expr",
    "Kyn_AHR_myeloid_score_z",
]


def load_files(in_dir, samples=None):
    files = sorted(Path(in_dir).glob("*.scored.h5ad"))
    if samples:
        wanted = list(samples)
        file_map = {f.name.split(".")[0]: f for f in files}
        files = [file_map[s] for s in wanted if s in file_map]
    return files


def get_sample_id_from_file(f):
    return f.name.split(".")[0]


def get_obs_array(adata, col):
    return adata.obs[col].to_numpy()


def compute_global_ranges(files, scores, clip_q_low=0.01, clip_q_high=0.99):
    """
    对每个 score，在所有待比较样本上统一计算 vmin / vmax。
    """
    range_rows = []

    for score in scores:
        all_vals = []

        for f in files:
            adata = sc.read_h5ad(f, backed="r")
            if score not in adata.obs.columns:
                continue

            v = adata.obs[score].astype(float).to_numpy()
            v = v[np.isfinite(v)]
            if len(v):
                all_vals.append(v)

        if len(all_vals) == 0:
            range_rows.append({
                "score": score,
                "n_values": 0,
                "vmin": np.nan,
                "vmax": np.nan,
                "clip_q_low": clip_q_low,
                "clip_q_high": clip_q_high,
            })
            continue

        vals = np.concatenate(all_vals)

        if clip_q_low == 0 and clip_q_high == 1:
            vmin = float(np.nanmin(vals))
            vmax = float(np.nanmax(vals))
        else:
            vmin = float(np.quantile(vals, clip_q_low))
            vmax = float(np.quantile(vals, clip_q_high))

        if not np.isfinite(vmin):
            vmin = float(np.nanmin(vals))
        if not np.isfinite(vmax):
            vmax = float(np.nanmax(vals))

        if vmin == vmax:
            vmin = float(np.nanmin(vals))
            vmax = float(np.nanmax(vals))

        if vmin == vmax:
            vmin -= 1e-6
            vmax += 1e-6

        range_rows.append({
            "score": score,
            "n_values": len(vals),
            "vmin": vmin,
            "vmax": vmax,
            "clip_q_low": clip_q_low,
            "clip_q_high": clip_q_high,
        })

    range_df = pd.DataFrame(range_rows)
    range_map = {
        r["score"]: (r["vmin"], r["vmax"])
        for _, r in range_df.iterrows()
        if pd.notna(r["vmin"]) and pd.notna(r["vmax"])
    }
    return range_df, range_map


def compute_shared_xy_limits(adatas):
    xs = []
    ys = []
    for adata in adatas:
        if "x_center" not in adata.obs.columns or "y_center" not in adata.obs.columns:
            continue
        xs.append(adata.obs["x_center"].to_numpy(dtype=float))
        ys.append(adata.obs["y_center"].to_numpy(dtype=float))

    if not xs or not ys:
        return None

    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    return (
        float(np.nanmin(x_all)),
        float(np.nanmax(x_all)),
        float(np.nanmin(y_all)),
        float(np.nanmax(y_all)),
    )


def make_title(adata, sample_id, score, title_mode="sample_only"):
    group = str(adata.obs["group"].iloc[0]) if "group" in adata.obs.columns else ""
    stage = str(adata.obs["stage"].iloc[0]) if "stage" in adata.obs.columns else ""

    if title_mode == "sample_only":
        return sample_id
    elif title_mode == "sample_group_stage":
        return f"{sample_id}\n{group} | {stage}"
    elif title_mode == "full":
        return f"{sample_id}\n{group} | {stage}\n{score}"
    else:
        return sample_id


def plot_single(
    adata,
    score,
    out_png,
    vmin,
    vmax,
    point_size=2.0,
    xlim=None,
    ylim=None,
    cmap="viridis",
    title_mode="full",
):
    if score not in adata.obs.columns:
        print(f"[skip] missing score: {score}")
        return

    sample_id = str(adata.obs["sample_id"].iloc[0]) if "sample_id" in adata.obs.columns else out_png.stem

    x = adata.obs["x_center"].to_numpy(dtype=float)
    y = adata.obs["y_center"].to_numpy(dtype=float)
    val = adata.obs[score].to_numpy(dtype=float)

    finite = np.isfinite(val)
    if finite.sum() == 0:
        print(f"[skip] all nan: {score}")
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    sca = ax.scatter(
        x,
        y,
        c=val,
        s=point_size,
        vmin=vmin,
        vmax=vmax,
        linewidths=0,
        cmap=cmap,
    )
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.axis("off")

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        # 注意 y 轴已反转，因此这里仍然给正常范围
        ax.set_ylim(ylim[::-1])

    ax.set_title(make_title(adata, sample_id, score, title_mode=title_mode))
    fig.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=220)
    plt.close(fig)


def plot_panel(
    adatas,
    sample_ids,
    score,
    out_png,
    vmin,
    vmax,
    point_size=2.0,
    share_axis_limits=True,
    cmap="viridis",
    title_mode="sample_group_stage",
    suptitle=True,
):
    valid = []
    valid_ids = []
    for adata, sid in zip(adatas, sample_ids):
        if score in adata.obs.columns:
            valid.append(adata)
            valid_ids.append(sid)

    if len(valid) == 0:
        print(f"[skip] no valid sample for panel score={score}")
        return

    n = len(valid)

    # 关键修改 1：缩小整体宽度 + 明确减小子图间距
    fig, axes = plt.subplots(
        1,
        n,
        figsize=(3.6 * n + 0.8, 5),   # 原来 5*n 太宽了
        gridspec_kw={"wspace": 0.02}  # 明确控制子图间横向间隔
    )

    if n == 1:
        axes = [axes]

    xy_limits = compute_shared_xy_limits(valid) if share_axis_limits else None
    xlim = None
    ylim = None
    if xy_limits is not None:
        xlim = (xy_limits[0], xy_limits[1])
        ylim = (xy_limits[2], xy_limits[3])

    last_sca = None
    for ax, adata, sid in zip(axes, valid, valid_ids):
        x = adata.obs["x_center"].to_numpy(dtype=float)
        y = adata.obs["y_center"].to_numpy(dtype=float)
        val = adata.obs[score].to_numpy(dtype=float)

        last_sca = ax.scatter(
            x,
            y,
            c=val,
            s=point_size,
            vmin=vmin,
            vmax=vmax,
            linewidths=0,
            cmap=cmap,
        )
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.axis("off")

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim[::-1])

        ax.set_title(make_title(adata, sid, score, title_mode=title_mode), fontsize=13)

    # 关键修改 2：colorbar 更贴近图
    cbar = fig.colorbar(
        last_sca,
        ax=axes,
        fraction=0.022,
        pad=0.008   # 原来更大，这里缩小
    )
    cbar.ax.set_ylabel(score, rotation=270, labelpad=16)

    if suptitle:
        fig.suptitle(score, fontsize=18)

    # 关键修改 3：不要再用 tight_layout，改用手工调 spacing
    fig.subplots_adjust(
        left=0.02,
        right=0.92,   # 给 colorbar 留位置
        top=0.84,
        bottom=0.03,
        wspace=0.03   # 再压缩一点子图间距
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_dir", default=str(project_path("data", "spatial_h5ad", "STT0000127", "bin50_tissuecut_scored")))
    ap.add_argument("--out_dir", default=str(project_path("results", "spatial_maps", "STT0000127", "bin50_tissuecut")))
    ap.add_argument("--scores", nargs="*", default=DEFAULT_SCORES)
    ap.add_argument("--samples", nargs="*", default=None)
    ap.add_argument("--point_size", type=float, default=2.0)
    ap.add_argument("--clip_q_low", type=float, default=0.01)
    ap.add_argument("--clip_q_high", type=float, default=0.99)
    ap.add_argument("--cmap", default="viridis")

    # 输出控制
    ap.add_argument("--plot_single", action="store_true", help="输出单样本图")
    ap.add_argument("--plot_panel", action="store_true", help="输出拼图")
    ap.add_argument("--panel_name", default=None, help="拼图名称，如 Control_CVB3d6_IVIGd6")
    ap.add_argument("--share_axis_limits", action="store_true", help="拼图时共享坐标范围")
    ap.add_argument("--title_mode_single", default="full", choices=["sample_only", "sample_group_stage", "full"])
    ap.add_argument("--title_mode_panel", default="sample_group_stage", choices=["sample_only", "sample_group_stage", "full"])

    args = ap.parse_args()

    # 默认两个都画
    if not args.plot_single and not args.plot_panel:
        args.plot_single = True
        args.plot_panel = True

    files = load_files(args.in_dir, args.samples)
    if len(files) == 0:
        print("[error] no input files found")
        return

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sample_ids = [get_sample_id_from_file(f) for f in files]
    print("[samples]", sample_ids)

    # 统一色条范围
    range_df, range_map = compute_global_ranges(
        files=files,
        scores=args.scores,
        clip_q_low=args.clip_q_low,
        clip_q_high=args.clip_q_high,
    )
    range_tsv = out_dir / "global_color_ranges.tsv"
    range_df.to_csv(range_tsv, sep="\t", index=False)
    print(f"[done] global ranges -> {range_tsv}")
    print(range_df.to_string(index=False))

    # 预先载入 adata，避免同一轮重复读太多次
    adata_map = {}
    for f, sid in zip(files, sample_ids):
        print(f"[read] {sid}")
        adata_map[sid] = sc.read_h5ad(f)

    # 单图
    if args.plot_single:
        for sid in sample_ids:
            adata = adata_map[sid]
            sample_out = out_dir / sid
            for score in args.scores:
                if score not in range_map:
                    print(f"[skip] no global range for {score}")
                    continue
                vmin, vmax = range_map[score]
                out_png = sample_out / f"{sid}.{score}.png"
                plot_single(
                    adata=adata,
                    score=score,
                    out_png=out_png,
                    vmin=vmin,
                    vmax=vmax,
                    point_size=args.point_size,
                    xlim=None,
                    ylim=None,
                    cmap=args.cmap,
                    title_mode=args.title_mode_single,
                )

    # 拼图
    if args.plot_panel:
        panel_name = args.panel_name if args.panel_name else "__".join(sample_ids)
        panel_out_dir = out_dir / "_panels" / panel_name
        adatas = [adata_map[sid] for sid in sample_ids]

        for score in args.scores:
            if score not in range_map:
                print(f"[skip] no global range for {score}")
                continue
            vmin, vmax = range_map[score]
            out_png = panel_out_dir / f"{panel_name}.{score}.panel.png"
            plot_panel(
                adatas=adatas,
                sample_ids=sample_ids,
                score=score,
                out_png=out_png,
                vmin=vmin,
                vmax=vmax,
                point_size=args.point_size,
                share_axis_limits=args.share_axis_limits,
                cmap=args.cmap,
                title_mode=args.title_mode_panel,
                suptitle=True,
            )

    print(f"[done] maps -> {out_dir}")


if __name__ == "__main__":
    main()
