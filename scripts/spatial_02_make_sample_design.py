#!/usr/bin/env python
from pathlib import Path
import pandas as pd

from project_paths import project_path

INV = project_path("results", "spatial_inventory", "STT0000127_preferred_gem_profile.tsv")
OUT = project_path("results", "spatial_inventory", "STT0000127_spatial_sample_design.tsv")


def order_stage(sample_id):
    if sample_id == "Control":
        return 0

    import re
    m = re.search(r"d(\d+)", str(sample_id))
    if m:
        return int(m.group(1))

    return 999


def main():
    df = pd.read_csv(INV, sep="\t")

    df["stage_order"] = df["sample_id"].map(order_stage)

    # profile 文件里路径列叫 path，这里统一改成 preferred_path。
    if "preferred_path" not in df.columns:
        df["preferred_path"] = df["path"]

    keep = [
        "sample_id",
        "group",
        "stage",
        "stage_order",
        "preferred_role",
        "preferred_path",
        "size_gb",
        "n_records_scanned",
        "unique_genes_scanned",
        "x_min",
        "x_max",
        "y_min",
        "y_max",
        "total_counts_scanned",
        "mean_counts_per_record",
    ]

    for c in keep:
        if c not in df.columns:
            df[c] = ""

    out = df[keep].sort_values(["group", "stage_order", "sample_id"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False)

    print(f"[done] {OUT}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
