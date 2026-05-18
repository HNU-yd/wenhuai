#!/usr/bin/env python
from pathlib import Path
import argparse
import pandas as pd

from project_paths import project_root


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(project_root()))
    ap.add_argument("--datasets", nargs="*", default=[
        "GSE183716",
        "CNP0005824",
        "GSE166489",
        "GSE167029",
    ])
    args = ap.parse_args()

    root = Path(args.root)

    print("==== h5ad_per_sample ====")
    for ds in args.datasets:
        n = len(list((root / "data/h5ad_per_sample" / ds).glob("*.h5ad")))
        print(f"{ds}\t{n}")

    print("\n==== merged h5ad ====")
    for ds in args.datasets:
        f = root / "data/h5ad_merged" / f"{ds}.raw_merged.h5ad"
        print(f"{ds}\t{f.exists()}\t{f.stat().st_size if f.exists() else 0}\t{f}")

    print("\n==== v0 final ====")
    for ds in args.datasets:
        f = root / "results/v0" / ds / f"{ds}.final_v0.h5ad"
        print(f"{ds}\t{f.exists()}\t{f.stat().st_size if f.exists() else 0}\t{f}")

    print("\n==== Kyn-AHR high summaries ====")
    for ds in args.datasets:
        f = root / "results/v0" / ds / f"{ds}.top10pct_Kyn_AHR_high_summary.tsv"
        if not f.exists():
            print(f"[missing] {ds}: {f}")
            continue
        print(f"\n## {ds}")
        df = pd.read_csv(f, sep="\t")
        print(df.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
