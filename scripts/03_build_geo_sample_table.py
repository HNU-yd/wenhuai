#!/usr/bin/env python
import argparse
import gzip
import shutil
from pathlib import Path

import sys
from project_paths import add_src_to_path, project_path

add_src_to_path()

from kyn_ahr_v0.metadata_rules import build_geo_sample_table


def gunzip_soft_if_needed(data_root: Path):
    for f in data_root.rglob("*family.soft.gz"):
        out = f.with_suffix("")
        if out.exists() and out.stat().st_size > 0:
            continue

        print(f"[gunzip soft] {f} -> {out}")
        with gzip.open(f, "rb") as fin, out.open("wb") as fout:
            shutil.copyfileobj(fin, fout)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data_root", default=str(project_path("data")))
    ap.add_argument("--out", default=str(project_path("data", "metadata", "geo_sample_table.tsv")))
    args = ap.parse_args()

    data_root = Path(args.data_root)
    gunzip_soft_if_needed(data_root)

    df = build_geo_sample_table(str(data_root))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"[done] GEO sample table: {args.out}")
    print(f"[rows] {len(df)}")
    if len(df):
        print(df.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
