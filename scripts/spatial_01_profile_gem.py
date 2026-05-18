#!/usr/bin/env python
import argparse
import gzip
import math
from pathlib import Path

import pandas as pd

from project_paths import project_path


def guess_sep(line: str):
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return None


def split_line(line: str):
    sep = guess_sep(line)
    if sep is None:
        return line.split()
    return line.rstrip("\n").split(sep)


def is_header(fields):
    low = [x.lower() for x in fields]
    if any("gene" in x for x in low) and any("x" == x or x.endswith("x") for x in low):
        return True
    if any("mid" in x or "count" in x or "umi" in x for x in low):
        return True
    return False


def infer_col_indices(header_fields):
    low = [x.lower() for x in header_fields]

    gene_i = None
    x_i = None
    y_i = None
    count_i = None

    for i, c in enumerate(low):
        if gene_i is None and ("gene" in c):
            gene_i = i
        if x_i is None and c in {"x", "xcoord", "x_coord", "xpos", "x_pos"}:
            x_i = i
        if y_i is None and c in {"y", "ycoord", "y_coord", "ypos", "y_pos"}:
            y_i = i
        if count_i is None and ("count" in c or "mid" in c or "umi" in c):
            count_i = i

    # GEM 常见无 header 顺序：geneID x y MIDCounts
    if gene_i is None:
        gene_i = 0
    if x_i is None:
        x_i = 1
    if y_i is None:
        y_i = 2
    if count_i is None:
        count_i = 3

    return gene_i, x_i, y_i, count_i


def profile_one(path: Path, max_records: int):
    n_records = 0
    n_bad = 0
    genes = set()

    x_min = math.inf
    x_max = -math.inf
    y_min = math.inf
    y_max = -math.inf
    total_counts = 0.0

    header = None
    col_idx = None
    first_data_line = ""

    with gzip.open(path, "rt", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = split_line(line)
            if len(fields) < 4:
                continue

            if header is None and is_header(fields):
                header = fields
                col_idx = infer_col_indices(header)
                continue

            if col_idx is None:
                header = ["gene", "x", "y", "count"]
                col_idx = infer_col_indices(header)

            gene_i, x_i, y_i, count_i = col_idx

            try:
                gene = fields[gene_i]
                x = int(float(fields[x_i]))
                y = int(float(fields[y_i]))
                count = float(fields[count_i])
            except Exception:
                n_bad += 1
                continue

            if not first_data_line:
                first_data_line = line

            genes.add(gene)
            x_min = min(x_min, x)
            x_max = max(x_max, x)
            y_min = min(y_min, y)
            y_max = max(y_max, y)
            total_counts += count
            n_records += 1

            if max_records > 0 and n_records >= max_records:
                break

    return {
        "path": str(path),
        "filename": path.name,
        "n_records_scanned": n_records,
        "n_bad_lines": n_bad,
        "unique_genes_scanned": len(genes),
        "x_min": None if x_min is math.inf else x_min,
        "x_max": None if x_max == -math.inf else x_max,
        "y_min": None if y_min is math.inf else y_min,
        "y_max": None if y_max == -math.inf else y_max,
        "total_counts_scanned": total_counts,
        "mean_counts_per_record": total_counts / n_records if n_records else None,
        "header": "|".join(header) if header else "",
        "first_data_line": first_data_line,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preferred_table", default=str(project_path("results", "spatial_inventory", "STT0000127_preferred_gem_files.tsv")))
    ap.add_argument("--out", default=str(project_path("results", "spatial_inventory", "STT0000127_preferred_gem_profile.tsv")))
    ap.add_argument("--max_records", type=int, default=2000000)
    args = ap.parse_args()

    pref = pd.read_csv(args.preferred_table, sep="\t")
    rows = []

    for _, r in pref.iterrows():
        path = Path(r["preferred_path"])
        if not path.exists():
            print(f"[skip missing] {path}")
            continue

        print(f"[profile] {r['sample_id']} {r['preferred_role']} {path}")
        prof = profile_one(path, max_records=args.max_records)
        prof.update({
            "sample_id": r["sample_id"],
            "group": r["group"],
            "stage": r["stage"],
            "preferred_role": r["preferred_role"],
            "size_gb": round(path.stat().st_size / 1024**3, 5),
        })
        rows.append(prof)

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"[done] {args.out}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
