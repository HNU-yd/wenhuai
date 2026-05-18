#!/usr/bin/env python
import argparse
import gzip
import re
from pathlib import Path

import pandas as pd

from project_paths import project_path


def classify_gem(path: Path) -> str:
    name = path.name

    if name.startswith("Final.") and "TissueCut" in name and name.endswith(".gem.gz"):
        return "tissuecut_gem"

    if name.endswith(".bin50.gem.gz"):
        return "bin50_gem"

    if name.endswith(".gem.gz"):
        return "raw_gem"

    return "other"


def infer_sample_id(path: Path) -> str:
    name = path.name

    # Final.Control.TissueCut.gem.gz -> Control
    m = re.match(r"Final\.(.+?)\.TissueCut\.gem\.gz$", name)
    if m:
        return m.group(1)

    # Control.bin50.gem.gz -> Control
    m = re.match(r"(.+?)\.bin50\.gem\.gz$", name)
    if m:
        return m.group(1)

    # Control.gem.gz -> Control
    m = re.match(r"(.+?)\.gem\.gz$", name)
    if m:
        return m.group(1)

    return name


def infer_group(sample_id: str) -> str:
    if sample_id == "Control":
        return "control"

    if sample_id.startswith("CVB3"):
        return "CVB3_myocarditis"

    if sample_id.startswith("IVIG"):
        return "IVIG_treated"

    return "unknown"


def infer_stage(sample_id: str) -> str:
    if sample_id == "Control":
        return "control"

    m = re.search(r"d(\d+)", sample_id)
    if m:
        return f"day{m.group(1)}"

    return "unknown"


def read_first_lines(path: Path, n: int = 20):
    lines = []
    try:
        with gzip.open(path, "rt", errors="ignore") as f:
            for i, line in enumerate(f):
                if i >= n:
                    break
                lines.append(line.rstrip("\n"))
    except Exception as e:
        lines.append(f"[ERROR] {repr(e)}")

    return lines


def guess_separator(line: str):
    if "\t" in line:
        return "\t"
    if "," in line:
        return ","
    return None


def split_line(line: str):
    sep = guess_separator(line)
    if sep is None:
        return line.split()
    return line.split(sep)


def guess_columns(lines):
    for line in lines:
        if not line or line.startswith("#"):
            continue

        fields = split_line(line)
        low = [x.lower() for x in fields]

        # header-like
        if any("gene" in x for x in low) and any(x in {"x", "y"} or x.endswith("x") or x.endswith("y") for x in low):
            return fields

        if len(fields) >= 4:
            # likely no header
            return ["gene", "x", "y", "count"]

    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stt_root", default=str(project_path("data", "STT0000127")))
    ap.add_argument("--out_dir", default=str(project_path("results", "spatial_inventory")))
    ap.add_argument("--peek_lines", type=int, default=30)
    args = ap.parse_args()

    stt_root = Path(args.stt_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []

    gem_files = sorted(stt_root.rglob("*.gem.gz"))

    for f in gem_files:
        sample_id = infer_sample_id(f)
        file_role = classify_gem(f)
        first_lines = read_first_lines(f, n=args.peek_lines)
        guessed_cols = guess_columns(first_lines)

        rows.append({
            "dataset": "STT0000127",
            "sample_id": sample_id,
            "group": infer_group(sample_id),
            "stage": infer_stage(sample_id),
            "file_role": file_role,
            "path": str(f),
            "filename": f.name,
            "size_bytes": f.stat().st_size,
            "size_gb": round(f.stat().st_size / 1024**3, 5),
            "guessed_columns": "|".join(guessed_cols),
            "first_nonempty_line": next((x for x in first_lines if x.strip()), ""),
        })

        peek_file = out_dir / f"peek_{sample_id}_{file_role}.txt"
        with peek_file.open("w") as w:
            w.write("\n".join(first_lines) + "\n")

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "STT0000127_gem_inventory.tsv", sep="\t", index=False)

    # 每个 sample 的文件矩阵
    if len(df):
        matrix = (
            df.pivot_table(
                index=["sample_id", "group", "stage"],
                columns="file_role",
                values="path",
                aggfunc="first",
            )
            .reset_index()
        )
        matrix.to_csv(out_dir / "STT0000127_sample_file_matrix.tsv", sep="\t", index=False)

        # 推荐优先使用 tissuecut，其次 raw，最后 bin50。
        pref_rows = []
        for _, row in matrix.iterrows():
            preferred = ""
            preferred_role = ""

            for role in ["tissuecut_gem", "raw_gem", "bin50_gem"]:
                if role in matrix.columns and pd.notna(row.get(role, None)):
                    preferred = row[role]
                    preferred_role = role
                    break

            pref_rows.append({
                "sample_id": row["sample_id"],
                "group": row["group"],
                "stage": row["stage"],
                "preferred_role": preferred_role,
                "preferred_path": preferred,
            })

        pd.DataFrame(pref_rows).to_csv(
            out_dir / "STT0000127_preferred_gem_files.tsv",
            sep="\t",
            index=False,
        )

    print(f"[done] gem files: {len(gem_files)}")
    print(f"[out] {out_dir}/STT0000127_gem_inventory.tsv")
    print(f"[out] {out_dir}/STT0000127_sample_file_matrix.tsv")
    print(f"[out] {out_dir}/STT0000127_preferred_gem_files.tsv")


if __name__ == "__main__":
    main()
