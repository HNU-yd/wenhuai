#!/usr/bin/env bash
set -euo pipefail

ROOT="${WENHUAI_MR_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
OUT="${ROOT}/results/full1400/accessions/chen1400_accessions.txt"
mkdir -p "$(dirname "$OUT")"

: > "$OUT"

for n in $(seq 90199621 90201020); do
  printf "GCST%08d\n" "$n" >> "$OUT"
done

echo "[DONE] accession list:"
wc -l "$OUT"
head "$OUT"
tail "$OUT"
