#!/usr/bin/env bash
set -euo pipefail

ROOT="${WENHUAI_MR_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
ACC_LIST="${ROOT}/results/full1400/accessions/chen1400_accessions.txt"
OUT="${ROOT}/exposure_full"
LOGDIR="${ROOT}/logs/full1400/download"
JOBS="${1:-16}"

mkdir -p "$OUT" "$LOGDIR"

range_for_acc () {
  local acc="$1"
  local n="${acc#GCST}"
  local start=$(( ((10#$n - 1) / 1000) * 1000 + 1 ))
  local end=$(( start + 999 ))
  printf "GCST%08d-GCST%08d" "$start" "$end"
}

download_one () {
  local acc="$1"
  local root="${ROOT}"
  local out="${root}/exposure_full/${acc}"
  local log="${root}/logs/full1400/download/${acc}.log"
  mkdir -p "$out"

  local range
  range=$(range_for_acc "$acc")

  local base="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/${range}/${acc}"

  local final="${out}/${acc}_buildGRCh38.tsv.gz"

  if [[ -s "$final" ]] && gzip -t "$final" >/dev/null 2>&1; then
    echo "[SKIP OK] $acc" > "$log"
    return 0
  fi

  rm -f "${out}/${acc}"*.tsv.gz "${out}/${acc}"*.tsv "${out}/index.html" "${out}/robots.txt" 2>/dev/null || true

  {
    echo "[ACC] $acc"
    echo "[RANGE] $range"

    for fname in \
      "${acc}_buildGRCh38.tsv.gz" \
      "${acc}.tsv.gz" \
      "${acc}.h.tsv.gz" \
      "${acc}_buildGRCh37.tsv.gz"
    do
      url="${base}/${fname}"
      echo "[TRY] $url"

      if wget --spider -q "$url"; then
        wget -c -O "${out}/${fname}" "$url"

        if [[ "$fname" == *.gz ]]; then
          gzip -t "${out}/${fname}"
        fi

        if [[ "$fname" != "${acc}_buildGRCh38.tsv.gz" ]]; then
          ln -sf "${fname}" "$final"
        fi

        echo "[DONE] $acc -> ${out}/${fname}"
        exit 0
      fi
    done

    echo "[FAILED] $acc no known file matched"
    exit 2
  } > "$log" 2>&1
}

export -f download_one
export -f range_for_acc

echo "[START] jobs=$JOBS"
cat "$ACC_LIST" | xargs -I{} -P "$JOBS" bash -lc 'download_one "$@"' _ {}

echo "[DONE] downloaded:"
find "$OUT" -name "*.tsv.gz" | wc -l

echo "[FAILED logs]"
grep -L "\[DONE\]\|\[SKIP OK\]" "$LOGDIR"/*.log 2>/dev/null || true
