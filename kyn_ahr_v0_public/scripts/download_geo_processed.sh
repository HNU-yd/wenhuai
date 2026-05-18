#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-$(cd "$(dirname "$0")/.." && pwd)}
RAW="${ROOT}/raw"
mkdir -p "${RAW}"/{GSE167029,GSE166489,GSE183716}

fetch() {
  local url="$1"
  local out="$2"
  if command -v aria2c >/dev/null 2>&1; then
    aria2c -x 8 -s 8 -c -o "$(basename "$out")" -d "$(dirname "$out")" "$url"
  else
    wget -c "$url" -O "$out"
  fi
}

# Processed GEO files; enough for V0. SRA FASTQ is not needed for the baseline.
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167029/suppl/GSE167029_RAW.tar" \
      "${RAW}/GSE167029/GSE167029_RAW.tar"

# Optional: the all-sample Seurat object is very large. Set DOWNLOAD_RDS=1 if you want it.
if [[ "${DOWNLOAD_RDS:-0}" == "1" ]]; then
  fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167029/suppl/GSE167029_SeuratObject_SC_AllSamples.rds.gz" \
        "${RAW}/GSE167029/GSE167029_SeuratObject_SC_AllSamples.rds.gz"
fi

fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166489/suppl/GSE166489_RAW.tar" \
      "${RAW}/GSE166489/GSE166489_RAW.tar"

for i in 1 2 3 4; do
  fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE183nnn/GSE183716/suppl/GSE183716_Sample${i}_GEXFB_filtered_feature_bc_matrix.h5" \
        "${RAW}/GSE183716/GSE183716_Sample${i}_GEXFB_filtered_feature_bc_matrix.h5"
done

if [[ "${EXTRACT:-1}" == "1" ]]; then
  for f in "${RAW}"/GSE167029/*.tar "${RAW}"/GSE166489/*.tar; do
    [[ -f "$f" ]] || continue
    outdir="$(dirname "$f")/extracted"
    mkdir -p "$outdir"
    tar -xvf "$f" -C "$outdir"
  done
fi

echo "Downloaded to: ${RAW}"
