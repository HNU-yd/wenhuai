#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/data1/wenhuai/GNN"
LOG="${ROOT}/logs"

mkdir -p "${ROOT}" "${LOG}"

download_geo_series () {
  local GSE="$1"
  local PREFIX="$2"

  local OUT="${ROOT}/GEO/${GSE}"
  mkdir -p "${OUT}/supp" "${OUT}/matrix" "${OUT}/soft" "${OUT}/logs"

  echo "============================================================"
  echo "[GEO] Download ${GSE}"
  echo "============================================================"

  echo "[1] Supplementary / processed files: ${GSE}"
  wget -e robots=off \
    -r -np -nH --cut-dirs=5 \
    -R "index.html*,*.png,*.gif,*.svg,*.css,*.js" \
    --reject-regex "/assets/" \
    --tries=0 --timeout=60 --waitretry=5 \
    -c \
    -P "${OUT}/supp" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${GSE}/suppl/" \
    2>&1 | tee "${OUT}/logs/download_${GSE}_supp.log"

  echo "[2] Series matrix files: ${GSE}"
  wget -e robots=off \
    -r -np -nH --cut-dirs=5 \
    -R "index.html*,*.png,*.gif,*.svg,*.css,*.js" \
    --tries=0 --timeout=60 --waitretry=5 \
    -c \
    -P "${OUT}/matrix" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${GSE}/matrix/" \
    2>&1 | tee "${OUT}/logs/download_${GSE}_matrix.log" || true

  echo "[3] SOFT metadata file: ${GSE}"
  wget -c \
    -O "${OUT}/soft/${GSE}_family.soft.gz" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${GSE}/soft/${GSE}_family.soft.gz" \
    2>&1 | tee "${OUT}/logs/download_${GSE}_soft.log" || true

  echo "[DONE] ${GSE}"
}

mkdir -p "${ROOT}/GEO"

# 1. 儿童 MIS-C with severe myocarditis PBMC scRNA-seq
download_geo_series "GSE167029" "GSE167nnn"

# 2. MIS-C severe/moderate/recovered + pediatric healthy, GEX/TCR/BCR/ADT
download_geo_series "GSE166489" "GSE166nnn"

# 3. MIS-C acute to recovery PBMC scRNA-seq/CITE-seq
download_geo_series "GSE183716" "GSE183nnn"

# 4. ICI myocarditis PBMC scRNA-seq/TCR/feature barcoding
download_geo_series "GSE180045" "GSE180nnn"

echo "============================================================"
echo "[CNGB] Download CNP0005824 snRNA-seq"
echo "============================================================"

mkdir -p "${ROOT}/CNP0005824_snRNA"

wget -e robots=off \
  -r -np -nH --cut-dirs=4 \
  -R "index.html*,*.png,*.gif,*.svg,*.css,*.js" \
  --reject-regex "/assets/" \
  --tries=0 --timeout=60 --waitretry=5 \
  -c \
  -P "${ROOT}/CNP0005824_snRNA" \
  "https://ftp.cngb.org/pub/CNSA/data5/CNP0005824/" \
  2>&1 | tee "${LOG}/download_CNP0005824.log"

echo "============================================================"
echo "[STOmics] Download STT0000127 Stereo-seq"
echo "============================================================"

mkdir -p "${ROOT}/STT0000127_stereo"

wget -e robots=off \
  -r -np -nH --cut-dirs=3 \
  -R "index.html*,*.png,*.gif,*.svg,*.css,*.js" \
  --reject-regex "/assets/" \
  --tries=0 --timeout=60 --waitretry=5 \
  -c \
  -P "${ROOT}/STT0000127_stereo" \
  "https://ftp.cngb.org/pub/stomics/STT0000127/" \
  2>&1 | tee "${LOG}/download_STT0000127.log"

echo "============================================================"
echo "[SUMMARY] Disk usage"
echo "============================================================"

du -sh \
  "${ROOT}/GEO/GSE167029" \
  "${ROOT}/GEO/GSE166489" \
  "${ROOT}/GEO/GSE183716" \
  "${ROOT}/GEO/GSE180045" \
  "${ROOT}/CNP0005824_snRNA" \
  "${ROOT}/STT0000127_stereo" || true

echo "[DONE] All downloads finished."
