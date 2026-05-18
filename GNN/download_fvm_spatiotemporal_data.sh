#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/data1/wenhuai/GNN"

mkdir -p \
  "${ROOT}/supplement" \
  "${ROOT}/CNP0005824_snRNA" \
  "${ROOT}/STT0000127_stereo" \
  "${ROOT}/logs"

cd "${ROOT}"

echo "[1/4] Download supplementary tables and figures from Nature/Springer..."

wget -c -O "${ROOT}/supplement/41392_2025_2143_MOESM2_tables_1_10.xlsx" \
  "https://static-content.springer.com/esm/art%3A10.1038%2Fs41392-025-02143-9/MediaObjects/41392_2025_2143_MOESM2_ESM.xlsx"

wget -c -O "${ROOT}/supplement/41392_2025_2143_MOESM1_figures_1_14.pdf" \
  "https://static-content.springer.com/esm/art%3A10.1038%2Fs41392-025-02143-9/MediaObjects/41392_2025_2143_MOESM1_ESM.pdf"

echo "[2/4] Probe CNSA public FTP for CNP0005824..."

wget -qO "${ROOT}/logs/CNSA_data6_index.html" \
  "https://ftp.cngb.org/pub/CNSA/data6/" || true

if grep -q "CNP0005824" "${ROOT}/logs/CNSA_data6_index.html"; then
  echo "[FOUND] CNP0005824 visible in public FTP. Start downloading..."
  wget -r -np -nH --cut-dirs=4 \
    -R "index.html*" \
    -c \
    -P "${ROOT}/CNP0005824_snRNA" \
    "https://ftp.cngb.org/pub/CNSA/data6/CNP0005824/"
else
  echo "[WARN] CNP0005824 is NOT visible in public CNSA FTP data6."
  echo "[WARN] You likely need to open CNGB/CNSA project page, login if needed, then copy the official download command."
fi

echo "[3/4] Probe STOmics public FTP for STT0000127..."

if wget --spider -q "https://ftp.cngb.org/pub/stomics/STT0000127/"; then
  echo "[FOUND] STT0000127 public FTP directory exists. Start downloading..."
  wget -r -np -nH --cut-dirs=3 \
    -R "index.html*" \
    -c \
    -P "${ROOT}/STT0000127_stereo" \
    "https://ftp.cngb.org/pub/stomics/STT0000127/"
else
  echo "[WARN] STT0000127 FTP directory was not directly visible by wget."
  echo "[WARN] Open STOmicsDB project STT0000127 and copy the official download command."
fi

echo "[4/4] Current downloaded files:"
find "${ROOT}" -maxdepth 3 -type f | sed "s#${ROOT}/##" | head -100

echo "[DONE]"
