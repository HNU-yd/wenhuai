#!/usr/bin/env bash
set -euo pipefail

BASE="${AHR_GWAS_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
LOG="${BASE}/logs/01_download_outcome_gwas.log"

mkdir -p "${BASE}/logs"

echo "[start] $(date)" | tee "${LOG}"

echo "============================================================" | tee -a "${LOG}"
echo "[download] Sakaue2021_BBJ_Myocarditis / GCST90018662" | tee -a "${LOG}"

wget -c \
  -O "${BASE}/data/raw/outcome/sakaue_2021_BBJ/hum0197.v3.BBJ.Myo.v1.zip" \
  "https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.Myo.v1.zip" \
  2>&1 | tee -a "${LOG}"

echo "============================================================" | tee -a "${LOG}"
echo "[download] Sakaue2021_EUR_Myocarditis / GCST90018882" | tee -a "${LOG}"

wget -c \
  -O "${BASE}/data/raw/outcome/sakaue_2021_EUR/hum0197.v3.EUR.Myo.v1.zip" \
  "https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.EUR.Myo.v1.zip" \
  2>&1 | tee -a "${LOG}"

echo "============================================================" | tee -a "${LOG}"
echo "[download] FinnGen_R12_I9_MYOCARD" | tee -a "${LOG}"

wget -c \
  -O "${BASE}/data/raw/outcome/finngen_R12/finngen_R12_I9_MYOCARD.gz" \
  "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_MYOCARD.gz" \
  2>&1 | tee -a "${LOG}"

wget -c \
  -O "${BASE}/data/raw/outcome/finngen_R12/finngen_R12_I9_MYOCARD.gz.tbi" \
  "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_MYOCARD.gz.tbi" \
  2>&1 | tee -a "${LOG}"

echo "============================================================" | tee -a "${LOG}"
echo "[check] downloaded files" | tee -a "${LOG}"

ls -lh \
  "${BASE}/data/raw/outcome/sakaue_2021_BBJ/" \
  "${BASE}/data/raw/outcome/sakaue_2021_EUR/" \
  "${BASE}/data/raw/outcome/finngen_R12/" \
  2>&1 | tee -a "${LOG}"

echo "[done] $(date)" | tee -a "${LOG}"
