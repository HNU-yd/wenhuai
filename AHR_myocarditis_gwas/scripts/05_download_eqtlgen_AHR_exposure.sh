#!/usr/bin/env bash
set -euo pipefail

BASE="${AHR_GWAS_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
RAW_DIR="${BASE}/data/raw/exposure/eqtlgen"
LOG="${BASE}/logs/05_download_eqtlgen_AHR_exposure.log"

mkdir -p "${RAW_DIR}" "${BASE}/logs"

echo "[start] $(date)" | tee "${LOG}"

cd "${RAW_DIR}"

echo "============================================================" | tee -a "${LOG}"
echo "[download] eQTLGen FDR0.05 significant cis-eQTL" | tee -a "${LOG}"

wget -c \
  "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
  -O "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
  2>&1 | tee -a "${LOG}"

echo "============================================================" | tee -a "${LOG}"
echo "[download] eQTLGen full cis-eQTL" | tee -a "${LOG}"

wget -c \
  "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
  -O "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
  2>&1 | tee -a "${LOG}"

echo "============================================================" | tee -a "${LOG}"
echo "[download] eQTLGen allele frequency file" | tee -a "${LOG}"

wget -c \
  "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz" \
  -O "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz" \
  2>&1 | tee -a "${LOG}" || true

echo "============================================================" | tee -a "${LOG}"
echo "[check] downloaded files" | tee -a "${LOG}"

ls -lh "${RAW_DIR}" | tee -a "${LOG}"
file "${RAW_DIR}"/* | tee -a "${LOG}" || true

echo "[done] $(date)" | tee -a "${LOG}"
