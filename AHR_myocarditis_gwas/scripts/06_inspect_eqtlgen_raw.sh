#!/usr/bin/env bash
set -euo pipefail

BASE="${AHR_GWAS_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
RAW_DIR="${BASE}/data/raw/exposure/eqtlgen"
OUT="${BASE}/results/qc/eqtlgen_raw_header_check.txt"
MANIFEST="${BASE}/results/qc/eqtlgen_file_manifest.tsv"
LOG="${BASE}/logs/06_inspect_eqtlgen_raw.log"

mkdir -p "${BASE}/results/qc" "${BASE}/logs"

echo "[start] $(date)" | tee "${LOG}"

printf "file_key\tpath\tsize_bytes\tsha256\tfile_type\n" > "${MANIFEST}"
: > "${OUT}"

inspect_one () {
    local file_key="$1"
    local f="$2"

    {
        echo "============================================================"
        echo "FILE_KEY: ${file_key}"
        echo "FILE: ${f}"
    } >> "${OUT}"

    if [[ ! -s "${f}" ]]; then
        echo "[missing] ${f}" | tee -a "${LOG}"
        echo "[missing] file does not exist or is empty" >> "${OUT}"
        echo "" >> "${OUT}"
        return 0
    fi

    local size
    local sha
    local ftype

    size=$(stat -c%s "${f}")
    sha=$(sha256sum "${f}" | awk '{print $1}')
    ftype=$(file -b "${f}")

    printf "%s\t%s\t%s\t%s\t%s\n" "${file_key}" "${f}" "${size}" "${sha}" "${ftype}" >> "${MANIFEST}"

    {
        echo "SIZE_BYTES: ${size}"
        echo "SHA256: ${sha}"
        echo "FILE_TYPE: ${ftype}"
        echo ""
        echo "[HEAD]"
    } >> "${OUT}"

    zcat "${f}" | head -n 8 >> "${OUT}" || true

    {
        echo ""
        echo "[HEADER_ONLY]"
    } >> "${OUT}"

    zcat "${f}" | head -n 1 | tr '\t' '\n' | nl -ba >> "${OUT}" || true

    echo "" >> "${OUT}"
}

inspect_one \
  "eqtlgen_significant_FDR0.05" \
  "${RAW_DIR}/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"

inspect_one \
  "eqtlgen_full_cis_eQTL" \
  "${RAW_DIR}/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"

inspect_one \
  "eqtlgen_allele_frequency" \
  "${RAW_DIR}/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"

echo "[done] header check -> ${OUT}" | tee -a "${LOG}"
echo "[done] manifest -> ${MANIFEST}" | tee -a "${LOG}"
echo "[done] $(date)" | tee -a "${LOG}"
