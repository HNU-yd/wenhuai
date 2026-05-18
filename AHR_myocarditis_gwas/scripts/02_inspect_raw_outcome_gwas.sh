#!/usr/bin/env bash
set -euo pipefail

BASE="${AHR_GWAS_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
OUT="${BASE}/results/qc/raw_outcome_gwas_header_check.txt"
MANIFEST="${BASE}/docs/file_manifest.tsv"
LOG="${BASE}/logs/02_inspect_raw_outcome_gwas.log"

mkdir -p "${BASE}/results/qc" "${BASE}/docs" "${BASE}/logs"

echo "[start] $(date)" | tee "${LOG}"
echo "[out] ${OUT}" | tee -a "${LOG}"
echo "[manifest] ${MANIFEST}" | tee -a "${LOG}"

printf "dataset_key\tpath\tsize_bytes\tsha256\tfile_type\n" > "${MANIFEST}"
: > "${OUT}"

inspect_file () {
    local dataset_key="$1"
    local source_id="$2"
    local population="$3"
    local f="$4"

    {
        echo "============================================================"
        echo "DATASET: ${dataset_key}"
        echo "SOURCE_ID: ${source_id}"
        echo "POPULATION: ${population}"
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

    printf "%s\t%s\t%s\t%s\t%s\n" "${dataset_key}" "${f}" "${size}" "${sha}" "${ftype}" >> "${MANIFEST}"

    {
        echo "SIZE_BYTES: ${size}"
        echo "SHA256: ${sha}"
        echo "FILE_TYPE: ${ftype}"
        echo ""
        echo "[HEAD_OR_ZIP_CONTENT]"
    } >> "${OUT}"

    if [[ "${f}" == *.zip ]]; then
        {
            echo "[zip file list]"
            unzip -l "${f}"
            echo ""
        } >> "${OUT}"

        inner=$(unzip -Z1 "${f}" | head -n 1)
        echo "[first inner file] ${inner}" >> "${OUT}"

        if [[ "${inner}" == *.gz ]]; then
            unzip -p "${f}" "${inner}" | gzip -dc | head -n 12 >> "${OUT}" || true
        else
            unzip -p "${f}" "${inner}" | head -n 12 >> "${OUT}" || true
        fi

    elif [[ "${f}" == *.gz ]]; then
        zcat "${f}" | head -n 12 >> "${OUT}" || true

    else
        head -n 12 "${f}" >> "${OUT}" || true
    fi

    echo "" >> "${OUT}"
}

inspect_file \
  "Sakaue2021_BBJ_Myocarditis" \
  "GCST90018662" \
  "East Asian / BBJ" \
  "${BASE}/data/raw/outcome/sakaue_2021_BBJ/hum0197.v3.BBJ.Myo.v1.zip"

inspect_file \
  "Sakaue2021_EUR_Myocarditis" \
  "GCST90018882" \
  "European" \
  "${BASE}/data/raw/outcome/sakaue_2021_EUR/hum0197.v3.EUR.Myo.v1.zip"

inspect_file \
  "FinnGen_R12_I9_MYOCARD" \
  "FINNGEN_R12_I9_MYOCARD" \
  "Finnish / European" \
  "${BASE}/data/raw/outcome/finngen_R12/finngen_R12_I9_MYOCARD.gz"

echo "[done] $(date)" | tee -a "${LOG}"
echo "[done] header check -> ${OUT}" | tee -a "${LOG}"
echo "[done] manifest -> ${MANIFEST}" | tee -a "${LOG}"
