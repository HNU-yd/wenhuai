#!/usr/bin/env bash
set -euo pipefail

ROOT="${WENHUAI_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$ROOT"

mkdir -p logs data/metadata data/h5ad_per_sample data/h5ad_merged results/v0

CORE_DATASETS="GSE183716 CNP0005824 GSE166489 GSE167029"

echo "[1/6] standardize prefixed 10x by copy"
python scripts/02_standardize_prefixed_10x_copy.py \
  --datasets CNP0005824 GSE166489 GSE167029 \
  2>&1 | tee logs/02_standardize_prefixed_10x_copy.log

echo "[2/6] discover local data"
python scripts/01_discover_local_data.py \
  --datasets CNP0005824 GSE166489 GSE167029 GSE183716 \
  2>&1 | tee logs/01_discover_local_data.log

echo "[3/6] build GEO sample table"
python scripts/03_build_geo_sample_table.py \
  2>&1 | tee logs/03_build_geo_sample_table.log

echo "[4/6] prepare h5ad per sample"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

python scripts/04_prepare_h5ad_from_inventory.py \
  --datasets ${CORE_DATASETS} \
  --workers 12 \
  2>&1 | tee logs/04_prepare_core_h5ad.log

echo "[5/6] concat h5ad by dataset"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8

python scripts/05_concat_h5ad_on_disk.py \
  --datasets ${CORE_DATASETS} \
  --overwrite \
  2>&1 | tee logs/05_concat_core_h5ad.log

echo "[6/6] run V0 analysis"

for ds in GSE183716 GSE166489 CNP0005824 GSE167029; do
  f="${ROOT}/data/h5ad_merged/${ds}.raw_merged.h5ad"
  if [ ! -s "$f" ]; then
    echo "[skip] missing merged h5ad: $f"
    continue
  fi

  if [ "$ds" = "GSE167029" ]; then
    threads=12
    hvg=3000
  elif [ "$ds" = "CNP0005824" ]; then
    threads=12
    hvg=3000
  else
    threads=8
    hvg=3000
  fi

  echo "[V0] $ds threads=$threads hvg=$hvg"
  python scripts/06_run_v0_scanpy.py \
    --input "$f" \
    --dataset "$ds" \
    --threads "$threads" \
    --min_genes 200 \
    --max_mito 20 \
    --n_hvg "$hvg" \
    2>&1 | tee "logs/06_run_v0_${ds}.log"
done

python scripts/07_check_v0_outputs.py \
  --datasets GSE183716 GSE166489 CNP0005824 GSE167029 \
  2>&1 | tee logs/07_check_v0_outputs.log

echo "[done] V0 core pipeline finished"
