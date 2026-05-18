# Kynurenine-AHR V0 public dataset baseline

This package implements the V0 Scanpy baseline for GSE167029, GSE166489, and GSE183716.

## 1. Create environment

```bash
mamba create -n kynahr_v0 python=3.10 -y
mamba activate kynahr_v0
pip install -r requirements.txt
```

## 2. Download processed GEO data

```bash
bash scripts/download_geo_processed.sh /path/to/kyn_ahr_v0_public
```

Optional large Seurat RDS for GSE167029:

```bash
DOWNLOAD_RDS=1 bash scripts/download_geo_processed.sh /path/to/kyn_ahr_v0_public
```

## 3. Build metadata

```bash
python scripts/build_metadata.py
```

The default output is `metadata/metadata_public_v0.tsv` inside the project folder.

## 4. Convert matrices to a combined raw h5ad

```bash
python scripts/convert_public_to_h5ad.py \
  --root /path/to/kyn_ahr_v0_public
```

Output:

```text
/path/to/kyn_ahr_v0_public/processed/public_mis_c_combined_raw.h5ad
/path/to/kyn_ahr_v0_public/processed/loaded_samples.tsv
```

## 5. Run V0 Scanpy baseline

```bash
python scripts/run_v0_scanpy.py \
  --in_h5ad /path/to/kyn_ahr_v0_public/processed/public_mis_c_combined_raw.h5ad \
  --out_h5ad /path/to/kyn_ahr_v0_public/processed/public_mis_c_v0_scanpy.h5ad \
  --figdir /path/to/kyn_ahr_v0_public/figures \
  --tables /path/to/kyn_ahr_v0_public/tables
```

Key outputs:

- `public_mis_c_v0_scanpy.h5ad`
- `public_mis_c_v0_myeloid.h5ad`
- `cell_metadata_with_v0_scores.tsv`
- `celltype_fraction_by_sample.tsv`
- `score_summary_by_group_celltype.tsv`
- UMAP/violin/dotplot figures

## Notes

- V0.0 uses marker-score cell typing and pathway score analysis.
- AHR regulon is currently represented by `AHR_regulon_score_proxy = AHR_response_score` until pySCENIC/DoRothEA is added.
- For manuscript-level analysis, verify GSE183716 clinical labels against the sample sheet/paper before final reporting.
