# Project Structure and Path Policy

## Repository Scope

This Git repository is for portable project assets:

- analysis scripts in `scripts/`, `src/`, `AHR_myocarditis_gwas/scripts/`, `myocarditis_mr_reproduce/scripts/`, and `kyn_ahr_v0_public/scripts/`
- documentation and run notes
- lightweight summary tables and figures, especially `report_files/`
- lightweight final MR/AHR summary outputs under the two MR subprojects

The following are deliberately kept out of Git because they are raw downloads or large generated intermediates:

- `data/`
- top-level `results/`
- `logs/`
- raw GNN/STT/CNP/GEO downloads
- `*.h5ad`, `*.rds`, `*.gem.gz`, LD reference files, large archives
- GWAS raw data and downloaded exposure/outcome files

## Path Handling

Python scripts in the top-level `scripts/` directory use `scripts/project_paths.py` to infer the repository root from the script location. You can override this with:

```bash
export WENHUAI_ROOT=/home/data1/wenhuai
```

The AHR GWAS shell scripts infer their own subproject root. Their R scripts use:

```bash
export AHR_GWAS_ROOT=/home/data1/wenhuai/AHR_myocarditis_gwas
```

The metabolome MR R scripts use:

```bash
export WENHUAI_MR_ROOT=/home/data1/wenhuai/myocarditis_mr_reproduce
```

If those environment variables are not set, the R scripts default to the current working directory, so run them after `cd` into the corresponding subproject.

## GitHub Push Policy

Do not commit raw data or full intermediate objects. The workspace contains files larger than GitHub's normal file limit and hundreds of GB of data. Use external storage or Git LFS only if a future workflow explicitly needs versioned large files.
