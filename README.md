# Wenhuai myocarditis KYN/AHR analysis workspace

This repository tracks the portable code, documentation, configs, and lightweight result tables/figures for the Wenhuai myocarditis KYN/AHR project.

The local analysis workspace is currently:

```bash
/home/data1/wenhuai
```

Most Python scripts now infer the project root from their own location. Set `WENHUAI_ROOT` only when running scripts from an unusual copied layout:

```bash
export WENHUAI_ROOT=/home/data1/wenhuai
```

Large raw datasets and intermediate objects are intentionally excluded from Git. See [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) for the path policy and what is tracked.

Main documentation:

- [README_wenhuai_total.md](README_wenhuai_total.md): full MR + scRNA + spatial workflow notes.
- [read_2_v0.md](read_2_v0.md): single-cell V0 analysis notes.
- [read_3_v0.md](read_3_v0.md): STT0000127 spatial transcriptomics notes.
- [AHR_myocarditis_gwas/README.md](AHR_myocarditis_gwas/README.md): AHR eQTL MR/coloc subproject.
- [myocarditis_mr_reproduce/readme.md](myocarditis_mr_reproduce/readme.md): metabolome MR reproduction subproject.
