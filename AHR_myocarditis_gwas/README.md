
# AHR 与心肌炎 GWAS 分析项目

## 1. 项目目的

本项目用于分析 AHR 遗传调控信号与心肌炎 GWAS 风险之间的关系。分析对象包括三套并列主分析心肌炎 GWAS 数据：

1. `Sakaue2021_BBJ_Myocarditis`
   - 对应 GWAS Catalog: `GCST90018662`
   - 对应 OpenGWAS: `ebi-a-GCST90018662`
   - 人群：East Asian / BBJ
   - 疾病：Myocarditis

2. `Sakaue2021_EUR_Myocarditis`
   - 对应 GWAS Catalog: `GCST90018882`
   - 对应 OpenGWAS: `ebi-a-GCST90018882`
   - 人群：European
   - 疾病：Myocarditis

3. `FinnGen_R12_I9_MYOCARD`
   - 对应 FinnGen endpoint: `I9_MYOCARD`
   - 对应 OpenGWAS: `finn-b-I9_MYOCARD`
   - 人群：Finnish / European
   - 疾病：Myocarditis

三套数据均作为并列主分析数据库，不区分主分析和补充分析。

---

## 2. 目录结构

```text
AHR_myocarditis_gwas/
├── config/
│   └── outcome_gwas_datasets.tsv
├── scripts/
│   ├── 01_download_outcome_gwas.sh
│   └── 02_inspect_raw_outcome_gwas.sh
├── docs/
│   ├── analysis_log.md
│   └── file_manifest.tsv
├── logs/
├── data/
│   ├── raw/
│   │   └── outcome/
│   │       ├── sakaue_2021_BBJ/
│   │       ├── sakaue_2021_EUR/
│   │       └── finngen_R12/
│   └── formatted/
│       └── outcome/
└── results/
    ├── qc/
    ├── mr/
    ├── coloc/
    ├── smr/
    └── figures/
````

---

## 3. `加项目.docx` 补充分析

补充任务围绕 AHR expression GWAS `eqtl-a-ENSG00000106546`，只分析两个指定显著 SNP：

- `rs17643734`
- `rs59291726`

可复现脚本：

```bash
conda run -n wenhuai python AHR_myocarditis_gwas/scripts/15_add_docx_ahr_mediation_figures.py
```

主要输出：

- `results/add_project/AHR_two_significant_snp_wald_ratios.tsv`
- `results/add_project/AHR_two_snp_mediation_status.tsv`
- `results/add_project/AHR_two_snp_coloc_posterior_context.tsv`
- `results/add_project/AHR_add_project_report.md`
- `results/add_project/AHR_add_project_word_report.docx`
- `results/figures/AHR_two_snp_mediation_model.png`
- `results/figures/AHR_eqtlgen_AHR_locus_manhattan_two_snp.png`
- `results/figures/AHR_myocarditis_locus_manhattan_two_snp.png`
- `results/figures/AHR_coloc_summary_with_two_snp_context.png`
- `results/figures/AHR_two_snp_wald_ratio_forest.png`

说明：本地已按 `加项目.docx` 完成 AHR expression -> myocarditis 的两 SNP Wald ratio、coloc 汇总、Manhattan/locus 图和中介模式图。本节主交付只纳入文档指定的 `rs17643734`、`rs59291726`。

Word 报告可复现生成：

```bash
conda run -n wenhuai python AHR_myocarditis_gwas/scripts/18_build_add_project_word_report.py
```
