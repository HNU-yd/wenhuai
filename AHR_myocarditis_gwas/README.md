
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