# 数据下载清单与当前状态

检查日期：2026-05-18

## Conda 环境

已创建本地 conda 环境：

```bash
conda activate wenhuai
```

环境定义文件位于项目根目录：

```bash
environment.wenhuai.yml
```

已验证可用的核心工具：

- Python: `pandas`, `numpy`, `matplotlib`, `scipy`
- R: `data.table`, `coloc`, `ggplot2`
- 命令行工具：`plink`, `tabix`, `aria2c`, `wget`, `curl`, `unzip`

`TwoSampleMR` 未能在当前机器联网安装：GitHub、MRCIEU r-universe 和 CRAN 访问均超时；conda-forge/bioconda 也没有 `r-twosamplemr` 包。网络可用后可执行：

```bash
conda run -n wenhuai Rscript -e "install.packages('TwoSampleMR', repos=c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))"
```

或：

```bash
conda run -n wenhuai Rscript -e "remotes::install_github('MRCIEU/TwoSampleMR')"
```

## 本地已存在的主要数据

### 心肌炎 outcome GWAS

| 数据 | 本地路径 | 下载网址 | 当前状态 |
| --- | --- | --- | --- |
| Sakaue2021 BBJ myocarditis | `AHR_myocarditis_gwas/data/raw/outcome/sakaue_2021_BBJ/hum0197.v3.BBJ.Myo.v1.zip` | `https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.BBJ.Myo.v1.zip` | 本地已存在；本机 HEAD 检查超时 |
| Sakaue2021 EUR myocarditis | `AHR_myocarditis_gwas/data/raw/outcome/sakaue_2021_EUR/hum0197.v3.EUR.Myo.v1.zip` | `https://humandbs.dbcls.jp/files/hum0197/hum0197.v3.EUR.Myo.v1.zip` | 本地已存在；本机 HEAD 检查 200 OK |
| FinnGen R12 I9_MYOCARD | `AHR_myocarditis_gwas/data/raw/outcome/finngen_R12/finngen_R12_I9_MYOCARD.gz` | `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_MYOCARD.gz` | 本地已存在；本机 HEAD 检查 200 OK |
| FinnGen R12 I9_MYOCARD index | `AHR_myocarditis_gwas/data/raw/outcome/finngen_R12/finngen_R12_I9_MYOCARD.gz.tbi` | `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_MYOCARD.gz.tbi` | 本地已存在；本机 HEAD 检查 200 OK |

可重下脚本：

```bash
conda run -n wenhuai bash AHR_myocarditis_gwas/scripts/01_download_outcome_gwas.sh
```

### AHR eQTLGen exposure

| 数据 | 本地路径 | 下载网址 | 当前状态 |
| --- | --- | --- | --- |
| eQTLGen FDR0.05 significant cis-eQTL | `AHR_myocarditis_gwas/data/raw/exposure/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz` | `https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz` | 本地已存在；本机 DNS 解析失败 |
| eQTLGen full cis-eQTL | `AHR_myocarditis_gwas/data/raw/exposure/eqtlgen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz` | `https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz` | 本地已存在；本机 DNS 解析失败 |
| eQTLGen allele frequency | `AHR_myocarditis_gwas/data/raw/exposure/eqtlgen/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz` | `https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz` | 本地已存在；本机 DNS 解析失败 |

可重下脚本：

```bash
conda run -n wenhuai bash AHR_myocarditis_gwas/scripts/05_download_eqtlgen_AHR_exposure.sh
```

如果 `molgenis26.gcc.rug.nl` 在当前网络不可解析，可尝试镜像：

- `https://tf.lisanwanglab.org/GADB/FILER2/Annotationtracks/Downloads/eQTL_gen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz`
- `https://tf.lisanwanglab.org/GADB/FILER2/Annotationtracks/Downloads/eQTL_gen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz`
- `https://tf.lisanwanglab.org/GADB/FILER2/Annotationtracks/Downloads/eQTL_gen/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz`

### eQTLGen trans-eQTL exposure（历史下载记录，非 add-project 主交付）

用途：此前按数据下载请求检查 Kynurenine 工具 SNP 在 AHR trans-expression summary statistics 中的覆盖情况。该排查不属于 `加项目.docx` 当前主交付；当前 add-project 主交付只使用文档指定的 `rs17643734`、`rs59291726`。

| 数据 | 本地路径 | 下载网址 | 当前状态 |
| --- | --- | --- | --- |
| eQTLGen trans-eQTL summary | `AHR_myocarditis_gwas/data/raw/exposure/eqtlgen/trans/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz` | `https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz` | 已下载；`gzip -t` 通过；大小约 5.7G |

用户提供的镜像 URL：

```text
https://tf.lisanwanglab.org/GADB/FILER2/Annotationtracks/Downloads/eQTL_gen/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz
```

该镜像在本机 HEAD 检查可返回 200，但实际传输时 `wget`/`curl` 中断。已改用 eQTLGen/Molgenis 官方下载地址；该站证书当前过期，因此下载时需要关闭证书校验。

可复用下载命令：

```bash
mkdir -p AHR_myocarditis_gwas/data/raw/exposure/eqtlgen/trans
conda run -n wenhuai aria2c \
  --check-certificate=false \
  --continue=true \
  --max-connection-per-server=8 \
  --split=8 \
  -d AHR_myocarditis_gwas/data/raw/exposure/eqtlgen/trans \
  -o 2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz \
  https://molgenis26.gcc.rug.nl/downloads/eqtlgen/trans-eqtl/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz
gzip -t AHR_myocarditis_gwas/data/raw/exposure/eqtlgen/trans/2018-09-04-trans-eQTLsFDR-CohortInfoRemoved-BonferroniAdded.txt.gz
```

历史检查结果：

| SNP | 是否在 eQTLGen trans 文件中出现 | 是否有 AHR/ENSG00000106546 记录 |
| --- | --- | --- |
| `rs4843270` | 否 | 否 |
| `rs61825638` | 否 | 否 |
| `rs3184504` | 是 | 是 |
| `rs6540080` | 否 | 否 |
| `rs10216901` | 否 | 否 |

筛选结果保留在未跟踪的原始数据目录：

```bash
AHR_myocarditis_gwas/data/raw/exposure/eqtlgen/trans/kynurenine_5snp_AHR_trans_eqtlgen_hits.tsv
```

命中内容为 `rs3184504 -> AHR` 一条记录：`Pvalue=0.03585218`, `Zscore=2.0987`, `NrSamples=28204`, `FDR=0.9126160448390098`, `BonferroniP=1`。

### 1000 Genomes LD reference

| 数据 | 本地路径 | 下载网址 | 当前状态 |
| --- | --- | --- | --- |
| 1000G LD reference archive | `AHR_myocarditis_gwas/data/ref/1000G/1kg.v3.tgz` | `http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz` | 本地已存在 |
| EUR/EAS/AFR/AMR/SAS PLINK files | `AHR_myocarditis_gwas/data/ref/1000G/{POP}.bed/.bim/.fam` | 解压自 `1kg.v3.tgz` | 本地已存在 |

如需重下：

```bash
mkdir -p AHR_myocarditis_gwas/data/ref/1000G
wget -c http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz -O AHR_myocarditis_gwas/data/ref/1000G/1kg.v3.tgz
tar -xzf AHR_myocarditis_gwas/data/ref/1000G/1kg.v3.tgz -C AHR_myocarditis_gwas/data/ref/1000G
```

## 当前仍缺或不能直接联网取得的内容

### OpenGWAS AHR expression 查询（历史排查，非 add-project 主交付）

用途：此前用于排查 Kynurenine 工具 SNP 在 AHR expression GWAS 中的覆盖情况。当前 add-project 主交付只按 `加项目.docx` 指定的 `rs17643734`、`rs59291726` 处理。

需要的数据：

- OpenGWAS dataset ID: `eqtl-a-ENSG00000106546`
- 至少需要 Kynurenine 工具 SNP 在该 AHR expression GWAS 中的 `beta`, `se`, `effect_allele`, `other_allele`, `eaf`, `pval`。
- 建议本地放置路径：`AHR_myocarditis_gwas/data/raw/exposure/opengwas/eqtl-a-ENSG00000106546.associations.tsv.gz`

相关网址：

- 数据集检索页：`https://opengwas.io/datasets/`
- API 文档：`https://api.opengwas.io/api/`
- AHR dataset 查询：`https://api.opengwas.io/api/gwasinfo?id=eqtl-a-ENSG00000106546`

限制：

- OpenGWAS 绝大多数 API endpoint 当前要求 JWT。
- 需要在 `https://api.opengwas.io/profile/` 登录后生成 token。
- 使用 R 包时可在环境变量中设置：`OPENGWAS_JWT=<your_token>`。

2026-05-18 使用 JWT 检查后的结果：

- `/api/associations` direct 查询已可用，并已保存到本地：
  `AHR_myocarditis_gwas/data/raw/exposure/opengwas/eqtl-a-ENSG00000106546.associations.tsv.gz`
- direct 查询 5 个 Kynurenine 工具 SNP 只返回 `rs3184504`。
- `proxies=1` 查询只额外返回 `rs4843270` 的 proxy 信息，仍缺 `rs61825638`、`rs6540080`、`rs10216901`。
- 已通过 `/api/gwasinfo/files` 下载完整 OpenGWAS VCF 到本地：
  `AHR_myocarditis_gwas/data/raw/exposure/opengwas/eqtl-a-ENSG00000106546.vcf.gz`
  和 `.tbi`。
- 完整 VCF 的 header 显示 `TotalVariants=19830`，基于 GRCh37/HG19。
- 完整 VCF 中 5 个 Kynurenine 工具 SNP 仍只有 `rs3184504` 存在，因此缺失不是 associations 接口漏查，而是该 OpenGWAS AHR eQTL dataset 本身不包含另外 4 个 direct SNP。

可复用脚本：

```bash
# 拉取指定 SNP association；JWT 可通过 OPENGWAS_JWT 或 stdin 传入
conda run -n wenhuai python AHR_myocarditis_gwas/scripts/16_fetch_opengwas_ahr_associations.py

# 下载 OpenGWAS 返回的完整 dataset VCF 和 tbi；JWT 可通过 OPENGWAS_JWT 或 stdin 传入
conda run -n wenhuai python AHR_myocarditis_gwas/scripts/17_fetch_opengwas_dataset_files.py
```

结论：不再需要额外下载该 OpenGWAS dataset 的“完整 VCF”，因为完整 VCF 已经下载并检查过；该历史排查不纳入当前 add-project 主报告。

2026-05-18 继续下载并检查 eQTLGen trans-eQTL summary 后，结果同样只有 `rs3184504` 命中 AHR/ENSG00000106546。该信息仅作为数据覆盖记录保留，不作为 `加项目.docx` 主分析结果。
