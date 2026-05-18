##1.复现流程
1. 下载 1400 个代谢物 GWAS summary statistics
2. 下载 FinnGen myocarditis GWAS
3. 下载 IEU/OpenGWAS myocarditis GWAS
4. 对每个代谢物筛选工具变量：
   p < 5e-8
   LD clumping: kb = 10000, r2 = 0.001
   F-statistic > 10
5. 与 myocarditis outcome harmonise
6. Steiger test 去除反向因果 SNP
7. MR 主分析：
   IVW
   MR-Egger
   Weighted median
   Weighted mode
8. FDR 校正
9. 筛选：
   IVW FDR p < 0.05
   四种 MR 方法 beta 方向一致
   无明显异质性
   无水平多效性
   MR-PRESSO 无 outlier 或 outlier 修正后仍显著
10. 对 5 个候选代谢物做：
   leave-one-out
   scatter plot
   funnel plot
   forest plot
   coloc
   IEU 验证集 meta-analysis

## 2. 重写 README.md

````bash
cd /home/data1/wenhuai/myocarditis_mr_reproduce

cat > README.md <<'MD'
# Myocarditis Blood Metabolome MR Reproduction

本项目用于复现论文：

**Association between human blood metabolome and risk of myocarditis: a Mendelian randomization study**

论文目标是使用双样本孟德尔随机化分析人类血液代谢物与心肌炎风险之间的因果关系。暴露数据来自 Chen et al. 的 1400 个 plasma metabolome GWAS，发现集结局为 FinnGen R10 myocarditis，后续可扩展到 IEU OpenGWAS 验证集、meta-analysis 和 coloc 分析。

当前项目已经完成论文 Table 1 中 5 个最终代谢物的 FinnGen MR 主流程复现。

---

## 1. Project path

```bash
/home/data1/wenhuai/myocarditis_mr_reproduce
````

目录结构：

```text
myocarditis_mr_reproduce/
├── exposure/          # GWAS Catalog 代谢物暴露数据
├── outcome/           # FinnGen myocarditis 结局数据
├── ld/                # 1000 Genomes EUR LD reference
├── scripts/           # 最终复现脚本
├── results/           # 输出表格
├── plots/             # 输出图
├── logs/              # 运行日志
└── README.md
```

---

## 2. Final scripts

当前保留的核心脚本如下：

```text
scripts/
├── 00_install_R_packages.R
├── download_final5_exposures_v2.sh
├── 01_rebuild_final5_exposure_instruments.R
├── 02_rerun_outcome_harmonise_mr.R
└── 02_get_ieu_download_link.R
```

其中：

```text
00_install_R_packages.R
```

用于安装 R 包。

```text
download_final5_exposures_v2.sh
```

用于下载 5 个最终代谢物的 GWAS Catalog 暴露数据。

```text
01_rebuild_final5_exposure_instruments.R
```

用于从原始暴露 GWAS 中筛选工具变量，执行 PLINK LD clumping，并生成：

```text
results/final5_exposure_instruments.tsv
results/final5_exposure_nsnp_counts.tsv
```

```text
02_rerun_outcome_harmonise_mr.R
```

用于读取 FinnGen outcome，执行 outcome matching、harmonise、MR 和 sensitivity analysis。

```text
02_get_ieu_download_link.R
```

预留给后续 IEU OpenGWAS 验证集下载。

---

## 3. Data

### 3.1 Exposure data

当前复现论文 Table 1 中的 5 个代谢物：

| Metabolite                                        | GWAS Catalog accession | Expected nSNP |
| ------------------------------------------------- | ---------------------: | ------------: |
| Kynurenine levels                                 |           GCST90199636 |             5 |
| 1-stearoyl-GPE (18:0) levels                      |           GCST90199772 |             7 |
| Deoxycarnitine levels                             |           GCST90199813 |             5 |
| X-25422 levels                                    |           GCST90200661 |             6 |
| 5-acetylamino-6-formylamino-3-methyluracil levels |           GCST90200680 |             4 |

暴露数据路径：

```bash
exposure/GCST*/GCST*_buildGRCh38.tsv.gz
```

下载命令：

```bash
bash scripts/download_final5_exposures_v2.sh 2>&1 | tee logs/download_final5_exposures_v2.log
```

检查暴露数据：

```bash
for f in exposure/GCST*/GCST*_buildGRCh38.tsv.gz; do
  echo "========== $f =========="
  ls -lh "$f"
  gzip -t "$f" && echo "gzip_ok"
  zcat "$f" | head -1
done
```

---

### 3.2 Outcome data

FinnGen R10 myocarditis：

```text
finngen_R10_I9_MYOCARD.gz
```

下载路径：

```bash
outcome/finngen_R10_I9_MYOCARD.gz
```

下载命令：

```bash
cd /home/data1/wenhuai/myocarditis_mr_reproduce/outcome

wget -c https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_MYOCARD.gz
wget -c https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_I9_MYOCARD.gz.tbi
```

---

### 3.3 LD reference

使用 1000 Genomes European LD reference：

```text
ld/EUR.bed
ld/EUR.bim
ld/EUR.fam
```

下载命令：

```bash
cd /home/data1/wenhuai/myocarditis_mr_reproduce/ld

wget -c http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz
tar -xzf 1kg.v3.tgz
```

---

## 4. Environment

推荐使用 conda 环境：

```bash
conda create -n mr43 -y -c conda-forge \
  r-base=4.3.1 \
  r-data.table r-dplyr r-stringr r-tibble r-readr r-jsonlite \
  r-ggplot2 r-remotes r-devtools r-optparse r-meta r-metafor

conda activate mr43
```

安装 R 包：

```bash
Rscript scripts/00_install_R_packages.R 2>&1 | tee logs/install_R_packages.log
```

---

## 5. Reproduction pipeline

### Step 1: rebuild exposure instruments

```bash
conda activate mr43
cd /home/data1/wenhuai/myocarditis_mr_reproduce

Rscript scripts/01_rebuild_final5_exposure_instruments.R \
  2>&1 | tee logs/01_rebuild_final5_exposure_instruments.log
```

检查工具变量数量：

```bash
cat results/final5_exposure_nsnp_counts.tsv | column -t -s $'\t'
```

期望结果：

```text
GCST90199636    5
GCST90199772    7
GCST90199813    5
GCST90200661    6
GCST90200680    4
```

重要说明：

本项目必须先使用原始 GWAS p 值执行 PLINK clumping，然后再进入 `TwoSampleMR::format_data()`。

错误流程：

```text
format_data -> pval.exposure -> PLINK clump
```

正确流程：

```text
raw GWAS p value -> PLINK clump -> F >= 10 -> format_data -> outcome -> harmonise -> MR
```

原因是 `TwoSampleMR::format_data()` 可能会截断极小 p 值，从而影响 PLINK clump 的排序和 index SNP 选择。对于 `GCST90200680`，如果先 `format_data()` 再 clump，会得到 3 个 SNP；使用原始 p 值 clump，则得到论文一致的 4 个 SNP。

---

### Step 2: outcome matching, harmonise, MR

```bash
Rscript scripts/02_rerun_outcome_harmonise_mr.R \
  2>&1 | tee logs/02_rerun_outcome_harmonise_mr.log
```

检查每一步 SNP 数：

```bash
cat results/clean_outcome_mr/08_nsnp_counts_all_steps.tsv | column -t -s $'\t'
```

期望结果：

```text
Kynurenine levels                                      5
1-stearoyl-GPE (18:0) levels                           7
Deoxycarnitine levels                                  5
X-25422 levels                                         6
5-acetylamino-6-formylamino-3-methyluracil levels      4
```

查看 IVW 主结果：

```bash
cat results/clean_outcome_mr/11_final5_ivw_only.tsv | column -t -s $'\t' | less -S
```

查看全部 MR 结果：

```bash
cat results/clean_outcome_mr/10_final5_mr_results.tsv | column -t -s $'\t' | less -S
```

---

## 6. Main output files

### Exposure instruments

```text
results/final5_exposure_instruments.tsv
results/final5_exposure_nsnp_counts.tsv
```

### Outcome matching and harmonise

```text
results/clean_outcome_mr/02_outcome_match_status.tsv
results/clean_outcome_mr/03_finngen_outcome_hits.tsv
results/clean_outcome_mr/06_harmonised_all.tsv
results/clean_outcome_mr/07_harmonised_mr_keep.tsv
results/clean_outcome_mr/08_nsnp_counts_all_steps.tsv
```

### MR results

```text
results/clean_outcome_mr/10_final5_mr_results.tsv
results/clean_outcome_mr/11_final5_ivw_only.tsv
```

### Sensitivity results

```text
results/clean_outcome_mr/12_heterogeneity.tsv
results/clean_outcome_mr/13_pleiotropy_egger_intercept.tsv
results/clean_outcome_mr/14_single_snp.tsv
results/clean_outcome_mr/15_leave_one_out.tsv
```

---

## 7. Key MR parameters

论文主流程对应参数：

```text
Exposure SNP threshold: p < 5e-8
LD clumping window: 10000 kb
LD r2 threshold: 0.001
Weak IV filtering: F-statistic >= 10
Primary MR method: IVW
Supplementary MR methods:
  - MR-Egger
  - weighted median
  - weighted mode
Multiple testing correction: FDR
Sensitivity analysis:
  - Cochran Q test
  - MR-Egger intercept
  - single SNP analysis
  - leave-one-out
```

---

## 8. CPU / GPU

该项目是 MR 统计分析，不是深度学习训练。

主要计算包括：

```text
PLINK LD clumping
GWAS summary statistics parsing
FinnGen outcome matching
harmonise
IVW / MR-Egger / weighted median / weighted mode
sensitivity analysis
```

这些任务主要使用 CPU 和磁盘 IO，基本不能有效使用 GPU。

如果后续扩展到 1400 个代谢物全量复现，优化方向应是 CPU 多进程并行，而不是 GPU 加速。

---

## 9. Next steps

当前已完成：

```text
1. 5 个最终代谢物的暴露工具变量重建
2. FinnGen outcome matching
3. harmonise
4. IVW / MR-Egger / weighted median / weighted mode
5. heterogeneity / pleiotropy / single SNP / leave-one-out
```

下一步可继续做：

```text
1. 生成论文风格 scatter / funnel / forest / leave-one-out 图
2. 接入 IEU OpenGWAS myocarditis 验证集
3. FinnGen + IEU meta-analysis
4. coloc 共定位分析
5. 扩展到 1400 个 metabolite / metabolite ratio 全量筛选
```

MD

````

---

## 3. 删除没用的老版本脚本

这个只删 `scripts/` 里的旧代码，不动数据和结果。

```bash
cd /home/data1/wenhuai/myocarditis_mr_reproduce

mkdir -p logs

find scripts -maxdepth 1 -type f \
  ! -name '00_install_R_packages.R' \
  ! -name 'download_final5_exposures_v2.sh' \
  ! -name '01_rebuild_final5_exposure_instruments.R' \
  ! -name '02_rerun_outcome_harmonise_mr.R' \
  ! -name '02_get_ieu_download_link.R' \
  -print | tee logs/deleted_old_scripts.txt | xargs -r rm -f
````

检查保留下来的脚本：

```bash
find scripts -maxdepth 1 -type f | sort
```

应该只剩：

```text
scripts/00_install_R_packages.R
scripts/01_rebuild_final5_exposure_instruments.R
scripts/02_get_ieu_download_link.R
scripts/02_rerun_outcome_harmonise_mr.R
scripts/download_final5_exposures_v2.sh
```

---

## 4. 可选：删除旧日志

如果你也想清理旧日志，只保留最终流程日志：

```bash
cd /home/data1/wenhuai/myocarditis_mr_reproduce

find logs -maxdepth 1 -type f \
  ! -name 'install_R_packages.log' \
  ! -name 'download_final5_exposures_v2.log' \
  ! -name '01_rebuild_final5_exposure_instruments.log' \
  ! -name '02_rerun_outcome_harmonise_mr.log' \
  ! -name 'deleted_old_scripts.txt' \
  -print | tee logs/deleted_old_logs.txt | xargs -r rm -f
```

---

## 5. 最后重新跑一遍标准流程

```bash
cd /home/data1/wenhuai/myocarditis_mr_reproduce
conda activate mr43

Rscript scripts/01_rebuild_final5_exposure_instruments.R \
  2>&1 | tee logs/01_rebuild_final5_exposure_instruments.log

Rscript scripts/02_rerun_outcome_harmonise_mr.R \
  2>&1 | tee logs/02_rerun_outcome_harmonise_mr.log

cat results/final5_exposure_nsnp_counts.tsv | column -t -s $'\t'
cat results/clean_outcome_mr/08_nsnp_counts_all_steps.tsv | column -t -s $'\t'
cat results/clean_outcome_mr/11_final5_ivw_only.tsv | column -t -s $'\t' | less -S
```
