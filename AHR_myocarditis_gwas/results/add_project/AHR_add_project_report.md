# 加项目.docx 补充分析报告

## 完成内容

根据 `加项目.docx`，本次围绕 AHR 表达 GWAS `eqtl-a-ENSG00000106546` 和两个指定显著 SNP（`rs17643734`、`rs59291726`）补充了中介 MR 相关材料、共定位结果汇总和 Manhattan/locus 图。

已完成的本地可计算部分：

- 使用本地 eQTLGen AHR cis-eQTL 表提取两个 SNP 的 AHR 表达效应。
- 在三套心肌炎 GWAS 的 AHR locus 数据中匹配两个 SNP，并计算 AHR expression -> myocarditis 的单 SNP Wald ratio。
- 生成中介模式图、AHR eQTL locus Manhattan 图、心肌炎 AHR locus Manhattan 图、coloc 汇总图和两 SNP forest 图。

说明：本报告只纳入 `加项目.docx` 指定的 `rs17643734` 和 `rs59291726`。其他未在文档中指定的 SNP 排查不属于本次主交付。

## 两 SNP Wald ratio 结果

| outcome_dataset | target_snp | status | match_method | beta_outcome_aligned | se_outcome | pval_outcome | ratio | ratio_se | ratio_pval | ratio_or | ratio_or_lci95 | ratio_or_uci95 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Sakaue2021_BBJ_Myocarditis | rs17643734 | missing_outcome | missing | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| Sakaue2021_BBJ_Myocarditis | rs59291726 | ok | rsID | 0.492 | 0.558 | 0.378 | 1.800 | 2.045 | 0.379 | 6.052 | 0.110 | 333.218 |
| Sakaue2021_EUR_Myocarditis | rs17643734 | ok | chr_pos_alleles | 0.131 | 0.115 | 0.257 | 0.274 | 0.242 | 0.257 | 1.315 | 0.819 | 2.113 |
| Sakaue2021_EUR_Myocarditis | rs59291726 | ok | chr_pos_alleles | 0.183 | 0.107 | 0.088 | 0.669 | 0.394 | 0.089 | 1.952 | 0.902 | 4.222 |
| FinnGen_R12_I9_MYOCARD | rs17643734 | ok | rsID | -0.004 | 0.063 | 0.945 | -0.009 | 0.131 | 0.945 | 0.991 | 0.766 | 1.282 |
| FinnGen_R12_I9_MYOCARD | rs59291726 | ok | rsID | -0.002 | 0.061 | 0.978 | -0.006 | 0.224 | 0.978 | 0.994 | 0.641 | 1.541 |

## 中介 MR 状态

| mediation_arm | status | method | nsnp | beta | se | pval | or | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| AHR expression -> Myocarditis | completed_two_snp_wald_where_available | Single-SNP Wald ratio | 2 | NA | NA | NA | NA | Docx-specified SNPs rs17643734 and rs59291726 were evaluated against the three myocarditis GWAS datasets; see AHR_two_significant_snp_wald_ratios.tsv. |
| AHR eQTL and Myocarditis coloc | completed | coloc.abf | 3 | NA | NA | NA | NA | Coloc summary and SNP-level posterior context were generated for the requested AHR SNPs. |

## Coloc 汇总

| outcome_dataset | n_overlap | PP.H3.abf | PP.H4.abf | coloc_interpretation |
| --- | --- | --- | --- | --- |
| Sakaue2021_BBJ_Myocarditis | 5192 | 0.309 | 0.062 | no_strong_colocalization |
| Sakaue2021_EUR_Myocarditis | 7383 | 0.356 | 0.048 | no_strong_colocalization |
| FinnGen_R12_I9_MYOCARD | 7171 | 0.284 | 0.021 | no_strong_colocalization |

三套心肌炎数据的 PP.H4 均低于常用强共定位参考阈值，既有结果解释为 `no_strong_colocalization`。

## 指定 SNP 的 coloc posterior

| outcome_dataset | snp | match_method | beta_gwas_aligned | se_gwas | pval_gwas | SNP.PP.H4 | outcome_variant_id |
| --- | --- | --- | --- | --- | --- | --- | --- |
| FinnGen_R12_I9_MYOCARD | rs17643734 | rsID | -0.004 | 0.063 | 0.945 | 1.000 | 7:17123258:A:G |
| FinnGen_R12_I9_MYOCARD | rs59291726 | rsID | -0.002 | 0.061 | 0.978 | 1.917e-94 | 7:17171454:C:T |
| Sakaue2021_BBJ_Myocarditis | rs59291726 | rsID | 0.492 | 0.558 | 0.378 | 1 | 7:17211078:C:T |
| Sakaue2021_EUR_Myocarditis | rs17643734 | chr_pos | 0.131 | 0.115 | 0.257 | 1.000 | 7:17162882:A:G |
| Sakaue2021_EUR_Myocarditis | rs59291726 | chr_pos | 0.183 | 0.107 | 0.088 | 3.513e-94 | 7:17211078:C:T |

注：`SNP.PP.H4` 来自既有 SNP-level posterior 文件，用于显示在 H4 假设下指定 SNP 的相对定位信息；整体共定位证据仍以主结果表的 `PP.H4.abf` 为准。

## 输出文件

- `results/add_project/AHR_two_significant_snp_wald_ratios.tsv`
- `results/add_project/AHR_two_snp_mediation_status.tsv`
- `results/add_project/AHR_two_snp_coloc_posterior_context.tsv`
- `results/figures/AHR_two_snp_mediation_model.png`
- `results/figures/AHR_eqtlgen_AHR_locus_manhattan_two_snp.png`
- `results/figures/AHR_myocarditis_locus_manhattan_two_snp.png`
- `results/figures/AHR_coloc_summary_with_two_snp_context.png`
- `results/figures/AHR_two_snp_wald_ratio_forest.png`

## 复现命令

```bash
conda run -n wenhuai python AHR_myocarditis_gwas/scripts/15_add_docx_ahr_mediation_figures.py
```
