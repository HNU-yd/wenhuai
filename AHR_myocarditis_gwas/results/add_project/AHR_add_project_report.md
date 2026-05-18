# 加项目.docx 补充分析报告

## 完成内容

根据 `加项目.docx`，本次围绕 AHR 表达 GWAS `eqtl-a-ENSG00000106546` 和两个指定显著 SNP（`rs17643734`、`rs59291726`）补充了中介 MR 相关材料、共定位结果汇总和 Manhattan/locus 图。

已完成的本地可计算部分：

- 使用本地 eQTLGen AHR cis-eQTL 表提取两个 SNP 的 AHR 表达效应。
- 在三套心肌炎 GWAS 的 AHR locus 数据中匹配两个 SNP，并计算 AHR expression -> myocarditis 的单 SNP Wald ratio。
- 使用已下载的 OpenGWAS AHR expression association 结果，计算 Kynurenine -> AHR expression 的探索性单 SNP Wald ratio。
- 复用既有 Kynurenine -> Myocarditis IVW MR 结果。
- 生成中介模式图、AHR eQTL locus Manhattan 图、心肌炎 AHR locus Manhattan 图、coloc 汇总图和两 SNP forest 图。

当前仍不能完成的部分：

- 完整 5-SNP Kynurenine -> AHR expression MR 仍不能完成：OpenGWAS AHR expression 完整 VCF 和 eQTLGen trans-eQTL 全量文件均已检查，5 个 Kynurenine 工具 SNP 中只有 `rs3184504` 有 AHR expression association。下面的 Kynurenine -> AHR expression 结果只能作为单 SNP 探索性结果，不能替代完整 5-SNP MR。

## Kynurenine -> AHR expression 探索性结果

| SNP | status | effect_allele_exposure | other_allele_exposure | beta_exposure | se_exposure | pval_exposure | effect_allele_outcome | other_allele_outcome | beta_outcome_aligned | se_outcome | pval_outcome | wald_beta | wald_se | wald_pval | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rs4843270 | missing_ahr_expression_association | A | C | 0.206 | 0.017 | 1.534e-35 | NA | NA | NA | NA | NA | NA | NA | NA | Not present in OpenGWAS AHR expression associations. |
| rs61825638 | missing_ahr_expression_association | T | C | 0.130 | 0.017 | 2.895e-14 | NA | NA | NA | NA | NA | NA | NA | NA | Not present in OpenGWAS AHR expression associations. |
| rs3184504 | ok | C | T | -0.108 | 0.015 | 2.839e-13 | C | T | 0.025 | 0.012 | 0.036 | -0.232 | 0.115 | 0.044 | Exploratory result: only this Kynurenine instrument is available for AHR expression; the complete 5-SNP MR remains not computable from the current AHR expression data. |
| rs6540080 | missing_ahr_expression_association | A | G | -0.108 | 0.016 | 5.012e-11 | NA | NA | NA | NA | NA | NA | NA | NA | Not present in OpenGWAS AHR expression associations. |
| rs10216901 | missing_ahr_expression_association | T | C | -0.093 | 0.015 | 1.461e-09 | NA | NA | NA | NA | NA | NA | NA | NA | Not present in OpenGWAS AHR expression associations. |

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
| Kynurenine -> Myocarditis | completed_existing_mr | Inverse variance weighted | 5 | 0.365 | 0.143 | 0.011 | 1.441 | Reused existing final5 MR result. |
| AHR expression -> Myocarditis | completed_two_snp_wald_where_available | Single-SNP Wald ratio | 2 | NA | NA | NA | NA | 5 outcome-SNP rows available across the three myocarditis datasets; see AHR_two_significant_snp_wald_ratios.tsv. |
| Kynurenine -> AHR expression | completed_exploratory_single_snp_not_5snp | Single-SNP Wald ratio | 1 | -0.232 | 0.115 | 0.044 | NA | Only rs3184504 was present in AHR expression data. OpenGWAS full VCF and eQTLGen trans-eQTL were checked; complete 5-SNP MR is still not computable. |
| Indirect effect: Kynurenine -> AHR -> Myocarditis | not_computed_for_full_mediation | Product of coefficients | NA | NA | NA | NA | NA | Full mediation requires the complete 5-SNP Kynurenine -> AHR expression estimate; the available one-SNP result is exploratory and not used as a definitive indirect effect. |

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
- `results/add_project/kynurenine_to_AHR_expression_single_snp_wald.tsv`
- `results/add_project/kynurenine_5snp_AHR_trans_eqtlgen_hits.tsv`
- `results/figures/AHR_two_snp_mediation_model.png`
- `results/figures/AHR_eqtlgen_AHR_locus_manhattan_two_snp.png`
- `results/figures/AHR_myocarditis_locus_manhattan_two_snp.png`
- `results/figures/AHR_coloc_summary_with_two_snp_context.png`
- `results/figures/AHR_two_snp_wald_ratio_forest.png`

## 复现命令

```bash
conda run -n wenhuai python AHR_myocarditis_gwas/scripts/15_add_docx_ahr_mediation_figures.py
```
