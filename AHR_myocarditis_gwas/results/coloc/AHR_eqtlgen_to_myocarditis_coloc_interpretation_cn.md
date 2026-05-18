# AHR eQTLGen 与心肌炎 GWAS 的 coloc 共定位结果解释

## 判读标准

- PP.H4.abf 高：支持 AHR eQTL 和心肌炎 GWAS 共享同一因果变异。
- PP.H3.abf 高：两个性状在该区域都有信号，但更可能是不同因果变异。
- PP.H0/H1/H2 高：说明共定位证据不足，可能只有一个性状有信号或两个性状都缺乏明显信号。

常用经验阈值：PP.H4.abf > 0.8 可视为较强共定位证据；0.5–0.8 为中等证据；低于 0.5 通常不支持强共定位。

## 本次输出文件

- summary: /home/data1/wenhuai/AHR_myocarditis_gwas/results/coloc/AHR_eqtlgen_to_myocarditis_coloc_summary.tsv
- QC: /home/data1/wenhuai/AHR_myocarditis_gwas/results/qc/AHR_eqtlgen_to_myocarditis_coloc_qc.tsv

## 注意

FinnGen 若显示 skipped_missing_ncase_ncontrol，说明其 locus overlap 已生成，但还缺少 FinnGen I9_MYOCARD 的病例数和对照数。补充 `--finngen_ncase` 与 `--finngen_ncontrol` 后可直接重跑。
