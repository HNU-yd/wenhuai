# AHR eQTLGen whole-blood cis-eQTL 与心肌炎 GWAS 的 MR 结果解释

## 核心结论

在三套并列主分析心肌炎 GWAS 中，AHR whole-blood cis-eQTL 工具变量未显示稳定且显著的遗传因果证据。

严格工具变量集合定义为：p < 5e-8，LD clumping r2 < 0.001，窗口 10000 kb。

在该严格集合下：

- BBJ 数据集仅匹配 1 个 SNP，只能进行 Wald ratio，结果不显著。
- Sakaue European 数据集可匹配 3 个 SNP，IVW OR 约为 1.10，结果不显著。
- FinnGen R12 数据集通过 rsID 匹配恢复 3 个 SNP，IVW OR 约为 0.92，结果不显著。

## 方向一致性

EUR 数据集多数设置下效应方向偏正，FinnGen 数据集多数设置下效应方向偏负，BBJ 数据集由于可用工具变量较少且置信区间较宽，估计不稳定。

因此，当前 MR 结果不能支持“遗传预测的 AHR 表达升高会显著增加或降低心肌炎风险”的结论。

## 后续分析重点

MR 结果阴性或方向不一致并不排除 AHR locus 在局部遗传调控中的作用。下一步应使用 full AHR cis-eQTL 和三套心肌炎 GWAS 的 AHR locus summary statistics 进行 coloc 共定位分析，判断 AHR eQTL 信号和心肌炎 GWAS 局部信号是否共享同一因果变异。
