# 加项目.docx 任务摘录

原 docx 要求的核心事项：

- 目的：跑中介 MR。
- AHR GWAS ID：`eqtl-a-ENSG00000106546`。
- NCBI Gene 参考：AHR gene ID `196`，Ensembl ID `ENSG00000106546`。
- OpenGWAS 数据集检索：`https://opengwas.io/datasets/`。
- 指定 SNP：`rs17643734`、`rs59291726`；只跑两个显著 SNP。
- 补充：中介模式图、coloc、Manhattan plot 图。

本地执行说明：

- AHR expression -> Myocarditis 已用本地 AHR cis-eQTL 和心肌炎 locus 数据完成两 SNP Wald ratio。
- Kynurenine -> AHR expression 需要 full AHR expression GWAS 或 OpenGWAS JWT；当前本地 cis-eQTL 文件不足以计算该路径。
