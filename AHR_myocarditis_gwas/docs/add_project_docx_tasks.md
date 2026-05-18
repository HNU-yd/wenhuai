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
- Kynurenine -> AHR expression 已检查 OpenGWAS AHR expression 完整 VCF 和 eQTLGen trans-eQTL 全量文件；5 个 Kynurenine 工具 SNP 中只有 `rs3184504` 有 AHR expression association，因此只能输出探索性单 SNP Wald ratio，不能完成完整 5-SNP MR。
