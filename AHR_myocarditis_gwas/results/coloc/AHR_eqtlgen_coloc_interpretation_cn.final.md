# AHR eQTLGen 与心肌炎 GWAS coloc 共定位结果总结

## 1. 核心结论

三套并列主分析心肌炎 GWAS 中，均未观察到 AHR whole-blood cis-eQTL 与心肌炎 GWAS 在 AHR locus 共享同一因果变异的强证据。

当前经验判读标准为：

- PP.H4.abf > 0.8：强共定位证据；
- PP.H4.abf = 0.5–0.8：中等共定位证据；
- PP.H4.abf < 0.5：通常不支持强共定位。

## 2. 三套数据库结果

### Sakaue2021_BBJ_Myocarditis

- status: ok
- n_overlap: 5192
- PP.H1.abf: 0.6287
- PP.H3.abf: 0.3095
- PP.H4.abf: 0.06186
- interpretation: no_strong_colocalization

### Sakaue2021_EUR_Myocarditis

- status: ok
- n_overlap: 7383
- PP.H1.abf: 0.5961
- PP.H3.abf: 0.3557
- PP.H4.abf: 0.04827
- interpretation: no_strong_colocalization

### FinnGen_R12_I9_MYOCARD

- status: ok
- n_overlap: 7171
- PP.H1.abf: 0.6956
- PP.H3.abf: 0.2836
- PP.H4.abf: 0.02084
- interpretation: no_strong_colocalization

## 3. 生物学解释

AHR 在 eQTLGen whole blood 中具有非常强的 cis-eQTL 信号，但这些调控 AHR 表达的遗传变异并未与心肌炎 GWAS 的 AHR 区域局部信号形成强共定位。

因此，当前遗传证据不支持“外周血 AHR 表达遗传调控本身是心肌炎风险的直接驱动因子”。

这并不否定 AHR 通路在心肌炎中的作用。更合理的解释是：AHR 可能更多体现为疾病状态下的免疫代谢响应或下游通路活动，而不是由 AHR 基因表达 cis-eQTL 单独驱动的遗传易感机制。

## 4. 与 MR 结果的关系

MR 分析未观察到稳定显著的 AHR 表达遗传因果效应；coloc 分析进一步显示 AHR eQTL 与心肌炎 GWAS 局部信号缺乏强共定位。两者方向一致，均提示 AHR 表达遗传调控证据较弱。

后续文章中应避免写成“AHR 表达升高导致心肌炎风险升高”。更合适的表述是：遗传层面未支持 AHR 表达本身作为心肌炎风险的直接因果因子，但 AHR 通路仍可结合犬尿氨酸代谢 MR、单细胞和空间转录组结果作为疾病相关免疫代谢通路进行讨论。
