
---

## 原始心肌炎 GWAS 文件检查

### 新增脚本

scripts/02_inspect_raw_outcome_gwas.sh

### 脚本目的

该脚本用于检查三套并列主分析心肌炎 GWAS 原始 summary statistics 文件：

1. Sakaue2021_BBJ_Myocarditis
2. Sakaue2021_EUR_Myocarditis
3. FinnGen_R12_I9_MYOCARD

本步骤不修改原始 GWAS 文件，只读取文件信息和前几行内容。

### 输入文件

data/raw/outcome/sakaue_2021_BBJ/hum0197.v3.BBJ.Myo.v1.zip  
data/raw/outcome/sakaue_2021_EUR/hum0197.v3.EUR.Myo.v1.zip  
data/raw/outcome/finngen_R12/finngen_R12_I9_MYOCARD.gz  

### 运行方式

bash scripts/02_inspect_raw_outcome_gwas.sh

### 输出文件

results/qc/raw_outcome_gwas_header_check.txt  
docs/file_manifest.tsv  
logs/02_inspect_raw_outcome_gwas.log  

### 输出文件怎么看

results/qc/raw_outcome_gwas_header_check.txt 用于查看每个原始 GWAS 的表头和前几行。重点关注变异 ID、染色体、位置、effect allele、other allele、beta、SE、P value、allele frequency、样本量、case 数和 control 数。

docs/file_manifest.tsv 用于记录原始文件完整性，包括文件路径、文件大小、sha256 校验值和文件类型。后续如果文件被替换、损坏或重复下载，可以通过 sha256 追踪。
