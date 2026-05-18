# Kynurenine–AHR 第二层 V0 单细胞 baseline 项目说明

## 1. 这个 README 描述的是什么

本 README 描述 `/home/data1/wenhuai` 项目下第二层 V0 单细胞分析流程，包括：

1. 项目在做什么。
2. 当前数据结构是什么。
3. `/src` 和 `/scripts` 下每个文件的作用。
4. 每个脚本的输入是什么。
5. 每个脚本输出什么文件。
6. 输出文件里面有什么。
7. 输出文件后续用来干什么。
8. 如何从原始公开数据跑到 V0 结果。
9. 如何检查最终结果。

本 README 对应当前版本：

```text
V0-alpha / V0-baseline
````

---

## 2. 这个项目在干什么

本项目是儿童心肌炎 / MIS-C / 病毒性心肌炎机制研究的第二层分析。

整体科学问题是：

```text
为什么有些儿童心肌炎患儿可以快速恢复，而有些会进展为 severe / fulminant myocarditis？
```

第二层 V0 当前要回答的是：

```text
公开 PBMC 单细胞数据中，是否存在 Kynurenine–AHR-high inflammatory myeloid state？
这种状态是否主要集中在 monocyte / DC / inflammatory myeloid cells？
这种状态是否在 severe / MIS-C_MYO / acute / CVB3 myocarditis 相关组别中富集？
```

注意：

```text
scRNA-seq 不能直接测量 Kynurenine 浓度。
```

所以本项目中所有 Kynurenine 相关结果都表述为：

```text
Kynurenine metabolism potential
Kyn-AHR pathway activity
AHR response / regulon proxy
```

而不是：

```text
Kynurenine level
```

---

## 3. 当前数据集定位

当前进入第二层 V0 的核心数据集是：

```text
GSE183716
CNP0005824
GSE166489
GSE167029
```

扩展数据集：

```text
GSE180045
```

暂不进入第二层 V0 的空间数据：

```text
STT0000127
```

### 3.1 核心数据集

| 数据集        | 当前格式                                   | 当前状态       | V0 作用                                            |
| ---------- | -------------------------------------- | ---------- | ------------------------------------------------ |
| GSE183716  | 4 个 10x h5                             | 可直接读入      | smoke test，小规模 acute/recovery 验证                 |
| CNP0005824 | 14 个标准 10x mtx 目录                      | 已通过复制方式标准化 | CVB3 / IVIG / control 模型验证                       |
| GSE166489  | 18 个 10x filtered_feature_bc_matrix 目录 | 已解包        | MIS-C severe / moderate / recovered / healthy 验证 |
| GSE167029  | 27 个 raw_feature_bc_matrix h5          | 已解包        | MIS-C_MYO 主数据                                    |

### 3.2 扩展和空间数据

| 数据集        | 格式                             | 当前处理策略                      |
| ---------- | ------------------------------ | --------------------------- |
| GSE180045  | 10x mtx 目录                     | 后续单独作为 ICI myocarditis 扩展分析 |
| STT0000127 | GEM / TissueCut / bin50.gem.gz | 留到第三层空间转录组分析，不进第二层 PBMC V0  |

---

## 4. 项目文件结构

项目根目录：

```text
/home/data1/wenhuai
```

目录结构：

```text
/home/data1/wenhuai
├── README.md
├── scripts/
│   ├── 01_discover_local_data.py
│   ├── 02_standardize_prefixed_10x_copy.py
│   ├── 03_build_geo_sample_table.py
│   ├── 04_prepare_h5ad_from_inventory.py
│   ├── 05_concat_h5ad_on_disk.py
│   ├── 06_run_v0_scanpy.py
│   ├── 07_check_v0_outputs.py
│   └── run_v0_core.sh
│
├── src/
│   └── kyn_ahr_v0/
│       ├── __init__.py
│       ├── gene_sets.py
│       └── metadata_rules.py
│
├── data/
│   ├── CNP0005824/
│   ├── GSE166489/
│   ├── GSE167029/
│   ├── GSE180045/
│   ├── GSE183716/
│   ├── STT0000127/
│   ├── metadata/
│   ├── h5ad_per_sample/
│   └── h5ad_merged/
│
├── results/
│   └── v0/
│       ├── GSE183716/
│       ├── CNP0005824/
│       ├── GSE166489/
│       └── GSE167029/
│
└── logs/
```

---

## 5. 数据结构说明

### 5.1 CNP0005824

原始数据是前缀式 10x mtx 文件：

```text
data/CNP0005824/Single_Cell/CSE0000419/
├── Control_matrix.mtx.gz
├── Control_barcodes.tsv.gz
├── Control_features.tsv.gz
├── CVB3d1_matrix.mtx.gz
├── CVB3d1_barcodes.tsv.gz
├── CVB3d1_features.tsv.gz
└── ...
```

标准化后复制到：

```text
data/CNP0005824/tenx_standard/
├── Control/
│   ├── matrix.mtx.gz
│   ├── barcodes.tsv.gz
│   └── features.tsv.gz
├── CVB3d1/
│   ├── matrix.mtx.gz
│   ├── barcodes.tsv.gz
│   └── features.tsv.gz
├── CVB3d2/
└── ...
```

后续 Scanpy 只读取：

```text
data/CNP0005824/tenx_standard/<sample>/
```

---

### 5.2 GSE166489

解包后结构：

```text
data/GSE166489/suppl_unpacked/
├── GSM5073055_P1.1_filtered_feature_bc_matrix/
├── GSM5073056_P2.1_filtered_feature_bc_matrix/
├── ...
├── GSM5073064_C.HD1_filtered_feature_bc_matrix/
└── ...
```

每个目录内应有：

```text
matrix.mtx.gz
barcodes.tsv.gz
features.tsv.gz
```

---

### 5.3 GSE167029

解包后结构：

```text
data/GSE167029/suppl_unpacked/
├── GSM5091035_C24_raw_feature_bc_matrix.h5
├── GSM5091036_C26_raw_feature_bc_matrix.h5
├── ...
└── GSM5091061_P55_raw_feature_bc_matrix.h5
```

这是 10x h5 格式，直接由 Scanpy 读取。

---

### 5.4 GSE183716

结构：

```text
data/GSE183716/suppl/
├── GSE183716_Sample1_GEXFB_filtered_feature_bc_matrix.h5
├── GSE183716_Sample2_GEXFB_filtered_feature_bc_matrix.h5
├── GSE183716_Sample3_GEXFB_filtered_feature_bc_matrix.h5
└── GSE183716_Sample4_GEXFB_filtered_feature_bc_matrix.h5
```

这是 10x h5 格式，适合 smoke test。

---

## 6. `/src` 文件说明

---

### 6.1 `src/kyn_ahr_v0/__init__.py`

它的作用是声明 `kyn_ahr_v0` 是一个 Python package。

内容很简单：

```text
__version__ = "0.2.0"
```

用途：

```text
让 scripts/ 下的脚本可以通过 sys.path 调用 src/kyn_ahr_v0 中的模块。
```

输出文件：

```text
无直接输出文件。
```

---

### 6.2 `src/kyn_ahr_v0/gene_sets.py`

它的作用是定义 V0 分析中使用的基因集。

主要包含：

```text
GENE_SETS
CELLTYPE_MARKERS
MYELOID_CELLTYPES
```

#### 1. `GENE_SETS`

用于计算 pathway score。

包括：

```text
KYN_metabolism_score
AHR_response_score
AHR_regulon_proxy_score
myeloid_inflammation_score
chemotaxis_score
antigen_presentation_score
```

这些 score 最后会写入：

```text
adata.obs
results/v0/<dataset>/<dataset>.obs_after_scores.tsv
results/v0/<dataset>/<dataset>.score_pseudobulk_by_sample_celltype.tsv
```

#### 2. `CELLTYPE_MARKERS`

用于粗粒度 PBMC 细胞注释。

包括：

```text
T_cell
CD8_T
NK
B_cell
Plasmablast
Monocyte
CD16_Monocyte
DC
pDC
Platelet
```

粗注释结果写入：

```text
adata.obs["cell_type_v0"]
```

#### 3. `MYELOID_CELLTYPES`

定义哪些 cell type 属于髓系细胞：

```text
Monocyte
CD16_Monocyte
DC
pDC
```

后续用于输出髓系专属表：

```text
<dataset>.myeloid_score_pseudobulk.tsv
<dataset>.DE_myeloid_disease_vs_control.wilcoxon.tsv
```

输出文件：

```text
无直接输出文件。
```

它只提供基因集，被 `06_run_v0_scanpy.py` 调用。

---

### 6.3 `src/kyn_ahr_v0/metadata_rules.py`

它的作用是统一推断样本 metadata。

主要功能：

```text
1. 清洗 sample_id。
2. 从路径中提取 GSM 编号。
3. 解析 GEO family.soft 文件。
4. 根据文件名和 GEO metadata 推断 group / severity / stage。
5. 生成每个样本统一 metadata。
```

它会生成这些字段：

```text
dataset
source_file
gsm
sample_id
patient_id
group
severity
stage
batch
geo_title
geo_source_name
geo_characteristics
```

这些字段会写入：

```text
data/h5ad_per_sample/<dataset>/<sample>.h5ad 的 adata.obs
data/h5ad_per_sample/<dataset>/<sample>.metadata.tsv
data/metadata/sample_sheet.auto.tsv
```

输出文件：

```text
无直接输出文件。
```

但它会被下面脚本调用：

```text
03_build_geo_sample_table.py
04_prepare_h5ad_from_inventory.py
```

---

## 7. `/scripts` 文件说明

---

# 7.1 `01_discover_local_data.py`

## 这个文件是做什么的

这个脚本用于扫描 `data/` 目录，发现哪些文件可以进入单细胞分析流程。

它会识别：

```text
h5ad
h5_or_hdf5
seurat_rds
10x_mtx_dir
10x_mtx_matrix
10x_mtx_barcodes
10x_mtx_features
archive_tar
archive_zip
geo_soft
geo_series_matrix
```

它默认扫描：

```text
CNP0005824
GSE166489
GSE167029
GSE180045
GSE183716
```

不默认扫描：

```text
STT0000127
```

因为 STT0000127 是空间 GEM 数据，不属于第二层 PBMC V0。

## 怎么运行

```bash
python scripts/01_discover_local_data.py \
  --datasets CNP0005824 GSE166489 GSE167029 GSE183716
```

## 输入是什么

输入是本地数据目录：

```text
data/CNP0005824/
data/GSE166489/
data/GSE167029/
data/GSE183716/
```

## 输出到哪里

核心输出：

```text
data/metadata/local_data_inventory.tsv
```

## 输出文件解析

### `data/metadata/local_data_inventory.tsv`

这是本地数据清单表。

每一行代表一个被识别到的文件或目录。

字段包括：

```text
dataset
record_type
path
parent
name
file_type
size_bytes
size_gb
```

#### 字段解释

`dataset`

```text
数据集名称，例如 GSE166489、GSE167029、CNP0005824。
```

`record_type`

```text
记录类型。
file 表示这是一个文件。
directory 表示这是一个目录，例如标准 10x mtx 目录。
```

`path`

```text
文件或目录的完整路径。
后续 04_prepare_h5ad_from_inventory.py 会根据这个 path 读取表达矩阵。
```

`parent`

```text
该文件或目录的上级目录。
```

`name`

```text
文件名或目录名。
```

`file_type`

```text
脚本识别出的文件类型。
```

常见值：

```text
10x_mtx_dir      可以直接用 scanpy.read_10x_mtx 读取
h5_or_hdf5       可以尝试用 scanpy.read_10x_h5 读取
h5ad             已经是 AnnData 格式
seurat_rds       Seurat RDS，当前 V0 Python 流程暂不直接读取
geo_soft         GEO soft metadata
geo_series_matrix GEO series matrix metadata
archive_tar      压缩包
```

`size_bytes` / `size_gb`

```text
文件大小。
目录类型通常是 0。
```

#### 这个文件用来干什么

它是后续转换 h5ad 的输入索引。

`04_prepare_h5ad_from_inventory.py` 会从这里筛选：

```text
10x_mtx_dir
h5_or_hdf5
h5ad
```

然后转换成单样本 h5ad。

---

# 7.2 `02_standardize_prefixed_10x_copy.py`

## 这个文件是做什么的

这个脚本用于处理 CNP0005824 这类前缀式 10x mtx 文件。

原始文件是：

```text
CVB3d1_matrix.mtx.gz
CVB3d1_barcodes.tsv.gz
CVB3d1_features.tsv.gz
```

但 Scanpy 默认需要标准 10x 目录：

```text
matrix.mtx.gz
barcodes.tsv.gz
features.tsv.gz
```

所以这个脚本会复制并统一命名。

## 怎么运行

```bash
python scripts/02_standardize_prefixed_10x_copy.py \
  --datasets CNP0005824
```

也可以对多个数据集运行：

```bash
python scripts/02_standardize_prefixed_10x_copy.py \
  --datasets CNP0005824 GSE166489 GSE167029
```

## 输入是什么

主要输入：

```text
data/CNP0005824/Single_Cell/CSE0000419/
```

里面有：

```text
*_matrix.mtx.gz
*_barcodes.tsv.gz
*_features.tsv.gz
```

## 输出到哪里

输出到：

```text
data/CNP0005824/tenx_standard/<sample>/
```

## 输出文件解析

### `data/CNP0005824/tenx_standard/<sample>/matrix.mtx.gz`

这是表达矩阵。

内容是：

```text
基因 × 细胞 的 sparse count matrix。
```

用途：

```text
被 scanpy.read_10x_mtx() 读取，构建 AnnData.X。
```

---

### `data/CNP0005824/tenx_standard/<sample>/barcodes.tsv.gz`

这是细胞 barcode 文件。

内容是：

```text
每一行一个 cell barcode。
```

用途：

```text
作为 AnnData.obs_names 的基础。
```

---

### `data/CNP0005824/tenx_standard/<sample>/features.tsv.gz`

这是基因注释文件。

内容通常包含：

```text
gene_id
gene_symbol
feature_type
```

用途：

```text
作为 AnnData.var 的基础。
```

---

### 标准化后的目录有什么用

标准化后，CNP0005824 就可以像普通 10x 数据一样被读取：

```python
scanpy.read_10x_mtx("data/CNP0005824/tenx_standard/CVB3d1")
```

这个脚本不使用软链接，不破坏原始文件，只复制一份标准输入。

---

# 7.3 `03_build_geo_sample_table.py`

## 这个文件是做什么的

这个脚本用于解析 GEO 的 `family.soft` 文件，提取每个 GSM 样本的注释信息。

它会读取：

```text
data/<dataset>/soft/*family.soft.gz
```

并自动解压为：

```text
data/<dataset>/soft/*family.soft
```

然后提取：

```text
GSM 编号
Sample title
source_name
characteristics
```

## 怎么运行

```bash
python scripts/03_build_geo_sample_table.py
```

## 输入是什么

输入是 GEO soft 文件：

```text
data/GSE166489/soft/GSE166489_family.soft.gz
data/GSE167029/soft/GSE167029_family.soft.gz
data/GSE183716/soft/GSE183716_family.soft.gz
```

## 输出到哪里

核心输出：

```text
data/metadata/geo_sample_table.tsv
```

## 输出文件解析

### `data/metadata/geo_sample_table.tsv`

这是 GEO 样本注释表。

字段包括：

```text
dataset
gsm
title
source_name
characteristics
```

#### 字段解释

`dataset`

```text
数据集名称。
```

`gsm`

```text
GEO 样本编号，例如 GSM5073055。
```

`title`

```text
GEO 中的 Sample title。
通常包含样本名、组别或患者编号。
```

`source_name`

```text
样本来源，例如 PBMC、whole blood、control 等。
```

`characteristics`

```text
GEO 中的样本属性信息。
可能包含 disease status、severity、timepoint、age、sex 等。
```

#### 这个文件用来干什么

它被 `04_prepare_h5ad_from_inventory.py` 调用，用于辅助推断：

```text
group
severity
stage
patient_id
sample_id
```

如果文件名中没有完整分组信息，就尽量从 GEO soft metadata 中补充。

---

# 7.4 `04_prepare_h5ad_from_inventory.py`

## 这个文件是做什么的

这个脚本根据 `local_data_inventory.tsv` 把每个样本转换成单样本 h5ad。

它支持：

```text
10x_mtx_dir
h5_or_hdf5
h5ad
```

不处理：

```text
GEM
soft
series_matrix
tar
zip
```

## 怎么运行

```bash
python scripts/04_prepare_h5ad_from_inventory.py \
  --datasets GSE183716 CNP0005824 GSE166489 GSE167029 \
  --workers 12
```

## 输入是什么

主要输入：

```text
data/metadata/local_data_inventory.tsv
data/metadata/geo_sample_table.tsv
```

其中：

```text
local_data_inventory.tsv 负责告诉脚本表达矩阵在哪里。
geo_sample_table.tsv 负责提供 GEO 样本注释。
```

## 输出到哪里

主要输出：

```text
data/h5ad_per_sample/<dataset>/<sample>.h5ad
data/h5ad_per_sample/<dataset>/<sample>.metadata.tsv
data/metadata/prepare_h5ad_from_inventory_report.tsv
data/metadata/sample_sheet.auto.tsv
```

---

## 输出文件解析

### `data/h5ad_per_sample/<dataset>/<sample>.h5ad`

这是每个样本单独转换出来的 AnnData 文件。

例如：

```text
data/h5ad_per_sample/GSE166489/GSM5073055_P1.1.h5ad
data/h5ad_per_sample/GSE167029/GSM5091035_C24.h5ad
```

里面包含：

```text
adata.X
adata.obs
adata.var
adata.uns["source_file"]
```

#### `adata.X`

这是原始 count matrix。

内容：

```text
细胞 × 基因 的表达矩阵。
```

用途：

```text
后续 QC、normalize、log1p、score、PCA 都基于它。
```

#### `adata.obs`

这是细胞 metadata。

包含：

```text
dataset
gsm
sample_id
patient_id
group
severity
stage
batch
geo_title
geo_source_name
geo_characteristics
n_genes_by_counts
total_counts
pct_counts_mt
```

用途：

```text
用于分组比较、batch correction、sample-level pseudobulk 和 patient-level 聚合。
```

#### `adata.var`

这是基因 metadata。

包含：

```text
gene symbols
gene ids
mt 标记
```

用途：

```text
用于识别 mitochondrial genes、pathway gene matching 和 feature selection。
```

#### `adata.uns["source_file"]`

记录该 h5ad 来源于哪个原始文件或目录。

用途：

```text
方便追踪样本来源。
```

---

### `data/h5ad_per_sample/<dataset>/<sample>.metadata.tsv`

这是单样本 metadata 表。

每个文件通常只有一行。

字段包括：

```text
dataset
source_file
gsm
sample_id
patient_id
group
severity
stage
batch
geo_title
geo_source_name
geo_characteristics
```

用途：

```text
用于人工检查某个样本的自动分组是否正确。
```

---

### `data/metadata/prepare_h5ad_from_inventory_report.tsv`

这是转换报告。

字段包括：

```text
status
dataset
file_type
path
out
n_obs
n_vars
error
```

#### 字段解释

`status`

```text
ok 表示转换成功。
skip 表示目标 h5ad 已存在，跳过。
error 表示转换失败。
```

`dataset`

```text
数据集名称。
```

`file_type`

```text
原始输入类型，例如 10x_mtx_dir 或 h5_or_hdf5。
```

`path`

```text
输入文件或目录。
```

`out`

```text
输出 h5ad 路径。
```

`n_obs`

```text
细胞数。
```

`n_vars`

```text
基因数。
```

`error`

```text
如果失败，记录失败原因。
```

用途：

```text
检查每个样本是否成功转成 h5ad。
```

---

### `data/metadata/sample_sheet.auto.tsv`

这是自动生成的样本总表。

字段包括：

```text
dataset
source_file
gsm
sample_id
patient_id
group
severity
stage
batch
geo_title
geo_source_name
geo_characteristics
```

用途：

```text
用于整体检查所有样本分组。
后续如果自动 metadata 不准确，可以基于这个表制作 sample_sheet.manual.tsv 进行覆盖。
```

---

# 7.5 `05_concat_h5ad_on_disk.py`

## 这个文件是做什么的

这个脚本把每个数据集内多个单样本 h5ad 合并成一个数据集级 h5ad。

例如：

```text
data/h5ad_per_sample/GSE166489/*.h5ad
```

会合并成：

```text
data/h5ad_merged/GSE166489.raw_merged.h5ad
```

## 怎么运行

推荐使用 memory 模式：

```bash
python scripts/05_concat_h5ad_on_disk.py \
  --datasets GSE183716 CNP0005824 GSE166489 GSE167029 \
  --overwrite \
  --mode memory
```

## 输入是什么

输入是：

```text
data/h5ad_per_sample/<dataset>/*.h5ad
```

## 输出到哪里

输出到：

```text
data/h5ad_merged/
```

---

## 输出文件解析

### `data/h5ad_merged/<dataset>.raw_merged.h5ad`

这是每个数据集合并后的 h5ad。

例如：

```text
data/h5ad_merged/GSE183716.raw_merged.h5ad
data/h5ad_merged/CNP0005824.raw_merged.h5ad
data/h5ad_merged/GSE166489.raw_merged.h5ad
data/h5ad_merged/GSE167029.raw_merged.h5ad
```

里面包含：

```text
所有样本合并后的 adata.X
所有细胞的 adata.obs
所有基因的 adata.var
sample_file 字段
```

用途：

```text
作为 06_run_v0_scanpy.py 的输入。
```

#### 为什么推荐 `--mode memory`

当前部分 h5ad 内部是 `csc` sparse matrix，而 `anndata.experimental.concat_on_disk` 对 csc 支持不好。

所以推荐使用：

```bash
--mode memory
```

脚本会自动：

```text
1. 逐个读取 per-sample h5ad。
2. 把 X 和 layers 转成 CSR。
3. 合并。
4. 先写临时文件。
5. 成功后再替换正式文件。
```

这样可以避免：

```text
Concat of following not supported: ['csc', ...]
unable to truncate a file which is already open
```

---

# 7.6 `06_run_v0_scanpy.py`

## 这个文件是做什么的

`06_run_v0_scanpy.py` 是 V0 主分析脚本。

它的作用是把 `05_concat_h5ad_on_disk.py` 合并出来的：

```text
data/h5ad_merged/<dataset>.raw_merged.h5ad
```

做完整 Scanpy V0 分析：

```text
QC
标准化
log1p
Kyn-AHR pathway scoring
PBMC marker scoring
cell_type_v0 粗注释
HVG
PCA
Harmony
neighbors
Leiden
UMAP
pseudobulk
top Kyn-AHR-high cells
快速差异分析
```

## 怎么运行

示例：

先跑小的：
```bash
python scripts/06_run_v0_scanpy.py \
  --input /home/data1/wenhuai/data/h5ad_merged/GSE183716.raw_merged.h5ad \
  --dataset GSE183716 \
  --threads 20 \
  2>&1 | tee logs/06_run_v0_GSE183716.log

python scripts/06_run_v0_scanpy.py \
  --input /home/data1/wenhuai/data/h5ad_merged/GSE166489.raw_merged.h5ad \
  --dataset GSE166489 \
  --threads 20 \
  2>&1 | tee logs/06_run_v0_GSE166489.log
```

再跑大的：
```bash
python scripts/06_run_v0_scanpy.py \
  --input /home/data1/wenhuai/data/h5ad_merged/CNP0005824.raw_merged.h5ad \
  --dataset CNP0005824 \
  --threads 20 \
  2>&1 | tee logs/06_run_v0_CNP0005824.log

python scripts/06_run_v0_scanpy.py \
  --input /home/data1/wenhuai/data/h5ad_merged/GSE167029.raw_merged.h5ad \
  --dataset GSE167029 \
  --threads 20 \
  2>&1 | tee logs/06_run_v0_GSE167029.log
```
```bash
python scripts/06_run_v0_scanpy.py \
  --input /home/data1/wenhuai/data/h5ad_merged/GSE166489.raw_merged.h5ad \
  --dataset GSE166489 \
  --threads 8 \
  --min_genes 200 \
  --max_mito 20 \
  --n_hvg 3000
```

## 输入是什么

输入：

```text
data/h5ad_merged/<dataset>.raw_merged.h5ad
```

例如：

```text
data/h5ad_merged/GSE166489.raw_merged.h5ad
data/h5ad_merged/GSE167029.raw_merged.h5ad
```

## 输出到哪里

`06_run_v0_scanpy.py` 跑完后，核心输出都在：

```text
results/v0/<dataset>/
```

例如：

```text
results/v0/GSE166489/
results/v0/GSE167029/
```

---

## 输出文件解析

### 1. h5ad 类输出

---

### `<dataset>.01_qc.h5ad`

这是 QC 后的数据。

已经做了：

```text
低基因数细胞过滤
高线粒体比例细胞过滤
低表达基因过滤
```

用途：

```text
看过滤后还剩多少细胞、多少基因。
如果这个文件细胞数掉得太多，说明 --min_genes 或 --max_mito 可能太严格。
```

---

### `<dataset>.02_v0_processed.h5ad`

这是完成主分析流程后的中间版本。

里面通常已经有：

```text
normalize/log1p 后的表达矩阵
KYN/AHR/myeloid 等 score
cell_type_v0 粗注释
PCA
Harmony embedding
neighbors
Leiden cluster
UMAP
```

用途：

```text
用于检查主流程是否正常。
例如 UMAP 有没有明显 batch effect、Leiden 分群是否合理、score 是否成功写入。
```

---

### `<dataset>.final_v0.h5ad`

这是 V0 的最终 h5ad。

它包含：

```text
QC 后细胞
标准化表达矩阵
pathway scores
marker scores
cell_type_v0
leiden_v0
PCA
Harmony
UMAP
neighbors graph
metadata
```

用途：

```text
这是后续 V1 scVI / MIL 或人工复查优先使用的标准结果文件。
```

如果后面要进入 V1，可以从这里继续提取：

```text
cell expression
cell metadata
patient_id
group
severity
stage
Kyn-AHR scores
cell_type_v0
```

---

### 2. score 和 metadata 类输出

---

### `<dataset>.obs_after_scores.tsv`

这是每个细胞的 metadata 和 score 表。

每一行是一个细胞。

主要字段包括：

```text
dataset
gsm
sample_id
patient_id
group
severity
stage
batch
cell_type_v0
KYN_metabolism_score
AHR_response_score
AHR_regulon_proxy_score
myeloid_inflammation_score
chemotaxis_score
antigen_presentation_score
Kyn_AHR_myeloid_score
```

用途：

```text
用于在细胞级别检查 Kyn-AHR score。
也可以按 sample / group / severity 自己重新聚合。
```

常见查看方式：

```bash
head -5 results/v0/GSE166489/GSE166489.obs_after_scores.tsv
```

---

### `<dataset>.score_gene_presence.tsv`

这是每个 score gene set 的基因命中报告。

字段包括：

```text
score
n_requested
n_present
genes_present
```

它说明：

```text
某个 gene set 原本需要多少基因。
当前数据集中实际匹配到了多少基因。
具体匹配到了哪些基因。
```

用途：

```text
判断 score 是否可靠。
如果 n_present 很低，说明该 score 在当前数据集解释力有限。
```

例如：

```text
KYN_metabolism_score 如果只命中 1 个基因，就不能过度解释。
AHR_response_score 如果命中多个 AHR target genes，则结果更可信。
```

---

### 3. 细胞比例类输出

---

### `<dataset>.celltype_proportions.tsv`

这是样本级细胞类型比例表。

字段包括：

```text
dataset
sample_id
patient_id
group
severity
stage
cell_type_v0
n_cells
sample_total_cells
fraction
```

它的含义是：

```text
每个 sample 中每种 cell_type_v0 占多少比例。
```

用途：

```text
比较不同组别中 Monocyte、CD16_Monocyte、DC、NK、T、B 等细胞比例是否变化。
```

例如可以看：

```text
MIS-C_MYO 是否有更多 inflammatory monocyte。
severe 是否有更多 myeloid cells。
CVB3 组是否有 monocyte/DC 富集。
```

---

### 4. pathway pseudobulk 类输出

---

### `<dataset>.score_pseudobulk_by_sample_celltype.tsv`

这是 V0 最重要的表之一。

它是按：

```text
sample_id + cell_type_v0
```

聚合后的 pathway score 平均值。

字段包括：

```text
dataset
sample_id
patient_id
group
severity
stage
cell_type_v0
KYN_metabolism_score
AHR_response_score
AHR_regulon_proxy_score
myeloid_inflammation_score
chemotaxis_score
antigen_presentation_score
Kyn_AHR_myeloid_score
```

用途：

```text
比较不同组、不同细胞类型中的 Kyn-AHR pathway activity。
```

重点看：

```text
Monocyte
CD16_Monocyte
DC
pDC
```

如果这些髓系细胞在 severe / MIS-C_MYO / acute / CVB3 组中：

```text
KYN_metabolism_score 高
AHR_response_score 高
myeloid_inflammation_score 高
Kyn_AHR_myeloid_score 高
```

就支持 Kyn-AHR-high inflammatory myeloid state 的存在。

---

### `<dataset>.myeloid_score_pseudobulk.tsv`

这是只保留髓系细胞的 pseudobulk score 表。

保留的 cell type 包括：

```text
Monocyte
CD16_Monocyte
DC
pDC
```

用途：

```text
专门判断 Kyn-AHR signal 是否集中在髓系细胞中。
```

这是后续写第二层结果时最应该优先看的表之一。

---

### 5. Kyn-AHR-high cells 类输出

---

### `<dataset>.top10pct_Kyn_AHR_high_cells.tsv`

这是 Kyn_AHR_myeloid_score 排名前 10% 的细胞列表。

每一行是一个细胞。

包含：

```text
dataset
sample_id
patient_id
group
severity
stage
cell_type_v0
leiden_v0
KYN_metabolism_score
AHR_response_score
AHR_regulon_proxy_score
myeloid_inflammation_score
chemotaxis_score
antigen_presentation_score
Kyn_AHR_myeloid_score
```

用途：

```text
定位 Kyn-AHR-high cells 来自哪些样本、哪些组、哪些细胞类型、哪些 Leiden cluster。
```

如果 top cells 大量来自：

```text
Monocyte
CD16_Monocyte
DC
```

并且集中于 disease/severe/MIS-C_MYO/CVB3 相关组，说明信号方向符合预期。

---

### `<dataset>.top10pct_Kyn_AHR_high_summary.tsv`

这是 top 10% Kyn-AHR-high cells 的汇总表。

字段包括：

```text
dataset
group
severity
stage
cell_type_v0
n_top10pct_cells
```

用途：

```text
快速判断 Kyn-AHR-high cells 是否主要来自 severe / MIS-C_MYO / acute / CVB3 相关组。
```

推荐首先查看这个文件：

```bash
head -40 results/v0/GSE166489/GSE166489.top10pct_Kyn_AHR_high_summary.tsv
head -40 results/v0/GSE167029/GSE167029.top10pct_Kyn_AHR_high_summary.tsv
```

---

### 6. 差异基因类输出

---

### `<dataset>.DE_disease_vs_control.wilcoxon.tsv`

这是所有细胞层面的 disease vs control 快速 Wilcoxon 差异基因结果。

字段通常包括：

```text
names
scores
logfoldchanges
pvals
pvals_adj
```

用途：

```text
快速筛查疾病组全局差异。
```

注意：

```text
这是 V0-alpha 快速 DEG，不是最终论文级 pseudobulk DESeq2 / edgeR。
```

---

### `<dataset>.DE_myeloid_disease_vs_control.wilcoxon.tsv`

这是只在髓系细胞中做 disease vs control 的快速 Wilcoxon 差异基因结果。

用途：

```text
重点检查髓系细胞中是否出现炎症和 AHR 相关基因变化。
```

重点关注：

```text
S100A8
S100A9
IL1B
TNF
CXCL8
CCL2
NFKBIA
STAT1
IRF1
AHR
AHRR
TIPARP
CYP1A1
CYP1B1
```

如果这些基因在 myeloid disease vs control 中靠前，说明髓系炎症状态明显。

---

### 7. 图像输出

---

### `figures/*.png`

UMAP 图保存在：

```text
results/v0/<dataset>/figures/
```

常见图包括：

```text
umap_<dataset>_group.png
umap_<dataset>_severity.png
umap_<dataset>_stage.png
umap_<dataset>_sample_id.png
umap_<dataset>_patient_id.png
umap_<dataset>_cell_type_v0.png
umap_<dataset>_leiden_v0.png
umap_<dataset>_KYN_metabolism_score.png
umap_<dataset>_AHR_response_score.png
umap_<dataset>_AHR_regulon_proxy_score.png
umap_<dataset>_myeloid_inflammation_score.png
umap_<dataset>_chemotaxis_score.png
umap_<dataset>_antigen_presentation_score.png
umap_<dataset>_Kyn_AHR_myeloid_score.png
```

用途：

```text
检查 UMAP 分群、batch effect、细胞类型分布、pathway score 空间分布。
```

如果 `sample_id` 图呈现明显分离，而 `group` 或 `cell_type_v0` 不清楚，可能说明 batch effect 仍然存在。

如果 `Kyn_AHR_myeloid_score` 高值集中在 myeloid cluster，则符合预期。

---

# 7.7 `07_check_v0_outputs.py`

## 这个文件是做什么的

这个脚本用于快速检查 V0 输出是否完整。

它会检查：

```text
每个数据集有多少 per-sample h5ad
是否生成了 merged h5ad
是否生成了 final_v0.h5ad
是否生成了 top10pct_Kyn_AHR_high_summary.tsv
```

## 怎么运行

```bash
python scripts/07_check_v0_outputs.py \
  --datasets GSE183716 GSE166489 CNP0005824 GSE167029
```

## 输入是什么

输入来自：

```text
data/h5ad_per_sample/
data/h5ad_merged/
results/v0/
```

## 输出到哪里

这个脚本主要输出到终端。

如果用 tee 保存：

```bash
python scripts/07_check_v0_outputs.py \
  --datasets GSE183716 GSE166489 CNP0005824 GSE167029 \
  2>&1 | tee logs/07_check_v0_outputs.log
```

则日志会保存到：

```text
logs/07_check_v0_outputs.log
```

## 输出文件解析

### `logs/07_check_v0_outputs.log`

这是检查日志。

里面会显示：

```text
h5ad_per_sample 每个数据集有多少 h5ad
h5ad_merged 是否存在
results/v0/<dataset>/<dataset>.final_v0.h5ad 是否存在
top10pct_Kyn_AHR_high_summary.tsv 前 20 行
```

用途：

```text
快速确认整个 V0 流程是否跑完。
```

---

# 7.8 `run_v0_core.sh`

## 这个文件是做什么的

这是核心 V0 的一键运行脚本。

它会依次执行：

```text
1. 复制标准化前缀式 10x。
2. 扫描本地数据。
3. 解析 GEO sample metadata。
4. 转换 per-sample h5ad。
5. 合并 dataset-level h5ad。
6. 分别运行 V0 Scanpy 分析。
7. 检查输出。
```

## 怎么运行

```bash
bash scripts/run_v0_core.sh
```

## 输入是什么

输入是当前项目下的原始数据：

```text
data/CNP0005824/
data/GSE166489/
data/GSE167029/
data/GSE183716/
```

## 输出到哪里

输出分布在：

```text
data/metadata/
data/h5ad_per_sample/
data/h5ad_merged/
results/v0/
logs/
```

## 输出文件解析

### `logs/run_v0_core.sh` 相关日志

如果用 tee 运行：

```bash
bash scripts/run_v0_core.sh 2>&1 | tee logs/run_v0_core.log
```

会得到：

```text
logs/run_v0_core.log
```

用途：

```text
记录完整运行过程，便于排查哪一步失败。
```

### 其他输出

这个脚本本身不产生新的独立数据格式，而是依次触发前面脚本产生：

```text
local_data_inventory.tsv
geo_sample_table.tsv
sample_sheet.auto.tsv
per-sample h5ad
merged h5ad
V0 result files
```

所以具体输出解释见前面各脚本说明。

---

## 8. 推荐运行顺序

### 8.1 准备环境

```bash
cd /home/data1/wenhuai
conda activate kynahr_v0

mkdir -p logs data/metadata data/h5ad_per_sample data/h5ad_merged results/v0
```

---

### 8.2 数据发现和转换

```bash
python scripts/02_standardize_prefixed_10x_copy.py \
  --datasets CNP0005824 GSE166489 GSE167029

python scripts/01_discover_local_data.py \
  --datasets CNP0005824 GSE166489 GSE167029 GSE183716

python scripts/03_build_geo_sample_table.py

python scripts/04_prepare_h5ad_from_inventory.py \
  --datasets GSE183716 CNP0005824 GSE166489 GSE167029 \
  --workers 12
```

---

### 8.3 合并 h5ad

推荐使用 memory 模式：

```bash
python scripts/05_concat_h5ad_on_disk.py \
  --datasets GSE183716 CNP0005824 GSE166489 GSE167029 \
  --overwrite \
  --mode memory
```

---

### 8.4 先跑 smoke test

```bash
python scripts/06_run_v0_scanpy.py \
  --input /home/data1/wenhuai/data/h5ad_merged/GSE183716.raw_merged.h5ad \
  --dataset GSE183716 \
  --threads 8 \
  --min_genes 200 \
  --max_mito 20 \
  --n_hvg 3000
```

---

### 8.5 跑核心数据集

```bash
for ds in GSE166489 CNP0005824 GSE167029; do
  python scripts/06_run_v0_scanpy.py \
    --input /home/data1/wenhuai/data/h5ad_merged/${ds}.raw_merged.h5ad \
    --dataset ${ds} \
    --threads 24 \
    --min_genes 200 \
    --max_mito 20 \
    --n_hvg 3000 \
    2>&1 | tee logs/06_run_v0_${ds}.log
done
```

---

### 8.6 一键运行

```bash
bash scripts/run_v0_core.sh
```

如果合并时遇到 sparse matrix 格式问题，就单独运行：

```bash
python scripts/05_concat_h5ad_on_disk.py \
  --datasets GSE183716 CNP0005824 GSE166489 GSE167029 \
  --overwrite \
  --mode memory
```

---

## 9. 结果怎么看

第一轮重点看：

```bash
python scripts/07_check_v0_outputs.py \
  --datasets GSE183716 GSE166489 CNP0005824 GSE167029
```

重点打开：

```text
results/v0/GSE166489/GSE166489.top10pct_Kyn_AHR_high_summary.tsv
results/v0/GSE167029/GSE167029.top10pct_Kyn_AHR_high_summary.tsv
results/v0/CNP0005824/CNP0005824.top10pct_Kyn_AHR_high_summary.tsv
results/v0/<dataset>/<dataset>.myeloid_score_pseudobulk.tsv
results/v0/<dataset>/<dataset>.score_pseudobulk_by_sample_celltype.tsv
```

核心判断问题：

```text
1. top 10% Kyn_AHR_myeloid_score 细胞是否主要是 Monocyte / CD16_Monocyte / DC / pDC？
2. GSE166489 中 MIS-C severe 或 acute 样本是否富集 Kyn-AHR-high myeloid cells？
3. GSE167029 中 MIS-C_MYO 或 patient 样本是否富集 Kyn-AHR-high myeloid cells？
4. CNP0005824 中 CVB3 或 IVIG 相关组别是否出现髓系炎症状态变化？
5. myeloid DE 中是否出现 S100A8 / S100A9 / IL1B / TNF / CXCL8 / CCL2 / STAT1 / IRF1？
```

---

## 10. 当前版本边界

当前版本已经完成：

```text
本地数据发现
前缀式 10x 标准化
GEO soft metadata 解析
per-sample h5ad 转换
dataset-level h5ad 合并
QC
标准化
HVG
PCA
Harmony
Leiden
UMAP
marker-based PBMC 粗注释
Kyn-AHR score
myeloid inflammation score
chemotaxis score
antigen presentation score
Kyn_AHR_myeloid_score
top 10% Kyn-AHR-high cells
pseudobulk score
快速 Wilcoxon DEG
```

暂未完成：

```text
正式 pySCENIC AHR regulon
decoupleR / DoRothEA regulon activity
CellChat / NicheNet 通讯分析
髓系细胞二次聚类和精细注释
正式 sample-level pseudobulk DESeq2 / edgeR
V1 scVI-Pathway-MIL
V2 Geneformer / scGPT encoder
STT0000127 空间 GEM 数据分析
```

---

## 11. 汇总

本项目当前阶段是第二层 V0 单细胞 baseline。

核心输入：

```text
GSE183716
CNP0005824
GSE166489
GSE167029
```

核心中间文件：

```text
data/metadata/local_data_inventory.tsv
data/metadata/geo_sample_table.tsv
data/metadata/sample_sheet.auto.tsv
data/h5ad_per_sample/<dataset>/*.h5ad
data/h5ad_merged/<dataset>.raw_merged.h5ad
```

核心输出：

```text
results/v0/<dataset>/<dataset>.final_v0.h5ad
results/v0/<dataset>/<dataset>.score_pseudobulk_by_sample_celltype.tsv
results/v0/<dataset>/<dataset>.myeloid_score_pseudobulk.tsv
results/v0/<dataset>/<dataset>.top10pct_Kyn_AHR_high_summary.tsv
results/v0/<dataset>/<dataset>.DE_myeloid_disease_vs_control.wilcoxon.tsv
```

最重要的结果判断依赖：

```text
Kyn_AHR_myeloid_score
AHR_response_score
AHR_regulon_proxy_score
myeloid_inflammation_score
cell_type_v0
group
severity
stage
```

如果 V0 结果成立，下一步进入：

```text
V1: scVI encoder + pathway-guided attention MIL
V2: Geneformer / scGPT encoder + MIL
第三层: Kyn-AHR-high myeloid signature 空间转录组迁移
```

