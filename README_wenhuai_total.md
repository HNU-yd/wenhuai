# 文怀项目总 README：MR + 单细胞 + 空间转录组

本项目围绕 **KYN/AHR-associated inflammatory myeloid state** 在心肌炎中的多层证据开展分析。当前三步已经走通：

```text
第一步：MR 遗传因果推断
第二步：公开单细胞测序 / PBMC 与小鼠心肌炎模型分析
第三步：STT0000127 空间转录组分析
```

当前最稳妥的总结论不是“三层全部一致阳性”，而是：

```text
MR 提供候选代谢轴线索；
小鼠 CVB3 单细胞数据 CNP0005824 支持 KYN/AHR-associated inflammatory myeloid state；
人类 PBMC 数据 GSE166489 / GSE167029 未稳定复现同方向结果，作为 boundary / heterogeneity analysis；
STT0000127 空间转录组定位到 CVB3 心肌炎组织中的局部 myeloid inflammation / AHR-responsive hotspots，但 KYN metabolism 模块本身空间信号较弱。
```

---

# 1. 工作目录

项目工作目录：

```bash
/home/data1/wenhuai
```

建议进入项目后再运行所有命令：

```bash
cd /home/data1/wenhuai
conda activate kynahr_v0
```

推荐环境变量：

```bash
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8
```

如果是逐样本转换，可以把 BLAS 线程设为 1，然后用 worker 并行：

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

---

# 2. 总体数据结构

## 2.1 原始数据目录

```text
/home/data1/wenhuai/data
├── CNP0005824/
│   ├── Single_Cell/
│   │   └── CSE0000419/
│   │       ├── Control_matrix.mtx.gz
│   │       ├── Control_barcodes.tsv.gz
│   │       ├── Control_features.tsv.gz
│   │       ├── CVB3d1_matrix.mtx.gz
│   │       ├── ...
│   │       └── IVIGd7_features.tsv.gz
│   └── tenx_standard/
│       ├── Control/
│       │   ├── matrix.mtx.gz
│       │   ├── barcodes.tsv.gz
│       │   └── features.tsv.gz
│       ├── CVB3d1/
│       └── ...
│
├── GSE166489/
│   ├── matrix/
│   ├── miniml/
│   ├── soft/
│   ├── suppl/
│   └── tenx_standard/ or suppl_unpacked/
│
├── GSE167029/
│   ├── matrix/
│   ├── miniml/
│   ├── soft/
│   └── suppl/
│       ├── *.h5
│       └── GSE167029_SeuratObject_SC_AllSamples.rds.gz
│
├── GSE183716/
│   ├── matrix/
│   ├── miniml/
│   ├── soft/
│   └── suppl/
│       ├── GSE183716_Sample1_GEXFB_filtered_feature_bc_matrix.h5
│       ├── GSE183716_Sample2_GEXFB_filtered_feature_bc_matrix.h5
│       ├── GSE183716_Sample3_GEXFB_filtered_feature_bc_matrix.h5
│       └── GSE183716_Sample4_GEXFB_filtered_feature_bc_matrix.h5
│
├── GSE180045/
│   ├── matrix/
│   ├── miniml/
│   ├── soft/
│   └── suppl/
│
└── STT0000127/
    ├── Analysis/
    │   ├── STSA0000921/STTS0001493/
    │   │   ├── Control.gem.gz
    │   │   ├── Control.bin50.gem.gz
    │   │   └── Final.Control.TissueCut.gem.gz
    │   ├── STSA0000922/STTS0001494/
    │   │   ├── CVB3d1.gem.gz
    │   │   ├── CVB3d1.bin50.gem.gz
    │   │   └── Final.CVB3d1.TissueCut.gem.gz
    │   └── ...
    └── supp/
```

## 2.2 中间结果目录

```text
data/
├── metadata/
│   ├── local_data_inventory.tsv
│   ├── geo_sample_table.tsv
│   ├── sample_sheet.auto.tsv
│   └── prepare_h5ad_from_inventory_report.tsv
│
├── h5ad_per_sample/
│   ├── CNP0005824/*.h5ad
│   ├── GSE166489/*.h5ad
│   ├── GSE167029/*.h5ad
│   └── GSE183716/*.h5ad
│
├── h5ad_merged/
│   ├── CNP0005824.raw_merged.h5ad
│   ├── GSE166489.raw_merged.h5ad
│   ├── GSE167029.raw_merged.h5ad
│   └── GSE183716.raw_merged.h5ad
│
└── spatial_h5ad/
    └── STT0000127/
        ├── bin50_tissuecut/
        │   ├── Control.bin50.tissuecut.h5ad
        │   ├── CVB3d1.bin50.tissuecut.h5ad
        │   └── ...
        └── bin50_tissuecut_scored/
            ├── Control.bin50.tissuecut.scored.h5ad
            ├── CVB3d1.bin50.tissuecut.scored.h5ad
            └── ...
```

## 2.3 结果目录

```text
results/
├── v0/
│   ├── CNP0005824/
│   ├── GSE166489/
│   ├── GSE167029/
│   └── GSE183716/
│
├── v0_summary/
│   ├── v0_all_myeloid_pseudobulk_merged.tsv
│   ├── v0_myeloid_group_means.tsv
│   ├── v0_pairwise_tests.tsv
│   ├── v0_boundary_flags.tsv
│   ├── GSE166489_paired_acute_recovery_detail.tsv
│   └── GSE166489_paired_acute_recovery_summary.tsv
│
├── v0_figures/
│   ├── Fig_scRNA_1_boundary_effect_heatmap.png
│   ├── Fig_scRNA_2_CNP0005824_Kyn_AHR_myeloid_score.png
│   ├── Fig_scRNA_2b_CNP0005824_myeloid_inflammation_score.png
│   ├── Fig_scRNA_3_human_boundary_Kyn_AHR_myeloid_score.png
│   └── Fig_scRNA_3b_human_boundary_myeloid_inflammation_score.png
│
├── spatial_inventory/
│   ├── STT0000127_gem_inventory.tsv
│   ├── STT0000127_sample_file_matrix.tsv
│   ├── STT0000127_preferred_gem_files.tsv
│   ├── STT0000127_preferred_gem_profile.tsv
│   └── STT0000127_spatial_sample_design.tsv
│
├── spatial_maps/
│   └── STT0000127/bin50_tissuecut/
│       ├── global_color_ranges.tsv
│       ├── Control/*.png
│       ├── CVB3d6/*.png
│       ├── IVIGd6/*.png
│       └── _panels/Control_CVB3d6_IVIGd6/*.panel.png
│
└── spatial_summary/
    └── STT0000127/
        ├── STT0000127_spatial_score_sample_summary.tsv
        ├── STT0000127_control_hotspot_thresholds.tsv
        ├── STT0000127_spatial_hotspot_fraction.tsv
        ├── STT0000127_CVB3_vs_IVIG_day_matched_delta.tsv
        └── STT0000127_CVB3_vs_IVIG_paired_wilcoxon.tsv
```

---

# 3. 脚本文件结构

```text
scripts/
├── 01_discover_local_data.py
├── 02_standardize_prefixed_10x_copy.py
├── 03_build_geo_sample_table.py
├── 04_prepare_h5ad_from_inventory.py
├── 05_concat_h5ad_on_disk.py
├── 06_run_v0_scanpy.py
├── 07_check_v0_outputs.py
├── 08_v0_score_boundary_analysis.py
├── 09_plot_scrna_v0_overview.py
├── run_v0_core.sh
│
├── spatial_00_inventory_stt.py
├── spatial_01_profile_gem.py
├── spatial_02_make_sample_design.py
├── spatial_03_gem_to_binned_h5ad.py
├── spatial_04_score_spatial_bins.py
├── spatial_05_plot_score_maps.py
└── spatial_06_compare_spatial_scores.py

src/
└── kyn_ahr_v0/
    ├── __init__.py
    ├── gene_sets.py
    └── metadata_rules.py
```

---

# 4. 第一层：MR 分析

当前 `scripts.zip` 中没有纳入 MR 脚本，因此 MR 部分在本 README 中作为结果归档规范说明。

建议 MR 结果统一放在：

```text
results/mr/
├── harmonised_data/
├── main_results/
├── sensitivity/
├── pleiotropy/
├── heterogeneity/
└── figures/
```

常见输出及用途：

| 文件类型 | 作用 | 怎么看 |
|---|---|---|
| `*_main_results.tsv` | IVW、MR-Egger、weighted median 等主结果 | 看 beta/OR、p value、方向是否一致 |
| `*_heterogeneity.tsv` | 异质性检验 | Q_pval 过低提示工具变量异质性强 |
| `*_pleiotropy.tsv` | MR-Egger intercept | p 值显著提示水平多效性风险 |
| `*_leave_one_out.tsv` | leave-one-out 稳健性 | 看是否由单个 SNP 驱动 |
| `figures/*.pdf/png` | forest/funnel/scatter/leave-one-out 图 | 用于报告第一层结果展示 |

MR 报告写法建议：

```text
MR 层主要用于提供候选代谢轴线索，而不是单独证明机制。若 KYN/tryptophan metabolism 相关暴露与心肌炎风险存在方向性结果，则作为后续单细胞和空间转录组聚焦 KYN/AHR 相关模块的依据。
```

---

# 5. 第二层：单细胞测序 V0 pipeline

第二层使用公开单细胞数据，主要分析髓系细胞中的 KYN/AHR-associated inflammatory state。

核心数据集：

| 数据集 | 类型 | 当前定位 |
|---|---|---|
| CNP0005824 | 小鼠 CVB3 心肌炎单细胞 | 最支持主线 |
| GSE166489 | 人类 MIS-C PBMC | boundary / heterogeneity analysis |
| GSE167029 | 人类 MIS-C / MIS-C_MYO PBMC | boundary / heterogeneity analysis |
| GSE183716 | 人类小规模 MIS-C 数据 | smoke / 辅助验证 |
| GSE180045 | ICI myocarditis 扩展 | 暂不进核心结论 |

## 5.1 `01_discover_local_data.py`

作用：

```text
扫描 data/ 下本地数据，识别 h5、h5ad、RDS、标准 10x mtx 目录、GEO soft、series matrix 等文件。
```

输入：

```text
data/<dataset>/...
```

输出：

```text
data/metadata/local_data_inventory.tsv
```

运行：

```bash
python scripts/01_discover_local_data.py \
  --datasets CNP0005824 GSE166489 GSE167029 GSE183716 \
  2>&1 | tee logs/01_discover_local_data.log
```

怎么看：

```bash
grep -E "h5ad|h5_or_hdf5|10x_mtx_dir|seurat_rds" \
  data/metadata/local_data_inventory.tsv \
  | column -t -s $'\t' | less -S
```

重点：

```text
CNP0005824 应该能识别出 10x_mtx_dir。
GSE166489 应该能识别出 10x_mtx_dir。
GSE167029 应该能识别出 h5_or_hdf5。
GSE183716 应该能识别出 h5_or_hdf5。
```

## 5.2 `02_standardize_prefixed_10x_copy.py`

作用：

```text
把 CNP0005824 这类前缀式 10x 文件复制成标准 10x 目录。
不使用软链接，不破坏原始文件。
```

输入示例：

```text
Control_matrix.mtx.gz
Control_barcodes.tsv.gz
Control_features.tsv.gz
```

输出示例：

```text
data/CNP0005824/tenx_standard/Control/
├── matrix.mtx.gz
├── barcodes.tsv.gz
└── features.tsv.gz
```

运行：

```bash
python scripts/02_standardize_prefixed_10x_copy.py \
  --datasets CNP0005824 GSE166489 GSE167029 \
  2>&1 | tee logs/02_standardize_prefixed_10x_copy.log
```

怎么看：

```bash
find data/CNP0005824/tenx_standard -maxdepth 2 -type f | sort | head -50
```

重点：

```text
每个样本目录必须有 matrix.mtx.gz、barcodes.tsv.gz、features.tsv.gz 三个文件。
```

## 5.3 `03_build_geo_sample_table.py`

作用：

```text
解析 GEO family.soft 文件，提取 GSM、Sample_title、source_name、characteristics，用于辅助自动 metadata 推断。
```

输出：

```text
data/metadata/geo_sample_table.tsv
```

运行：

```bash
python scripts/03_build_geo_sample_table.py \
  2>&1 | tee logs/03_build_geo_sample_table.log
```

怎么看：

```bash
column -t -s $'\t' data/metadata/geo_sample_table.tsv | less -S
```

## 5.4 `04_prepare_h5ad_from_inventory.py`

作用：

```text
根据 local_data_inventory.tsv，把 h5 / 10x mtx 目录转换成单样本 h5ad。
同时写入 sample_id、patient_id、group、severity、stage 等 metadata。
```

输入：

```text
data/metadata/local_data_inventory.tsv
data/metadata/geo_sample_table.tsv
```

输出：

```text
data/h5ad_per_sample/<dataset>/*.h5ad
data/metadata/prepare_h5ad_from_inventory_report.tsv
data/metadata/sample_sheet.auto.tsv
```

运行：

```bash
python scripts/04_prepare_h5ad_from_inventory.py \
  --datasets GSE183716 CNP0005824 GSE166489 GSE167029 \
  --workers 12 \
  2>&1 | tee logs/04_prepare_core_h5ad.log
```

怎么看：

```bash
cat data/metadata/prepare_h5ad_from_inventory_report.tsv | column -t -s $'\t' | less -S

for ds in GSE183716 CNP0005824 GSE166489 GSE167029; do
  echo "==== $ds ===="
  find data/h5ad_per_sample/$ds -name "*.h5ad" | wc -l
done
```

重点：

```text
status 应该是 ok 或 skip。
error 行需要检查。
```

## 5.5 `05_concat_h5ad_on_disk.py`

作用：

```text
把每个数据集的单样本 h5ad 合并成数据集级 h5ad。
当前版本会把 csc 稀疏矩阵转成 csr，并使用临时文件写入，避免 concat_on_disk 的 csc 不支持问题。
```

输入：

```text
data/h5ad_per_sample/<dataset>/*.h5ad
```

输出：

```text
data/h5ad_merged/<dataset>.raw_merged.h5ad
```

推荐运行：

```bash
python scripts/05_concat_h5ad_on_disk.py \
  --datasets GSE183716 CNP0005824 GSE166489 GSE167029 \
  --overwrite \
  --mode memory \
  2>&1 | tee logs/05_concat_core_h5ad.memory.log
```

怎么看：

```bash
ls -lh data/h5ad_merged/*.raw_merged.h5ad

python - <<'PY'
import scanpy as sc
for ds in ["GSE183716", "CNP0005824", "GSE166489", "GSE167029"]:
    f=f"/home/data1/wenhuai/data/h5ad_merged/{ds}.raw_merged.h5ad"
    adata=sc.read_h5ad(f, backed="r")
    print(ds, adata)
PY
```

## 5.6 `06_run_v0_scanpy.py`

作用：

```text
运行单细胞 V0 主流程：QC、normalize、log1p、模块打分、粗细胞类型注释、HVG、PCA、Harmony、UMAP、Leiden、pseudobulk、DE 和图。
```

输入：

```text
data/h5ad_merged/<dataset>.raw_merged.h5ad
```

输出目录：

```text
results/v0/<dataset>/
```

主要输出：

| 文件 | 作用 | 怎么看 |
|---|---|---|
| `<dataset>.01_qc.h5ad` | QC 后对象 | 后续排错用 |
| `<dataset>.02_v0_processed.h5ad` | 完成 UMAP/score 的对象 | 后续复用 |
| `<dataset>.final_v0.h5ad` | 最终对象 | 下游分析优先用 |
| `<dataset>.score_gene_presence.tsv` | 每个 score 命中多少基因 | 检查 gene set 是否可用 |
| `<dataset>.obs_after_scores.tsv` | 每个细胞的 metadata 和 score | 单细胞层面追踪 |
| `<dataset>.celltype_proportions.tsv` | 每个样本各细胞类型比例 | 看细胞组成 |
| `<dataset>.score_pseudobulk_by_sample_celltype.tsv` | sample × cell type 的所有 score 均值 | 主分析表之一 |
| `<dataset>.myeloid_score_pseudobulk.tsv` | 髓系细胞 pseudobulk score | 第二步核心表 |
| `<dataset>.top10pct_Kyn_AHR_high_summary.tsv` | top 10% 高分细胞来自哪些组和细胞类型 | 看高分细胞富集 |
| `<dataset>.DE_*.tsv` | 快速差异基因 | 仅作辅助 |
| `figures/*.png` | UMAP 图 | 看聚类和 score 分布 |

运行单个数据集：

```bash
python scripts/06_run_v0_scanpy.py \
  --input /home/data1/wenhuai/data/h5ad_merged/CNP0005824.raw_merged.h5ad \
  --dataset CNP0005824 \
  --threads 12 \
  --min_genes 200 \
  --max_mito 20 \
  --n_hvg 3000 \
  2>&1 | tee logs/06_run_v0_CNP0005824.log
```

核心全量：

```bash
for ds in GSE183716 GSE166489 CNP0005824 GSE167029; do
  python scripts/06_run_v0_scanpy.py \
    --input /home/data1/wenhuai/data/h5ad_merged/${ds}.raw_merged.h5ad \
    --dataset ${ds} \
    --threads 8 \
    --min_genes 200 \
    --max_mito 20 \
    --n_hvg 3000 \
    2>&1 | tee logs/06_run_v0_${ds}.log
done
```

## 5.7 `07_check_v0_outputs.py`

作用：

```text
检查 h5ad_per_sample、merged h5ad、final_v0.h5ad 和 Kyn-AHR high summary 是否存在。
```

运行：

```bash
python scripts/07_check_v0_outputs.py \
  --datasets GSE183716 GSE166489 CNP0005824 GSE167029 \
  2>&1 | tee logs/07_check_v0_outputs.log
```

## 5.8 `08_v0_score_boundary_analysis.py`

作用：

```text
整合各数据集 myeloid_score_pseudobulk.tsv，做组均值、pairwise test、boundary flags 和 GSE166489 paired acute-recovery 分析。
```

输入：

```text
results/v0/<dataset>/<dataset>.myeloid_score_pseudobulk.tsv
```

输出：

```text
results/v0_summary/v0_all_myeloid_pseudobulk_merged.tsv
results/v0_summary/v0_myeloid_group_means.tsv
results/v0_summary/v0_pairwise_tests.tsv
results/v0_summary/v0_boundary_flags.tsv
results/v0_summary/GSE166489_paired_acute_recovery_detail.tsv
results/v0_summary/GSE166489_paired_acute_recovery_summary.tsv
```

运行：

```bash
python scripts/08_v0_score_boundary_analysis.py \
  --datasets CNP0005824 GSE166489 GSE167029 GSE183716 \
  2>&1 | tee logs/08_v0_score_boundary_analysis.log
```

怎么看：

```bash
column -t -s $'\t' results/v0_summary/v0_boundary_flags.tsv | less -S
```

重点：

```text
CNP0005824: CVB3_myocarditis > control / IVIG_treated 应该多为 supports_expected_direction。
GSE166489: MIS-C > pediatric_healthy 通常不支持。
GSE167029: MIS-C_MYO > control 通常不支持。
```

## 5.9 `09_plot_scrna_v0_overview.py`

作用：

```text
生成第二步单细胞的一眼看懂图，包括 boundary heatmap、CNP 正向结果图、人类 PBMC boundary 图。
```

输入：

```text
results/v0/<dataset>/<dataset>.myeloid_score_pseudobulk.tsv
```

输出：

```text
results/v0_figures/Fig_scRNA_1_boundary_effect_heatmap.png
results/v0_figures/Fig_scRNA_2_CNP0005824_Kyn_AHR_myeloid_score.png
results/v0_figures/Fig_scRNA_2b_CNP0005824_myeloid_inflammation_score.png
results/v0_figures/Fig_scRNA_3_human_boundary_Kyn_AHR_myeloid_score.png
results/v0_figures/Fig_scRNA_3b_human_boundary_myeloid_inflammation_score.png
results/v0_figures/scrna_boundary_effect_table.tsv
results/v0_figures/scrna_myeloid_pseudobulk_merged.tsv
```

运行：

```bash
python scripts/09_plot_scrna_v0_overview.py \
  2>&1 | tee logs/09_plot_scrna_v0_overview.log
```

图怎么看：

### `Fig_scRNA_1_boundary_effect_heatmap.png`

这是第二步最重要的总览图。

```text
行 = 数据集比较 + 髓系细胞类型
列 = score 模块
颜色 = 目标组均值 - 参考组均值
红色 = 目标组更高
蓝色 = 目标组更低
```

推荐解读：

```text
CNP0005824 中 CVB3 myocarditis 相比 control 或 IVIG_treated 多数模块和髓系细胞类型呈正向升高。
GSE166489 和 GSE167029 中，人类 PBMC 的 MIS-C 或 MIS-C_MYO 未显示一致正向升高。
该图支持“小鼠 CVB3 模型强，人类 PBMC 是边界结果”的叙事。
```

### `Fig_scRNA_2_CNP0005824_Kyn_AHR_myeloid_score.png`

```text
展示 CNP0005824 中 control、IVIG_treated、CVB3_myocarditis 在不同髓系细胞中的分数分布。
如果 CVB3 高、IVIG 下降，说明小鼠模型支持主线。
```

### `Fig_scRNA_3_human_boundary_Kyn_AHR_myeloid_score.png`

```text
展示 GSE166489 和 GSE167029 中人类 PBMC 的结果。
主要作用不是证明阳性，而是展示人类 PBMC 没有稳定复现同方向升高。
```

## 5.10 `run_v0_core.sh`

作用：

```text
一键运行单细胞核心 pipeline：标准化 10x、发现文件、建 metadata、转 h5ad、合并、跑 V0、检查输出。
```

运行：

```bash
bash scripts/run_v0_core.sh
```

如果你已经有 `data/h5ad_merged/*.raw_merged.h5ad`，不建议重复全跑，可以只重跑：

```bash
for ds in GSE183716 GSE166489 CNP0005824 GSE167029; do
  python scripts/06_run_v0_scanpy.py \
    --input /home/data1/wenhuai/data/h5ad_merged/${ds}.raw_merged.h5ad \
    --dataset ${ds} \
    --threads 8

done

python scripts/08_v0_score_boundary_analysis.py \
  --datasets CNP0005824 GSE166489 GSE167029 GSE183716

python scripts/09_plot_scrna_v0_overview.py
```

---

# 6. 第三层：STT0000127 空间转录组 pipeline

空间层目标：

```text
用小鼠 CVB3 心肌炎空间转录组定位组织中的 inflammatory myeloid / AHR-responsive hotspots，并评估 IVIG 是否降低 hotspot burden。
```

当前空间层支持：

```text
CVB3d6 等阶段有明显 myeloid inflammation spatial hotspots。
IVIGd6 热点负荷下降。
AHR response / AHR regulon proxy 在高分位区域有趋势。
KYN metabolism 模块空间诱导有限。
完整 Kyn-AHR-myeloid composite 不稳定。
```

## 6.1 `spatial_00_inventory_stt.py`

作用：

```text
扫描 STT0000127 中所有 .gem.gz 文件，识别 raw_gem、bin50_gem、tissuecut_gem。
```

输出：

```text
results/spatial_inventory/STT0000127_gem_inventory.tsv
results/spatial_inventory/STT0000127_sample_file_matrix.tsv
results/spatial_inventory/STT0000127_preferred_gem_files.tsv
results/spatial_inventory/peek_<sample>_<file_role>.txt
```

运行：

```bash
python scripts/spatial_00_inventory_stt.py \
  2>&1 | tee logs/spatial_00_inventory_stt.log
```

怎么看：

```bash
column -t -s $'\t' results/spatial_inventory/STT0000127_gem_inventory.tsv | less -S
column -t -s $'\t' results/spatial_inventory/STT0000127_sample_file_matrix.tsv | less -S
```

重点：

```text
应识别到 45 个 GEM 文件。
每个样本应有 raw_gem、bin50_gem、tissuecut_gem。
preferred_role 应优先为 tissuecut_gem。
```

## 6.2 `spatial_01_profile_gem.py`

作用：

```text
对 preferred tissuecut GEM 做轻量 profile，检查列名、坐标范围、记录数、坏行数。
```

输出：

```text
results/spatial_inventory/STT0000127_preferred_gem_profile.tsv
```

运行：

```bash
python scripts/spatial_01_profile_gem.py \
  --max_records 2000000 \
  2>&1 | tee logs/spatial_01_profile_gem.log
```

怎么看：

```bash
column -t -s $'\t' results/spatial_inventory/STT0000127_preferred_gem_profile.tsv | less -S
```

重点：

```text
header 应为 geneID|x|y|MIDCount|ExonCount。
n_bad_lines 应接近 0。
```

## 6.3 `spatial_02_make_sample_design.py`

作用：

```text
生成正式空间样本设计表。
```

输出：

```text
results/spatial_inventory/STT0000127_spatial_sample_design.tsv
```

运行：

```bash
python scripts/spatial_02_make_sample_design.py \
  2>&1 | tee logs/spatial_02_make_sample_design.log
```

怎么看：

```bash
column -t -s $'\t' results/spatial_inventory/STT0000127_spatial_sample_design.tsv | less -S
```

重点：

```text
应包含 Control、CVB3d1-d7、IVIGd1-d7 共 15 个样本。
preferred_path 不能是空。
```

## 6.4 `spatial_03_gem_to_binned_h5ad.py`

作用：

```text
将 Final.*.TissueCut.gem.gz 按 50×50 空间 bin 聚合为 bin-level AnnData。
```

转换逻辑：

```text
geneID, x, y, MIDCount, ExonCount
→ x_bin = floor(x / bin_size)
→ y_bin = floor(y / bin_size)
→ spot = x_bin_y_bin
→ 聚合每个 spot 的 gene counts
→ 输出 h5ad
```

输出：

```text
data/spatial_h5ad/STT0000127/bin50_tissuecut/<sample>.bin50.tissuecut.h5ad
data/spatial_h5ad/STT0000127/bin50_tissuecut/STT0000127_bin_h5ad_build_report.tsv
```

单样本 smoke：

```bash
python scripts/spatial_03_gem_to_binned_h5ad.py \
  --samples Control \
  --bin_size 50 \
  --min_counts_per_bin 5 \
  --min_genes_per_bin 3 \
  2>&1 | tee logs/spatial_03_gem_to_binned_h5ad_Control.log
```

全量：

```bash
python scripts/spatial_03_gem_to_binned_h5ad.py \
  --bin_size 50 \
  --min_counts_per_bin 5 \
  --min_genes_per_bin 3 \
  2>&1 | tee logs/spatial_03_gem_to_binned_h5ad_all.log
```

怎么看：

```bash
column -t -s $'\t' data/spatial_h5ad/STT0000127/bin50_tissuecut/STT0000127_bin_h5ad_build_report.tsv | less -S
```

示例检查：

```bash
python - <<'PY'
import scanpy as sc
f="/home/data1/wenhuai/data/spatial_h5ad/STT0000127/bin50_tissuecut/Control.bin50.tissuecut.h5ad"
adata=sc.read_h5ad(f)
print(adata)
print(adata.obs.head())
print(adata.obs[["n_counts", "n_genes", "x_center", "y_center"]].describe())
PY
```

## 6.5 `spatial_04_score_spatial_bins.py`

作用：

```text
对每个空间 bin 计算 KYN metabolism、AHR response、AHR regulon proxy、myeloid inflammation 和 composite scores。
```

输出：

```text
data/spatial_h5ad/STT0000127/bin50_tissuecut_scored/<sample>.bin50.tissuecut.scored.h5ad
data/spatial_h5ad/STT0000127/bin50_tissuecut_scored/*.gene_presence.tsv
data/spatial_h5ad/STT0000127/bin50_tissuecut_scored/STT0000127_spatial_scoring_report.tsv
```

运行：

```bash
python scripts/spatial_04_score_spatial_bins.py \
  --overwrite \
  2>&1 | tee logs/spatial_04_score_spatial_bins_all.log
```

怎么看 gene presence：

```bash
cat data/spatial_h5ad/STT0000127/bin50_tissuecut_scored/*.gene_presence.tsv \
  > results/spatial_summary/STT0000127_gene_presence_all.tsv

grep "KYN_metabolism_score" results/spatial_summary/STT0000127_gene_presence_all.tsv
grep "AHR_response_score" results/spatial_summary/STT0000127_gene_presence_all.tsv
grep "myeloid_inflammation_score" results/spatial_summary/STT0000127_gene_presence_all.tsv
```

## 6.6 `spatial_05_plot_score_maps.py`

作用：

```text
绘制空间 score map 和论文用 panel 图。
同一个 score 在多个样本间统一色条强度，避免误导。
```

输入：

```text
data/spatial_h5ad/STT0000127/bin50_tissuecut_scored/*.scored.h5ad
```

输出：

```text
results/spatial_maps/STT0000127/bin50_tissuecut/global_color_ranges.tsv
results/spatial_maps/STT0000127/bin50_tissuecut/<sample>/*.png
results/spatial_maps/STT0000127/bin50_tissuecut/_panels/<panel_name>/*.panel.png
```

推荐主图命令：

```bash
python scripts/spatial_05_plot_score_maps.py \
  --samples Control CVB3d6 IVIGd6 \
  --scores myeloid_inflammation_score_expr AHR_response_score_expr AHR_regulon_proxy_score_expr KYN_metabolism_score_expr Kyn_AHR_myeloid_score_expr \
  --clip_q_low 0.01 \
  --clip_q_high 0.99 \
  --plot_panel \
  --plot_single \
  --panel_name Control_CVB3d6_IVIGd6 \
  --share_axis_limits \
  2>&1 | tee logs/spatial_05_plot_key_maps.panel.log
```

图怎么看：

```text
同一个 score 的所有样本必须用同一个 colorbar。
紫色表示低 score，黄色表示高 score。
局部黄色区域表示 hotspot。
Control 低、CVB3d6 高、IVIGd6 降低，说明 IVIG 可能降低空间炎症热点负荷。
```

统一色条范围记录在：

```text
results/spatial_maps/STT0000127/bin50_tissuecut/global_color_ranges.tsv
```

查看：

```bash
column -t -s $'\t' results/spatial_maps/STT0000127/bin50_tissuecut/global_color_ranges.tsv
```

关键图推荐：

| 图 | 用途 |
|---|---|
| `myeloid_inflammation_score_expr.panel.png` | 主图，展示 CVB3 炎症热点和 IVIG 缓解 |
| `AHR_response_score_expr.panel.png` | 辅助图，展示 AHR-responsive hotspot 趋势 |
| `AHR_regulon_proxy_score_expr.panel.png` | 辅助图，展示 AHR proxy hotspot |
| `KYN_metabolism_score_expr.panel.png` | boundary/negative 图，显示 KYN metabolism 空间信号不强 |
| `Kyn_AHR_myeloid_score_expr.panel.png` | composite 图，谨慎解读 |

## 6.7 `spatial_06_compare_spatial_scores.py`

作用：

```text
统计空间 score 的样本级分布、高分位、hotspot fraction，并做 CVB3 vs IVIG 的 day-matched 比较。
```

输出：

```text
results/spatial_summary/STT0000127/STT0000127_spatial_score_sample_summary.tsv
results/spatial_summary/STT0000127/STT0000127_control_hotspot_thresholds.tsv
results/spatial_summary/STT0000127/STT0000127_spatial_hotspot_fraction.tsv
results/spatial_summary/STT0000127/STT0000127_CVB3_vs_IVIG_day_matched_delta.tsv
results/spatial_summary/STT0000127/STT0000127_CVB3_vs_IVIG_paired_wilcoxon.tsv
```

运行：

```bash
python scripts/spatial_06_compare_spatial_scores.py \
  2>&1 | tee logs/spatial_06_compare_spatial_scores.log
```

怎么看：

```bash
column -t -s $'\t' results/spatial_summary/STT0000127/STT0000127_CVB3_vs_IVIG_paired_wilcoxon.tsv | less -S

grep "myeloid_inflammation_score_expr" \
  results/spatial_summary/STT0000127/STT0000127_spatial_hotspot_fraction.tsv \
  | column -t -s $'\t'

grep "myeloid_inflammation_score_expr" \
  results/spatial_summary/STT0000127/STT0000127_CVB3_vs_IVIG_day_matched_delta.tsv \
  | column -t -s $'\t'
```

重点解读：

```text
hotspot_fraction_control_q95：超过 Control q95 阈值的空间 bin 比例。
CVB3d6 hotspot_fraction 明显高于 IVIGd6，说明 CVB3d6 炎症热点面积更大。
n_positive_days 表示 7 个时间点中 CVB3 > IVIG 的天数。
```

---

# 7. 一键运行顺序

## 7.1 单细胞完整核心流程

```bash
cd /home/data1/wenhuai
conda activate kynahr_v0

bash scripts/run_v0_core.sh
```

如果已经完成 h5ad 合并，只想重跑结果和图：

```bash
for ds in GSE183716 GSE166489 CNP0005824 GSE167029; do
  python scripts/06_run_v0_scanpy.py \
    --input /home/data1/wenhuai/data/h5ad_merged/${ds}.raw_merged.h5ad \
    --dataset ${ds} \
    --threads 8 \
    --min_genes 200 \
    --max_mito 20 \
    --n_hvg 3000 \
    2>&1 | tee logs/06_run_v0_${ds}.log
done

python scripts/08_v0_score_boundary_analysis.py \
  --datasets CNP0005824 GSE166489 GSE167029 GSE183716 \
  2>&1 | tee logs/08_v0_score_boundary_analysis.log

python scripts/09_plot_scrna_v0_overview.py \
  2>&1 | tee logs/09_plot_scrna_v0_overview.log
```

## 7.2 空间转录组完整流程

```bash
cd /home/data1/wenhuai
conda activate kynahr_v0

python scripts/spatial_00_inventory_stt.py \
  2>&1 | tee logs/spatial_00_inventory_stt.log

python scripts/spatial_01_profile_gem.py \
  --max_records 2000000 \
  2>&1 | tee logs/spatial_01_profile_gem.log

python scripts/spatial_02_make_sample_design.py \
  2>&1 | tee logs/spatial_02_make_sample_design.log

python scripts/spatial_03_gem_to_binned_h5ad.py \
  --bin_size 50 \
  --min_counts_per_bin 5 \
  --min_genes_per_bin 3 \
  2>&1 | tee logs/spatial_03_gem_to_binned_h5ad_all.log

python scripts/spatial_04_score_spatial_bins.py \
  --overwrite \
  2>&1 | tee logs/spatial_04_score_spatial_bins_all.log

python scripts/spatial_05_plot_score_maps.py \
  --samples Control CVB3d6 IVIGd6 \
  --scores myeloid_inflammation_score_expr AHR_response_score_expr AHR_regulon_proxy_score_expr KYN_metabolism_score_expr Kyn_AHR_myeloid_score_expr \
  --clip_q_low 0.01 \
  --clip_q_high 0.99 \
  --plot_panel \
  --plot_single \
  --panel_name Control_CVB3d6_IVIGd6 \
  --share_axis_limits \
  2>&1 | tee logs/spatial_05_plot_key_maps.panel.log

python scripts/spatial_06_compare_spatial_scores.py \
  2>&1 | tee logs/spatial_06_compare_spatial_scores.log
```

---

# 8. 报告中推荐展示的图

## 8.1 MR 图

根据 MR 输出选择：

```text
MR forest plot
MR scatter plot
leave-one-out plot
```

用途：

```text
展示候选代谢轴与心肌炎风险的因果线索。
```

## 8.2 单细胞图

推荐主图：

```text
results/v0_figures/Fig_scRNA_1_boundary_effect_heatmap.png
```

用途：

```text
一眼看出 CNP0005824 支持主线，人类 PBMC 数据不稳定复现。
```

推荐辅助图：

```text
results/v0_figures/Fig_scRNA_2_CNP0005824_Kyn_AHR_myeloid_score.png
results/v0_figures/Fig_scRNA_3_human_boundary_Kyn_AHR_myeloid_score.png
```

## 8.3 空间图

推荐主图：

```text
results/spatial_maps/STT0000127/bin50_tissuecut/_panels/Control_CVB3d6_IVIGd6/Control_CVB3d6_IVIGd6.myeloid_inflammation_score_expr.panel.png
```

用途：

```text
显示 Control 低、CVB3d6 出现大范围炎症 hotspot、IVIGd6 hotspot 下降。
```

推荐辅助图：

```text
AHR_response_score_expr.panel.png
AHR_regulon_proxy_score_expr.panel.png
KYN_metabolism_score_expr.panel.png
Kyn_AHR_myeloid_score_expr.panel.png
```

---

# 9. 当前结果解释边界

可以写：

```text
小鼠 CVB3 单细胞数据支持 KYN/AHR-associated inflammatory myeloid state。
STT0000127 空间转录组支持 CVB3 心肌炎组织中的局部 myeloid inflammation / AHR-responsive hotspots，并显示 IVIG 后热点负荷下降。
```

不能强行写：

```text
所有人类 PBMC 数据都支持 KYN-AHR-myeloid signature 升高。
空间转录组证明 KYN metabolism 模块整体升高。
完整 KYN-AHR-myeloid composite 在空间层面稳定升高。
```

推荐最终表述：

```text
本研究通过 MR、单细胞和空间转录组构建了 KYN/AHR-associated inflammatory myeloid state 的多层证据框架。小鼠 CVB3 心肌炎单细胞数据为该状态提供了最强支持，空间转录组进一步显示 CVB3 心肌组织中存在局部 inflammatory myeloid / AHR-responsive hotspots，并在 IVIG 处理后降低。然而，人类 PBMC 数据未稳定复现同方向升高，且空间 KYN metabolism 模块本身诱导有限，提示该信号具有物种、疾病阶段和组织微环境依赖性。
```

---

# 10. 常见问题

## 10.1 为什么空间图一定要统一 colorbar？

如果每个样本单独自动缩放颜色，Control 的 0-0.5 和 CVB3d6 的 0-0.16 会看起来都很亮，导致无法比较。

现在 `spatial_05_plot_score_maps.py` 会对同一个 score 在所有指定样本上统一计算 vmin/vmax，并保存：

```text
results/spatial_maps/STT0000127/bin50_tissuecut/global_color_ranges.tsv
```

## 10.2 为什么 KYN metabolism 图看起来不强？

因为当前空间 bin-level scoring 中，KYN metabolism 模块没有显示稳定 CVB3 特异性增强。这是重要 boundary result，不应该删除。

## 10.3 为什么 composite score 不作为空间主图？

`Kyn_AHR_myeloid_score_expr` 同时混合 KYN metabolism、AHR response、AHR proxy 和 myeloid inflammation。当前 KYN metabolism 弱，会稀释或抵消 myeloid inflammation/AHR hotspot，因此 composite 不如模块拆分图清楚。

## 10.4 第二步最直观的图是哪一个？

```text
results/v0_figures/Fig_scRNA_1_boundary_effect_heatmap.png
```

这张图最适合作为单细胞总览图。

## 10.5 第三步最直观的图是哪一个？

```text
results/spatial_maps/STT0000127/bin50_tissuecut/_panels/Control_CVB3d6_IVIGd6/Control_CVB3d6_IVIGd6.myeloid_inflammation_score_expr.panel.png
```

这张图最适合作为空间转录组主图。
