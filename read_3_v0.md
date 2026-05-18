
# STT0000127 空间转录组分析 README

本 README 记录 STT0000127 空间转录组数据的目录结构、脚本结构、每个文件的作用、输出结果、以及图怎么看。

当前空间分析目标是：

```text
用 STT0000127 小鼠 CVB3 病毒性心肌炎空间转录组，验证心肌组织中是否存在局部 inflammatory myeloid / AHR-responsive spatial hotspots，并观察 IVIG 处理是否降低这些空间热点。
````

当前阶段的核心结论是：

```text
STT0000127 空间层面支持“局部髓系炎症热点增强”，尤其在 CVB3d6 等阶段明显；
但不支持完整 KYN-AHR-myeloid composite 或 KYN metabolism 模块在 CVB3 中整体一致升高。
```

---

# 1. 项目目录

工作目录：

```bash
/home/data1/wenhuai
```

主要目录：

```text
/home/data1/wenhuai
├── data/
│   ├── STT0000127/
│   └── spatial_h5ad/
│
├── scripts/
│   ├── spatial_00_inventory_stt.py
│   ├── spatial_01_profile_gem.py
│   ├── spatial_02_make_sample_design.py
│   ├── spatial_03_gem_to_binned_h5ad.py
│   ├── spatial_04_score_spatial_bins.py
│   ├── spatial_05_plot_score_maps.py
│   └── spatial_06_compare_spatial_scores.py
│
├── src/
│   └── kyn_ahr_v0/
│       └── gene_sets.py
│
├── results/
│   ├── spatial_inventory/
│   ├── spatial_maps/
│   └── spatial_summary/
│
└── logs/
```

---

# 2. 原始数据结构

STT0000127 原始数据在：

```bash
data/STT0000127
```

当前识别到的数据结构：

```text
data/STT0000127/
├── Analysis/
│   ├── STSA0000921/
│   │   └── STTS0001493/
│   │       ├── Control.gem.gz
│   │       ├── Control.bin50.gem.gz
│   │       └── Final.Control.TissueCut.gem.gz
│   │
│   ├── STSA0000922/
│   │   └── STTS0001494/
│   │       ├── CVB3d1.gem.gz
│   │       ├── CVB3d1.bin50.gem.gz
│   │       └── Final.CVB3d1.TissueCut.gem.gz
│   │
│   ├── ...
│   │
│   └── STSA0000935/
│       └── STTS0001507/
│           ├── IVIGd7.gem.gz
│           ├── IVIGd7.bin50.gem.gz
│           └── Final.IVIGd7.TissueCut.gem.gz
│
└── supp/
```

每个样本有三类 GEM 文件：

| 文件类型          | 示例                              | 当前用途             |
| ------------- | ------------------------------- | ---------------- |
| raw_gem       | `CVB3d6.gem.gz`                 | 原始 GEM，暂作备用      |
| bin50_gem     | `CVB3d6.bin50.gem.gz`           | 可用于快速浏览，但当前主流程不用 |
| tissuecut_gem | `Final.CVB3d6.TissueCut.gem.gz` | 主分析输入            |

当前主分析使用：

```text
Final.*.TissueCut.gem.gz
```

原因：

```text
TissueCut 文件已经做过组织区域裁剪，适合做组织内空间热点分析。
```

---

# 3. 样本设计

空间样本设计表在：

```bash
results/spatial_inventory/STT0000127_spatial_sample_design.tsv
```

样本一共 15 个：

| group            | samples                                                |
| ---------------- | ------------------------------------------------------ |
| control          | Control                                                |
| CVB3_myocarditis | CVB3d1, CVB3d2, CVB3d3, CVB3d4, CVB3d5, CVB3d6, CVB3d7 |
| IVIG_treated     | IVIGd1, IVIGd2, IVIGd3, IVIGd4, IVIGd5, IVIGd6, IVIGd7 |

字段说明：

| 字段                      | 含义                                           |
| ----------------------- | -------------------------------------------- |
| sample_id               | 样本名，例如 `CVB3d6`                              |
| group                   | 分组：control / CVB3_myocarditis / IVIG_treated |
| stage                   | 时间点，例如 day6                                  |
| stage_order             | 数字时间顺序，control 为 0，day1-day7 为 1-7           |
| preferred_role          | 当前使用的文件类型，目前是 tissuecut_gem                  |
| preferred_path          | 对应 `Final.*.TissueCut.gem.gz` 文件路径           |
| x_min/x_max/y_min/y_max | 抽样得到的空间坐标范围                                  |
| unique_genes_scanned    | 抽样统计到的基因数                                    |
| total_counts_scanned    | 抽样统计到的 MIDCount 总和                           |

查看：

```bash
column -t -s $'\t' results/spatial_inventory/STT0000127_spatial_sample_design.tsv | less -S
```

---

# 4. 脚本结构与作用

## 4.1 `scripts/spatial_00_inventory_stt.py`

作用：

```text
扫描 STT0000127 目录，识别所有 .gem.gz 文件，并判断它们是 raw_gem、bin50_gem 还是 tissuecut_gem。
```

输入：

```text
data/STT0000127/
```

输出：

```text
results/spatial_inventory/STT0000127_gem_inventory.tsv
results/spatial_inventory/STT0000127_sample_file_matrix.tsv
results/spatial_inventory/STT0000127_preferred_gem_files.tsv
results/spatial_inventory/peek_<sample>_<file_role>.txt
```

输出文件说明：

| 文件                                   | 作用                               |
| ------------------------------------ | -------------------------------- |
| `STT0000127_gem_inventory.tsv`       | 所有 GEM 文件清单，包含文件类型、大小、路径         |
| `STT0000127_sample_file_matrix.tsv`  | 每个样本对应 raw/bin50/tissuecut 的文件矩阵 |
| `STT0000127_preferred_gem_files.tsv` | 每个样本推荐使用的 GEM 文件，目前优先 tissuecut  |
| `peek_*.txt`                         | 每个 GEM 文件前若干行，用于检查列名和格式          |

运行：

```bash
python scripts/spatial_00_inventory_stt.py \
  2>&1 | tee logs/spatial_00_inventory_stt.log
```

怎么看：

```bash
column -t -s $'\t' results/spatial_inventory/STT0000127_gem_inventory.tsv | less -S
column -t -s $'\t' results/spatial_inventory/STT0000127_sample_file_matrix.tsv | less -S
column -t -s $'\t' results/spatial_inventory/STT0000127_preferred_gem_files.tsv | less -S
```

重点看：

```text
1. 是否识别到 45 个 GEM 文件。
2. 每个样本是否都有 tissuecut_gem。
3. preferred_role 是否都是 tissuecut_gem。
```

---

## 4.2 `scripts/spatial_01_profile_gem.py`

作用：

```text
对每个 preferred tissuecut GEM 做轻量 profile，不完整读入内存。
统计列名、坐标范围、抽样记录数、基因数、count 总量等。
```

输入：

```text
results/spatial_inventory/STT0000127_preferred_gem_files.tsv
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

输出字段说明：

| 字段                      | 含义                |
| ----------------------- | ----------------- |
| sample_id               | 样本名               |
| group                   | 分组                |
| stage                   | 时间点               |
| preferred_role          | 当前用 tissuecut_gem |
| path                    | GEM 文件路径          |
| n_records_scanned       | 抽样扫描的记录数          |
| n_bad_lines             | 解析失败的行数           |
| unique_genes_scanned    | 抽样中出现的唯一基因数       |
| x_min/x_max/y_min/y_max | 坐标范围              |
| total_counts_scanned    | 抽样 count 总量       |
| header                  | GEM 文件列名          |

怎么看：

```bash
column -t -s $'\t' results/spatial_inventory/STT0000127_preferred_gem_profile.tsv | less -S
```

重点看：

```text
1. n_bad_lines 是否为 0。
2. header 是否统一为 geneID|x|y|MIDCount|ExonCount。
3. 坐标范围是否合理。
```

---

## 4.3 `scripts/spatial_02_make_sample_design.py`

作用：

```text
把 profile 结果整理成正式空间样本设计表。
```

输入：

```text
results/spatial_inventory/STT0000127_preferred_gem_profile.tsv
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

重点看：

```text
Control、CVB3d1-d7、IVIGd1-d7 是否都存在。
preferred_path 是否不是空。
```

---

## 4.4 `scripts/spatial_03_gem_to_binned_h5ad.py`

作用：

```text
把 TissueCut GEM 文件按 bin_size=50 聚合成空间 bin-level AnnData/h5ad。
```

转换逻辑：

```text
geneID, x, y, MIDCount
→ x_bin = floor(x / 50)
→ y_bin = floor(y / 50)
→ spot_id = sample_id::x_bin_y_bin
→ 聚合每个 bin 的 gene counts
→ 输出 h5ad
```

输入：

```text
results/spatial_inventory/STT0000127_spatial_sample_design.tsv
Final.*.TissueCut.gem.gz
```

输出目录：

```text
data/spatial_h5ad/STT0000127/bin50_tissuecut/
```

主要输出：

```text
Control.bin50.tissuecut.h5ad
CVB3d1.bin50.tissuecut.h5ad
...
IVIGd7.bin50.tissuecut.h5ad

STT0000127_bin_h5ad_build_report.tsv
```

运行单样本 smoke test：

```bash
python scripts/spatial_03_gem_to_binned_h5ad.py \
  --samples Control \
  --bin_size 50 \
  --min_counts_per_bin 5 \
  --min_genes_per_bin 3 \
  2>&1 | tee logs/spatial_03_gem_to_binned_h5ad_Control.log
```

全量运行：

```bash
python scripts/spatial_03_gem_to_binned_h5ad.py \
  --bin_size 50 \
  --min_counts_per_bin 5 \
  --min_genes_per_bin 3 \
  2>&1 | tee logs/spatial_03_gem_to_binned_h5ad_all.log
```

输出 h5ad 内容：

```text
AnnData object
obs:
  spot
  sample_id
  group
  stage
  stage_order
  file_role
  x_bin
  y_bin
  x_center
  y_center
  n_counts
  n_genes

var:
  gene_symbols

uns:
  spatial_bin_info
```

怎么看：

```bash
python - <<'PY'
import scanpy as sc

f="/home/data1/wenhuai/data/spatial_h5ad/STT0000127/bin50_tissuecut/Control.bin50.tissuecut.h5ad"
adata=sc.read_h5ad(f)
print(adata)
print(adata.obs.head())
print(adata.var.head())
print(adata.obs[["n_counts","n_genes","x_center","y_center"]].describe())
PY
```

重点看：

```text
1. n_obs 是空间 bin 数。
2. n_vars 是基因数。
3. n_counts 和 n_genes 是否合理。
4. x_center/y_center 是否形成正常组织坐标范围。
```

---

## 4.5 `scripts/spatial_04_score_spatial_bins.py`

作用：

```text
对每个空间 bin 计算 KYN-AHR 相关模块分数。
```

输入：

```text
data/spatial_h5ad/STT0000127/bin50_tissuecut/*.h5ad
```

输出目录：

```text
data/spatial_h5ad/STT0000127/bin50_tissuecut_scored/
```

主要输出：

```text
Control.bin50.tissuecut.scored.h5ad
CVB3d1.bin50.tissuecut.scored.h5ad
...
IVIGd7.bin50.tissuecut.scored.h5ad

*.gene_presence.tsv
STT0000127_spatial_scoring_report.tsv
```

计算的 score：

| score                             | 含义                                         |
| --------------------------------- | ------------------------------------------ |
| `KYN_metabolism_score_expr`       | KYN metabolism 模块表达均值                      |
| `AHR_response_score_expr`         | AHR response 模块表达均值                        |
| `AHR_regulon_proxy_score_expr`    | AHR regulon proxy 表达均值                     |
| `myeloid_inflammation_score_expr` | 髓系炎症模块表达均值                                 |
| `Kyn_AHR_axis_score_expr`         | KYN + AHR 模块组合分数                           |
| `Kyn_AHR_myeloid_score_expr`      | KYN + AHR + myeloid inflammation composite |
| `*_z`                             | 样本内部 z-score 版本，更适合看空间局部相对高低               |

运行单样本：

```bash
python scripts/spatial_04_score_spatial_bins.py \
  --samples Control \
  --overwrite \
  2>&1 | tee logs/spatial_04_score_spatial_bins_Control.log
```

全量运行：

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

重点看：

```text
1. n_present 是否太低。
2. KYN metabolism 如果 n_present 不低但 score 很弱，说明不是注释问题，而是空间表达本身弱。
```

---

## 4.6 `scripts/spatial_05_plot_score_maps.py`

作用：

```text
绘制空间 score map。
这是看图的核心脚本。
```

这版脚本支持：

```text
1. 同一个 score 在多个样本之间统一色条范围。
2. 输出单样本图。
3. 输出 1×N panel 拼图。
4. panel 共享 colorbar。
5. panel 可共享 x/y 坐标范围。
```

输入：

```text
data/spatial_h5ad/STT0000127/bin50_tissuecut_scored/*.scored.h5ad
```

输出目录：

```text
results/spatial_maps/STT0000127/bin50_tissuecut/
```

主要输出：

```text
global_color_ranges.tsv

Control/*.png
CVB3d6/*.png
IVIGd6/*.png

_panels/Control_CVB3d6_IVIGd6/*.panel.png
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

输出 panel：

```bash
find results/spatial_maps/STT0000127/bin50_tissuecut/_panels/Control_CVB3d6_IVIGd6 \
  -name "*.png" | sort
```

## 图怎么看

以这张图为例：

```text
Control | CVB3d6 | IVIGd6
score = myeloid_inflammation_score_expr
```

解读方式：

```text
1. 三张图必须使用同一个 colorbar，否则不能比较。
2. 颜色越黄，说明该空间 bin 的 score 越高。
3. 紫色区域表示低表达或低 score。
4. 局部黄色区域表示 inflammatory hotspot。
5. 如果 CVB3d6 中黄色区域更多、更大片，而 IVIGd6 中减少，说明 IVIG 可能降低炎症空间热点负荷。
```

当前图中可见：

```text
Control 整体低。
CVB3d6 出现大面积高 myeloid inflammation score 区域。
IVIGd6 仍有热点，但整体弱于 CVB3d6。
```

因此该图支持：

```text
CVB3 myocarditis 组织中存在局部髓系炎症空间热点，IVIG 处理后热点负荷下降。
```

但不能用这张图证明：

```text
KYN metabolism 整体升高。
完整 KYN-AHR-myeloid composite 整体升高。
```

因为这些需要看对应 score 的图和统计表。

## 关于统一颜色强度

统一颜色范围记录在：

```text
results/spatial_maps/STT0000127/bin50_tissuecut/global_color_ranges.tsv
```

查看：

```bash
column -t -s $'\t' results/spatial_maps/STT0000127/bin50_tissuecut/global_color_ranges.tsv
```

字段：

| 字段          | 含义               |
| ----------- | ---------------- |
| score       | 当前 score         |
| n_values    | 用来计算色条范围的 bin 总数 |
| vmin        | 统一色条最小值          |
| vmax        | 统一色条最大值          |
| clip_q_low  | 低端裁剪分位数          |
| clip_q_high | 高端裁剪分位数          |

默认使用：

```text
q01 ~ q99
```

原因：

```text
避免极少数极端点拉爆色条，导致主体组织区域全黑。
```

如果想保留完整范围：

```bash
--clip_q_low 0.00 --clip_q_high 1.00
```

如果想突出主体区域：

```bash
--clip_q_low 0.02 --clip_q_high 0.98
```

---

## 4.7 `scripts/spatial_06_compare_spatial_scores.py`

作用：

```text
对空间 score 做样本级统计、热点比例统计、CVB3 vs IVIG day-matched 比较。
```

输入：

```text
data/spatial_h5ad/STT0000127/bin50_tissuecut_scored/*.scored.h5ad
```

输出目录：

```text
results/spatial_summary/STT0000127/
```

主要输出：

```text
STT0000127_spatial_score_sample_summary.tsv
STT0000127_control_hotspot_thresholds.tsv
STT0000127_spatial_hotspot_fraction.tsv
STT0000127_CVB3_vs_IVIG_day_matched_delta.tsv
STT0000127_CVB3_vs_IVIG_paired_wilcoxon.tsv
```

运行：

```bash
python scripts/spatial_06_compare_spatial_scores.py \
  2>&1 | tee logs/spatial_06_compare_spatial_scores.log
```

---

# 5. 统计输出怎么看

## 5.1 `STT0000127_spatial_score_sample_summary.tsv`

位置：

```text
results/spatial_summary/STT0000127/STT0000127_spatial_score_sample_summary.tsv
```

作用：

```text
每个样本、每个 score 的全组织空间分布摘要。
```

字段：

| 字段          | 含义                 |
| ----------- | ------------------ |
| sample_id   | 样本                 |
| group       | 分组                 |
| stage       | 时间点                |
| score       | 分数名                |
| n_bins      | 空间 bin 数           |
| mean        | 所有 bin 的均值         |
| median      | 所有 bin 的中位数        |
| q90/q95/q99 | 高分位数，反映热点强度        |
| top5_mean   | top 5% 高分 bin 的平均值 |
| top1_mean   | top 1% 高分 bin 的平均值 |
| max         | 最大值                |

怎么看：

```bash
grep "myeloid_inflammation_score_expr" \
  results/spatial_summary/STT0000127/STT0000127_spatial_score_sample_summary.tsv \
  | column -t -s $'\t'
```

重点：

```text
如果 median 不明显，但 q95/q99/top5_mean 明显升高，说明信号是局灶性 hotspot，而不是全组织均匀升高。
```

---

## 5.2 `STT0000127_control_hotspot_thresholds.tsv`

位置：

```text
results/spatial_summary/STT0000127/STT0000127_control_hotspot_thresholds.tsv
```

作用：

```text
用 Control 样本的 q90/q95/q99 作为 hotspot 阈值。
```

字段：

| 字段                    | 含义               |
| --------------------- | ---------------- |
| score                 | 分数名              |
| control_q90_threshold | Control 中 q90 阈值 |
| control_q95_threshold | Control 中 q95 阈值 |
| control_q99_threshold | Control 中 q99 阈值 |

怎么看：

```bash
column -t -s $'\t' results/spatial_summary/STT0000127/STT0000127_control_hotspot_thresholds.tsv
```

用途：

```text
判断 CVB3 或 IVIG 中有多少空间 bin 超过正常 Control 的高分阈值。
```

---

## 5.3 `STT0000127_spatial_hotspot_fraction.tsv`

位置：

```text
results/spatial_summary/STT0000127/STT0000127_spatial_hotspot_fraction.tsv
```

作用：

```text
计算每个样本中超过 Control q90/q95/q99 阈值的空间 bin 比例。
```

字段：

| 字段                 | 含义                     |
| ------------------ | ---------------------- |
| sample_id          | 样本                     |
| group              | 分组                     |
| stage              | 时间点                    |
| score              | 分数                     |
| threshold_level    | 使用 control_q90/q95/q99 |
| threshold          | 阈值                     |
| n_bins             | 总 bin 数                |
| n_hotspot_bins     | 超过阈值的 bin 数            |
| hotspot_fraction   | 热点面积比例                 |
| hotspot_mean_score | 热点区域平均 score           |

怎么看：

```bash
grep "myeloid_inflammation_score_expr" \
  results/spatial_summary/STT0000127/STT0000127_spatial_hotspot_fraction.tsv \
  | column -t -s $'\t'
```

关键例子：

```text
CVB3d6 control_q95 hotspot_fraction = 0.704
IVIGd6 control_q95 hotspot_fraction = 0.236
```

解释：

```text
CVB3d6 中约 70.4% 的空间 bin 超过 Control q95 炎症阈值；
IVIGd6 中约 23.6% 的空间 bin 超过该阈值。
这说明 CVB3d6 的组织炎症热点面积明显高于 IVIGd6。
```

---

## 5.4 `STT0000127_CVB3_vs_IVIG_day_matched_delta.tsv`

位置：

```text
results/spatial_summary/STT0000127/STT0000127_CVB3_vs_IVIG_day_matched_delta.tsv
```

作用：

```text
按 day1-day7 配对比较 CVB3dX 和 IVIGdX。
```

字段：

| 字段                    | 含义                                                     |
| --------------------- | ------------------------------------------------------ |
| metric                | 比较指标，如 mean/q95/top5_mean/hotspot_fraction_control_q95 |
| stage_order           | day 编号                                                 |
| score                 | score 名                                                |
| CVB3                  | CVB3dX 的值                                              |
| IVIG                  | IVIGdX 的值                                              |
| delta_CVB3_minus_IVIG | CVB3 - IVIG                                            |

怎么看：

```bash
grep "myeloid_inflammation_score_expr" \
  results/spatial_summary/STT0000127/STT0000127_CVB3_vs_IVIG_day_matched_delta.tsv \
  | column -t -s $'\t'
```

重点看：

```text
delta_CVB3_minus_IVIG > 0
```

表示：

```text
同一天 CVB3 高于 IVIG。
```

当前结果中，myeloid inflammation 多数指标在 6/7 天为正，说明 CVB3 相比 IVIG 有更强炎症热点趋势。

---

## 5.5 `STT0000127_CVB3_vs_IVIG_paired_wilcoxon.tsv`

位置：

```text
results/spatial_summary/STT0000127/STT0000127_CVB3_vs_IVIG_paired_wilcoxon.tsv
```

作用：

```text
把 day1-day7 当成配对点，对 CVB3-IVIG delta 做 Wilcoxon 检验。
```

字段：

| 字段                           | 含义              |
| ---------------------------- | --------------- |
| metric                       | 比较指标            |
| score                        | 分数              |
| n_days                       | 配对天数，通常是 7      |
| mean_delta_CVB3_minus_IVIG   | 7 天 delta 的平均值  |
| median_delta_CVB3_minus_IVIG | 7 天 delta 的中位数  |
| n_positive_days              | CVB3 > IVIG 的天数 |
| wilcoxon_p_value             | Wilcoxon p 值    |

怎么看：

```bash
column -t -s $'\t' \
  results/spatial_summary/STT0000127/STT0000127_CVB3_vs_IVIG_paired_wilcoxon.tsv \
  | less -S
```

重点指标：

```text
myeloid_inflammation_score_expr
AHR_response_score_expr
AHR_regulon_proxy_score_expr
KYN_metabolism_score_expr
Kyn_AHR_myeloid_score_expr
```

当前结论：

```text
myeloid_inflammation_score_expr：
  多数 hotspot 和高分位指标中 CVB3 > IVIG，通常 6/7 天为正，但 p 值约 0.109-0.156。
  表示方向稳定，但样本量有限，统计功效不足。

AHR_response_score_expr：
  q90/q95/top5/top1 等高分位指标较强，提示局部 AHR-responsive hotspot。

KYN_metabolism_score_expr：
  基本无稳定正向趋势，不能作为空间阳性证据。

Kyn_AHR_myeloid_score_expr：
  composite 不稳定，有些指标甚至负向，不适合作为空间主结果。
```

---

# 6. 推荐图和推荐结论

## 6.1 主图推荐

建议作为主图或主结果图：

```text
Control / CVB3d6 / IVIGd6
score = myeloid_inflammation_score_expr
```

命令：

```bash
python scripts/spatial_05_plot_score_maps.py \
  --samples Control CVB3d6 IVIGd6 \
  --scores myeloid_inflammation_score_expr \
  --clip_q_low 0.01 \
  --clip_q_high 0.99 \
  --plot_panel \
  --panel_name Control_CVB3d6_IVIGd6 \
  --share_axis_limits
```

解释：

```text
Control 低；
CVB3d6 出现大范围高炎症空间热点；
IVIGd6 仍有热点但明显弱于 CVB3d6。
```

## 6.2 辅助图推荐

建议作为辅助图：

```text
AHR_response_score_expr
AHR_regulon_proxy_score_expr
KYN_metabolism_score_expr
Kyn_AHR_myeloid_score_expr
```

其中：

```text
AHR_response_score_expr：
  支持局部 AHR-responsive signal。

AHR_regulon_proxy_score_expr：
  支持极端高分位区域有 AHR-related proxy signal。

KYN_metabolism_score_expr：
  显示空间层面 KYN metabolism 不强，是边界结果。

Kyn_AHR_myeloid_score_expr：
  composite 不稳定，不作为主图核心结论。
```

---

# 7. 推荐写法

## 7.1 中文结果描述

```text
在 STT0000127 空间转录组中，我们首先将 TissueCut GEM 文件按 50 × 50 空间 bin 聚合为 bin-level AnnData，并计算 KYN metabolism、AHR response、AHR regulon proxy 和 myeloid inflammation 等模块分数。空间图显示，CVB3d6 心肌组织中出现广泛的 myeloid inflammation 高分区域，而 IVIGd6 中该信号明显减弱。以 Control 样本 q95 作为 hotspot 阈值时，CVB3d6 中约 70.4% 的空间 bin 超过阈值，而 IVIGd6 中约 23.6% 超过阈值。Day-matched 分析显示，myeloid inflammation hotspot fraction 在多数时间点 CVB3 高于 IVIG，提示 IVIG 处理可能降低 CVB3 心肌炎中的空间炎症热点负荷。

相比之下，KYN metabolism 模块在空间 bin 层面未显示一致升高，完整 Kyn-AHR-myeloid composite 也未表现为稳定的 CVB3 特异性增强。因此，空间转录组结果更支持局部 inflammatory myeloid / AHR-responsive hotspots，而不是完整 KYN metabolism-AHR axis 的整体空间激活。
```

## 7.2 英文结果描述

```text
We converted the TissueCut GEM files from STT0000127 into 50 × 50 bin-level spatial AnnData objects and computed module-level scores for KYN metabolism, AHR response, AHR regulon proxy, and myeloid inflammation. Spatial score maps revealed extensive myeloid inflammation-high regions in CVB3d6 myocardium, whereas this signal was attenuated in the IVIGd6 sample. Using the 95th percentile of the control sample as a hotspot threshold, approximately 70.4% of spatial bins in CVB3d6 exceeded the threshold, compared with 23.6% in IVIGd6. Day-matched comparisons further showed that myeloid inflammation hotspot fractions were higher in CVB3 than IVIG in most time points, suggesting that IVIG treatment may reduce the spatial inflammatory hotspot burden in CVB3 myocarditis.

In contrast, the KYN metabolism module showed limited spatial induction, and the full KYN-AHR-myeloid composite score did not exhibit a consistent CVB3-specific increase. These findings indicate that the spatial transcriptomic data support localized inflammatory myeloid / AHR-responsive hotspots rather than a uniform activation of the complete KYN metabolism-AHR axis.
```

---

# 8. 一键复现空间分析

从已有 STT0000127 原始数据开始：

```bash
cd /home/data1/wenhuai
conda activate kynahr_v0

# 1. 文件扫描
python scripts/spatial_00_inventory_stt.py \
  2>&1 | tee logs/spatial_00_inventory_stt.log

# 2. GEM profile
python scripts/spatial_01_profile_gem.py \
  --max_records 2000000 \
  2>&1 | tee logs/spatial_01_profile_gem.log

# 3. 样本设计表
python scripts/spatial_02_make_sample_design.py \
  2>&1 | tee logs/spatial_02_make_sample_design.log

# 4. GEM 转 bin-level h5ad
python scripts/spatial_03_gem_to_binned_h5ad.py \
  --bin_size 50 \
  --min_counts_per_bin 5 \
  --min_genes_per_bin 3 \
  2>&1 | tee logs/spatial_03_gem_to_binned_h5ad_all.log

# 5. 空间 bin 打分
python scripts/spatial_04_score_spatial_bins.py \
  --overwrite \
  2>&1 | tee logs/spatial_04_score_spatial_bins_all.log

# 6. 空间图，统一色条 + panel
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

# 7. 空间统计比较
python scripts/spatial_06_compare_spatial_scores.py \
  2>&1 | tee logs/spatial_06_compare_spatial_scores.log
```

---

# 9. 当前分析边界

当前空间分析支持：

```text
1. TissueCut GEM 可以稳定转为 bin50 h5ad。
2. CVB3 心肌炎组织中存在局部 myeloid inflammation hotspots。
3. IVIG 处理组整体 hotspot burden 降低。
4. AHR response / AHR regulon proxy 在高分位或热点区域有一定趋势。
```

当前空间分析不支持：

```text
1. KYN metabolism 模块在 CVB3 中整体升高。
2. 完整 KYN-AHR-myeloid composite 在 CVB3 中稳定升高。
3. 仅凭 STT0000127 证明完整 KYN metabolism-AHR axis 被空间激活。
```

因此，空间层结论应该写成：

```text
空间转录组定位了 CVB3 心肌炎中的局部 inflammatory myeloid / AHR-responsive hotspots，并显示 IVIG 可能降低这些热点负荷；但 KYN metabolism 模块本身未显示一致空间诱导。
```

MD

````

保存后查看：

```bash
less -S README_spatial_STT0000127.md
````

这版 README 已经按你现在的真实流程写好了，尤其把**统一色条怎么看图**和**哪些结果能写、哪些不能写**单独说明了。
