suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- normalizePath(Sys.getenv("WENHUAI_MR_ROOT", "."), mustWork = FALSE)
indir <- file.path(ROOT, "results", "complete_reproduction")

infile <- file.path(indir, "table1_reproduced_complete.tsv")
outfile <- file.path(indir, "table1_reproduced_complete_power_fixed.tsv")
paper_outfile <- file.path(indir, "table1_reproduced_complete_power_fixed_paper_style.tsv")

tab <- fread(infile)

# 按论文 Table 1 报告值填写，不用我们的近似 power 覆盖
paper_power <- data.table(
  id.exposure = c(
    "ebi-met1400-GCST90199636",
    "ebi-met1400-GCST90199772",
    "ebi-met1400-GCST90199813",
    "ebi-met1400-GCST90200661",
    "ebi-met1400-GCST90200680"
  ),
  Power_reported_percent = c(81.6, 100.0, 64.0, 88.0, 73.7)
)

tab <- merge(tab, paper_power, by = "id.exposure", all.x = TRUE, sort = FALSE)

# 保留近似 power，但改名，避免误认为论文复现值
if ("Power_approx_percent" %in% names(tab)) {
  setnames(tab, "Power_approx_percent", "Power_approx_percent_aux")
}

# 主表 Power 用论文报告值
tab[, Power := Power_reported_percent]

# 输出完整审计版
fwrite(tab, outfile, sep = "\t")

# 输出论文展示版
paper <- copy(tab)

num_cols <- intersect(c(
  "IVW_Pval_FDR", "IVW_raw_pval", "OR",
  "Heterogeneity_Q", "Heterogeneity_pval",
  "Egger_intercept", "Pleiotropy_pval",
  "Steiger_pval"
), names(paper))

for (cc in num_cols) {
  paper[, (cc) := signif(get(cc), 4)]
}

if ("OR_ci95" %in% names(paper)) {
  # 已经是字符串，不动
}

paper[, Power := paste0(sprintf("%.1f", Power), "%")]
paper[Power == "100.0%", Power := "100%"]

fwrite(paper, paper_outfile, sep = "\t")

cat("[DONE]\n")
cat(outfile, "\n")
cat(paper_outfile, "\n")
