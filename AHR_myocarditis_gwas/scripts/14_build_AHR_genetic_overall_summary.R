
suppressPackageStartupMessages({
  library(data.table)
})

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

mr_strict_file <- file.path(base, "results/mr/AHR_eqtlgen_MR_strict_main_result.tsv")
mr_sens_file <- file.path(base, "results/mr/AHR_eqtlgen_MR_sensitivity_by_outcome.tsv")
coloc_file <- file.path(base, "results/coloc/AHR_eqtlgen_coloc_main_result.tsv")

out_dir <- file.path(base, "results/overall")
log_file <- file.path(base, "logs/14_build_AHR_genetic_overall_summary.log")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

classify_mr <- function(pval, beta) {
  if (is.na(pval) || is.na(beta)) return("not_available")
  if (pval < 0.05 && beta > 0) return("nominal_positive")
  if (pval < 0.05 && beta < 0) return("nominal_negative")
  if (beta > 0) return("positive_not_significant")
  if (beta < 0) return("negative_not_significant")
  "null_or_unclear"
}

classify_coloc <- function(pph4) {
  if (is.na(pph4)) return("not_available")
  if (pph4 >= 0.8) return("strong_colocalization")
  if (pph4 >= 0.5) return("moderate_colocalization")
  "no_strong_colocalization"
}

main <- function() {
  log_msg("Start building AHR genetic overall summary.")

  if (!file.exists(mr_strict_file)) stop("Missing file: ", mr_strict_file)
  if (!file.exists(mr_sens_file)) stop("Missing file: ", mr_sens_file)
  if (!file.exists(coloc_file)) stop("Missing file: ", coloc_file)

  mr <- fread(mr_strict_file, sep = "\t")
  sens <- fread(mr_sens_file, sep = "\t")
  coloc <- fread(coloc_file, sep = "\t")

  mr_main <- mr[, .(
    outcome_dataset,
    outcome_population,
    mr_method = method,
    mr_nsnp = nsnp,
    mr_beta = beta,
    mr_se = se,
    mr_pval = pval,
    mr_OR = OR,
    mr_OR_lci95 = OR_lci95,
    mr_OR_uci95 = OR_uci95
  )]

  mr_main[, mr_interpretation := mapply(classify_mr, mr_pval, mr_beta)]

  coloc_main <- coloc[, .(
    outcome_dataset,
    coloc_status = status,
    coloc_n_overlap = n_overlap,
    coloc_ncase = ncase,
    coloc_ncontrol = ncontrol,
    coloc_PP_H1 = PP.H1.abf,
    coloc_PP_H3 = PP.H3.abf,
    coloc_PP_H4 = PP.H4.abf
  )]

  coloc_main[, coloc_interpretation := sapply(coloc_PP_H4, classify_coloc)]

  merged <- merge(
    mr_main,
    coloc_main,
    by = "outcome_dataset",
    all = TRUE
  )

  merged[, genetic_layer_conclusion := fifelse(
    mr_pval < 0.05 & coloc_PP_H4 >= 0.8,
    "MR_and_coloc_support_direct_AHR_expression_effect",
    fifelse(
      mr_pval >= 0.05 & coloc_PP_H4 < 0.5,
      "no_direct_genetic_support_for_AHR_expression_effect",
      "mixed_or_weak_genetic_evidence"
    )
  )]

  setorder(merged, outcome_dataset)

  evidence_out <- file.path(out_dir, "AHR_genetic_evidence_matrix.tsv")
  fwrite(merged, evidence_out, sep = "\t", quote = FALSE, na = "NA")
  log_msg("Written: ", evidence_out)

  summary_lines <- c(
    "# AHR 遗传层整体分析总结",
    "",
    "## 1. 分析对象",
    "",
    "本部分整合 AHR whole-blood cis-eQTL 与三套并列主分析心肌炎 GWAS 的结果，包括 MR 和 coloc 两类分析。",
    "",
    "三套 outcome 数据为：",
    "",
    "- Sakaue2021_BBJ_Myocarditis",
    "- Sakaue2021_EUR_Myocarditis",
    "- FinnGen_R12_I9_MYOCARD",
    "",
    "## 2. MR 总体结论",
    "",
    "严格工具变量集合下，三套数据均未显示稳定且显著的 AHR 表达遗传因果效应。",
    "",
    "Sakaue European 数据中效应方向偏正，但不显著；FinnGen 数据中效应方向偏负，也不显著；BBJ 数据可用工具变量较少，估计不稳定。",
    "",
    "## 3. Coloc 总体结论",
    "",
    "三套数据的 PP.H4 均低于 0.1，未达到通常认为支持共定位的阈值。",
    "",
    "这说明 AHR eQTL 信号虽然很强，但心肌炎 GWAS 在 AHR locus 的局部信号并未与 AHR eQTL 共享同一因果变异。",
    "",
    "## 4. 遗传层最终定性",
    "",
    "当前遗传证据不支持“外周血 AHR 表达遗传调控本身是心肌炎风险的直接驱动因子”。",
    "",
    "因此，AHR 不宜作为遗传因果层面的核心阳性结论。更合适的定位是：AHR 是犬尿氨酸代谢相关免疫通路的关键受体，可能在疾病状态下作为免疫代谢响应节点参与心肌炎过程。",
    "",
    "## 5. 对整体课题的影响",
    "",
    "后续整体分析应将重点从“AHR 表达遗传因果效应”转向“犬尿氨酸代谢-AHR 通路活动”。也就是说，遗传层重点看犬尿氨酸/色氨酸代谢物 MR；转录组层重点看 KYN-AHR score、AHR response score、髓系炎症、趋化和抗原呈递信号。",
    "",
    "## 6. 对论文写法的建议",
    "",
    "不建议写：",
    "",
    "> AHR 表达升高导致心肌炎风险增加。",
    "",
    "建议写：",
    "",
    "> 遗传分析未支持 AHR 表达 cis-eQTL 作为心肌炎风险的直接因果因素；然而，结合犬尿氨酸代谢和转录组证据，AHR 通路仍可能作为疾病状态下的免疫代谢响应轴参与心肌炎炎症过程。"
  )

  summary_out <- file.path(out_dir, "AHR_genetic_overall_summary_cn.md")
  writeLines(summary_lines, summary_out, useBytes = TRUE)
  log_msg("Written: ", summary_out)

  log_msg("Done.")
}

main()
