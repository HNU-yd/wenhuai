
suppressPackageStartupMessages({
  library(data.table)
})

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

mr_file <- file.path(base, "results/mr/AHR_eqtlgen_to_myocarditis_hybrid_MR_results.tsv")
qc_file <- file.path(base, "results/qc/AHR_eqtlgen_to_myocarditis_hybrid_harmonise_qc.tsv")

out_dir <- file.path(base, "results/mr")
fig_dir <- file.path(base, "results/figures")
log_file <- file.path(base, "logs/11_summarize_AHR_eqtlgen_MR_results.log")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

setting_label <- function(x) {
  x <- gsub("p0p00000005", "p<5e-8", x)
  x <- gsub("p0p000001", "p<1e-6", x)
  x <- gsub("p0p00001", "p<1e-5", x)
  x <- gsub("_r2_0p001", ", r2<0.001", x)
  x <- gsub("_r2_0p01", ", r2<0.01", x)
  x <- gsub("_kb10000", ", 10Mb", x)
  x
}

interpret_direction <- function(beta, pval) {
  ifelse(
    is.na(beta) | is.na(pval),
    "NA",
    ifelse(
      pval < 0.05 & beta > 0,
      "nominal_positive",
      ifelse(
        pval < 0.05 & beta < 0,
        "nominal_negative",
        ifelse(beta > 0, "positive_not_significant", "negative_not_significant")
      )
    )
  )
}

format_ci <- function(or, lci, uci) {
  sprintf("%.3f (%.3f–%.3f)", or, lci, uci)
}

main <- function() {
  log_msg("Start summarizing AHR eQTLGen MR results.")

  if (!file.exists(mr_file)) stop("Missing MR file: ", mr_file)
  if (!file.exists(qc_file)) stop("Missing harmonise QC file: ", qc_file)

  mr <- fread(mr_file, sep = "\t", header = TRUE)
  qc <- fread(qc_file, sep = "\t", header = TRUE)

  log_msg("MR rows: ", nrow(mr))
  log_msg("QC rows: ", nrow(qc))

  primary_methods <- c("Wald ratio", "Inverse variance weighted")

  primary <- mr[method %in% primary_methods]

  primary[, setting_label := setting_label(instrument_setting)]
  primary[, OR_CI95 := format_ci(OR, OR_lci95, OR_uci95)]
  primary[, direction_summary := interpret_direction(beta, pval)]

  primary <- primary[, .(
    exposure,
    outcome_dataset,
    outcome_population,
    instrument_setting,
    setting_label,
    method,
    nsnp,
    beta,
    se,
    z,
    pval,
    OR,
    OR_lci95,
    OR_uci95,
    OR_CI95,
    Q,
    Q_df,
    Q_pval,
    direction_summary
  )]

  setorder(primary, outcome_dataset, instrument_setting, method)

  all_primary_out <- file.path(out_dir, "AHR_eqtlgen_MR_primary_methods_summary.tsv")
  fwrite(primary, all_primary_out, sep = "\t", quote = FALSE, na = "NA")
  log_msg("Written: ", all_primary_out)

  strict_setting <- "p0p00000005_r2_0p001_kb10000"

  strict <- primary[instrument_setting == strict_setting]

  strict_out <- file.path(out_dir, "AHR_eqtlgen_MR_strict_main_result.tsv")
  fwrite(strict, strict_out, sep = "\t", quote = FALSE, na = "NA")
  log_msg("Written: ", strict_out)

  sensitivity <- primary[, .(
    n_result_rows = .N,
    n_min = min(nsnp, na.rm = TRUE),
    n_max = max(nsnp, na.rm = TRUE),
    beta_min = min(beta, na.rm = TRUE),
    beta_max = max(beta, na.rm = TRUE),
    OR_min = min(OR, na.rm = TRUE),
    OR_max = max(OR, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE),
    n_p_lt_0p05 = sum(pval < 0.05, na.rm = TRUE),
    n_positive = sum(beta > 0, na.rm = TRUE),
    n_negative = sum(beta < 0, na.rm = TRUE)
  ), by = .(outcome_dataset, outcome_population)]

  sensitivity[, conclusion := fifelse(
    n_p_lt_0p05 == 0,
    "No nominally significant MR result across tested instrument settings",
    "At least one nominally significant MR result; inspect sensitivity and pleiotropy"
  )]

  sens_out <- file.path(out_dir, "AHR_eqtlgen_MR_sensitivity_by_outcome.tsv")
  fwrite(sensitivity, sens_out, sep = "\t", quote = FALSE, na = "NA")
  log_msg("Written: ", sens_out)

  qc_sum <- qc[, .(
    n_settings = .N,
    min_n_exposure_instruments = min(n_exposure_instruments, na.rm = TRUE),
    max_n_exposure_instruments = max(n_exposure_instruments, na.rm = TRUE),
    min_n_harmonised = min(n_harmonised, na.rm = TRUE),
    max_n_harmonised = max(n_harmonised, na.rm = TRUE),
    total_dropped_ambiguous_palindromic = sum(n_dropped_ambiguous_palindromic, na.rm = TRUE)
  ), by = .(outcome_dataset)]

  qc_sum_out <- file.path(out_dir, "AHR_eqtlgen_MR_harmonise_summary_by_outcome.tsv")
  fwrite(qc_sum, qc_sum_out, sep = "\t", quote = FALSE, na = "NA")
  log_msg("Written: ", qc_sum_out)

  # Forest plot: strict main setting
  strict_plot <- strict[order(outcome_dataset)]
  strict_fig <- file.path(fig_dir, "AHR_eqtlgen_MR_strict_main_forest.png")

  if (nrow(strict_plot) > 0) {
    png(strict_fig, width = 1600, height = 900, res = 160)

    par(mar = c(5, 12, 4, 2))

    y <- seq_len(nrow(strict_plot))
    x_min <- min(log(strict_plot$OR_lci95), na.rm = TRUE)
    x_max <- max(log(strict_plot$OR_uci95), na.rm = TRUE)
    pad <- 0.15 * (x_max - x_min)

    plot(
      NA,
      xlim = c(x_min - pad, x_max + pad),
      ylim = c(0.5, length(y) + 0.5),
      yaxt = "n",
      xlab = "MR estimate, log(OR)",
      ylab = "",
      main = "AHR expression eQTLGen whole blood → Myocarditis GWAS"
    )

    abline(v = 0, lty = 2)

    segments(
      x0 = log(strict_plot$OR_lci95),
      y0 = y,
      x1 = log(strict_plot$OR_uci95),
      y1 = y
    )

    points(log(strict_plot$OR), y, pch = 19)

    axis(
      side = 2,
      at = y,
      labels = paste0(
        strict_plot$outcome_dataset,
        " | ",
        strict_plot$method,
        " | nsnp=",
        strict_plot$nsnp
      ),
      las = 2
    )

    text(
      x = x_max + 0.05 * (x_max - x_min),
      y = y,
      labels = paste0("OR=", sprintf("%.3f", strict_plot$OR), ", P=", signif(strict_plot$pval, 3)),
      pos = 4,
      cex = 0.8
    )

    dev.off()
    log_msg("Written: ", strict_fig)
  }

  # Forest plot: all primary settings
  all_plot <- copy(primary)
  all_plot[, label := paste0(outcome_dataset, " | ", setting_label, " | ", method, " | nsnp=", nsnp)]
  all_plot <- all_plot[order(outcome_dataset, instrument_setting)]

  all_fig <- file.path(fig_dir, "AHR_eqtlgen_MR_all_primary_forest.png")

  if (nrow(all_plot) > 0) {
    png(all_fig, width = 2200, height = max(1200, 55 * nrow(all_plot)), res = 160)

    par(mar = c(5, 18, 4, 2))

    y <- seq_len(nrow(all_plot))
    x_min <- min(log(all_plot$OR_lci95), na.rm = TRUE)
    x_max <- max(log(all_plot$OR_uci95), na.rm = TRUE)
    pad <- 0.15 * (x_max - x_min)

    plot(
      NA,
      xlim = c(x_min - pad, x_max + pad),
      ylim = c(0.5, length(y) + 0.5),
      yaxt = "n",
      xlab = "MR estimate, log(OR)",
      ylab = "",
      main = "AHR eQTLGen MR sensitivity across instrument settings"
    )

    abline(v = 0, lty = 2)

    segments(
      x0 = log(all_plot$OR_lci95),
      y0 = y,
      x1 = log(all_plot$OR_uci95),
      y1 = y
    )

    points(log(all_plot$OR), y, pch = 19)

    axis(side = 2, at = y, labels = all_plot$label, las = 2, cex.axis = 0.65)

    dev.off()
    log_msg("Written: ", all_fig)
  }

  interpretation_file <- file.path(out_dir, "AHR_eqtlgen_MR_interpretation_cn.md")

  interpretation <- c(
    "# AHR eQTLGen whole-blood cis-eQTL 与心肌炎 GWAS 的 MR 结果解释",
    "",
    "## 核心结论",
    "",
    "在三套并列主分析心肌炎 GWAS 中，AHR whole-blood cis-eQTL 工具变量未显示稳定且显著的遗传因果证据。",
    "",
    "严格工具变量集合定义为：p < 5e-8，LD clumping r2 < 0.001，窗口 10000 kb。",
    "",
    "在该严格集合下：",
    "",
    "- BBJ 数据集仅匹配 1 个 SNP，只能进行 Wald ratio，结果不显著。",
    "- Sakaue European 数据集可匹配 3 个 SNP，IVW OR 约为 1.10，结果不显著。",
    "- FinnGen R12 数据集通过 rsID 匹配恢复 3 个 SNP，IVW OR 约为 0.92，结果不显著。",
    "",
    "## 方向一致性",
    "",
    "EUR 数据集多数设置下效应方向偏正，FinnGen 数据集多数设置下效应方向偏负，BBJ 数据集由于可用工具变量较少且置信区间较宽，估计不稳定。",
    "",
    "因此，当前 MR 结果不能支持“遗传预测的 AHR 表达升高会显著增加或降低心肌炎风险”的结论。",
    "",
    "## 后续分析重点",
    "",
    "MR 结果阴性或方向不一致并不排除 AHR locus 在局部遗传调控中的作用。下一步应使用 full AHR cis-eQTL 和三套心肌炎 GWAS 的 AHR locus summary statistics 进行 coloc 共定位分析，判断 AHR eQTL 信号和心肌炎 GWAS 局部信号是否共享同一因果变异。"
  )

  writeLines(interpretation, interpretation_file, useBytes = TRUE)
  log_msg("Written: ", interpretation_file)

  log_msg("Done.")
}

main()
