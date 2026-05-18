
suppressPackageStartupMessages({
  library(data.table)
})

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

summary_file <- file.path(base, "results/coloc/AHR_eqtlgen_to_myocarditis_coloc_summary.tsv")
coloc_dir <- file.path(base, "results/coloc")
fig_dir <- file.path(base, "results/figures")
log_file <- file.path(base, "logs/13_summarize_AHR_eqtlgen_coloc_results.log")

dir.create(coloc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

classify_coloc <- function(pph4, pph3) {
  if (is.na(pph4)) return("not_available")
  if (pph4 >= 0.8) return("strong_colocalization")
  if (pph4 >= 0.5) return("moderate_colocalization")
  if (!is.na(pph3) && pph3 > pph4 && pph3 >= 0.5) return("distinct_signals_likely")
  return("no_strong_colocalization")
}

main <- function() {
  log_msg("Start summarizing AHR coloc results.")

  if (!file.exists(summary_file)) {
    stop("Missing coloc summary file: ", summary_file)
  }

  s <- fread(summary_file, sep = "\t", header = TRUE)

  s[, PP.H0.abf := as.numeric(PP.H0.abf)]
  s[, PP.H1.abf := as.numeric(PP.H1.abf)]
  s[, PP.H2.abf := as.numeric(PP.H2.abf)]
  s[, PP.H3.abf := as.numeric(PP.H3.abf)]
  s[, PP.H4.abf := as.numeric(PP.H4.abf)]

  s[, coloc_interpretation := mapply(classify_coloc, PP.H4.abf, PP.H3.abf)]

  main_table <- s[, .(
    outcome_dataset,
    status,
    n_overlap,
    ncase,
    ncontrol,
    case_fraction = s,
    nsnps,
    PP.H0.abf,
    PP.H1.abf,
    PP.H2.abf,
    PP.H3.abf,
    PP.H4.abf,
    coloc_interpretation,
    overlap_file,
    snp_posterior_file
  )]

  main_out <- file.path(coloc_dir, "AHR_eqtlgen_coloc_main_result.tsv")
  fwrite(main_table, main_out, sep = "\t", quote = FALSE, na = "NA")
  log_msg("Written: ", main_out)

  top_list <- list()

  for (i in seq_len(nrow(main_table))) {
    dataset <- main_table$outcome_dataset[i]
    f <- main_table$snp_posterior_file[i]

    if (is.na(f) || !file.exists(f)) next

    x <- fread(f, sep = "\t", header = TRUE)

    pp_col <- grep("SNP.*PP.*H4|SNP.PP.H4", names(x), value = TRUE)

    if (length(pp_col) == 0) {
      log_msg("No SNP.PP.H4-like column found in: ", f)
      next
    }

    pp_col <- pp_col[1]

    x[, snp_pp_h4_for_sort := as.numeric(get(pp_col))]
    setorder(x, -snp_pp_h4_for_sort)

    keep_cols <- intersect(
      c(
        "outcome_dataset",
        "snp",
        "exp_SNP",
        "exp_chr",
        "exp_pos",
        "match_method",
        "exp_effect_allele",
        "exp_other_allele",
        "beta_eqtl",
        "se_eqtl",
        "pval_eqtl",
        "beta_gwas_aligned",
        "se_gwas",
        "pval_gwas",
        "outcome_variant_id",
        pp_col
      ),
      names(x)
    )

    top <- x[seq_len(min(20, .N)), ..keep_cols]
    top[, posterior_column := pp_col]
    top_list[[dataset]] <- top
  }

  if (length(top_list) > 0) {
    top_all <- rbindlist(top_list, fill = TRUE)
    top_out <- file.path(coloc_dir, "AHR_eqtlgen_coloc_top_snp_posterior.tsv")
    fwrite(top_all, top_out, sep = "\t", quote = FALSE, na = "NA")
    log_msg("Written: ", top_out)
  }

  fig_file <- file.path(fig_dir, "AHR_eqtlgen_coloc_posterior_barplot.png")

  plot_dt <- main_table[status == "ok"]

  if (nrow(plot_dt) > 0) {
    png(fig_file, width = 1600, height = 900, res = 160)

    mat <- t(as.matrix(plot_dt[, .(
      PP.H0.abf,
      PP.H1.abf,
      PP.H2.abf,
      PP.H3.abf,
      PP.H4.abf
    )]))

    colnames(mat) <- plot_dt$outcome_dataset

    barplot(
      mat,
      beside = TRUE,
      ylim = c(0, 1),
      las = 2,
      ylab = "Posterior probability",
      main = "AHR eQTLGen whole-blood cis-eQTL coloc with myocarditis GWAS"
    )

    legend(
      "topright",
      legend = rownames(mat),
      bty = "n",
      cex = 0.8
    )

    dev.off()
    log_msg("Written: ", fig_file)
  }

  interpretation_file <- file.path(coloc_dir, "AHR_eqtlgen_coloc_interpretation_cn.final.md")

  lines <- c(
    "# AHR eQTLGen 与心肌炎 GWAS coloc 共定位结果总结",
    "",
    "## 1. 核心结论",
    "",
    "三套并列主分析心肌炎 GWAS 中，均未观察到 AHR whole-blood cis-eQTL 与心肌炎 GWAS 在 AHR locus 共享同一因果变异的强证据。",
    "",
    "当前经验判读标准为：",
    "",
    "- PP.H4.abf > 0.8：强共定位证据；",
    "- PP.H4.abf = 0.5–0.8：中等共定位证据；",
    "- PP.H4.abf < 0.5：通常不支持强共定位。",
    "",
    "## 2. 三套数据库结果",
    ""
  )

  for (i in seq_len(nrow(main_table))) {
    row <- main_table[i]
    lines <- c(
      lines,
      paste0("### ", row$outcome_dataset),
      "",
      paste0("- status: ", row$status),
      paste0("- n_overlap: ", row$n_overlap),
      paste0("- PP.H1.abf: ", signif(row$PP.H1.abf, 4)),
      paste0("- PP.H3.abf: ", signif(row$PP.H3.abf, 4)),
      paste0("- PP.H4.abf: ", signif(row$PP.H4.abf, 4)),
      paste0("- interpretation: ", row$coloc_interpretation),
      ""
    )
  }

  lines <- c(
    lines,
    "## 3. 生物学解释",
    "",
    "AHR 在 eQTLGen whole blood 中具有非常强的 cis-eQTL 信号，但这些调控 AHR 表达的遗传变异并未与心肌炎 GWAS 的 AHR 区域局部信号形成强共定位。",
    "",
    "因此，当前遗传证据不支持“外周血 AHR 表达遗传调控本身是心肌炎风险的直接驱动因子”。",
    "",
    "这并不否定 AHR 通路在心肌炎中的作用。更合理的解释是：AHR 可能更多体现为疾病状态下的免疫代谢响应或下游通路活动，而不是由 AHR 基因表达 cis-eQTL 单独驱动的遗传易感机制。",
    "",
    "## 4. 与 MR 结果的关系",
    "",
    "MR 分析未观察到稳定显著的 AHR 表达遗传因果效应；coloc 分析进一步显示 AHR eQTL 与心肌炎 GWAS 局部信号缺乏强共定位。两者方向一致，均提示 AHR 表达遗传调控证据较弱。",
    "",
    "后续文章中应避免写成“AHR 表达升高导致心肌炎风险升高”。更合适的表述是：遗传层面未支持 AHR 表达本身作为心肌炎风险的直接因果因子，但 AHR 通路仍可结合犬尿氨酸代谢 MR、单细胞和空间转录组结果作为疾病相关免疫代谢通路进行讨论。"
  )

  writeLines(lines, interpretation_file, useBytes = TRUE)
  log_msg("Written: ", interpretation_file)

  log_msg("Done.")
}

main()
