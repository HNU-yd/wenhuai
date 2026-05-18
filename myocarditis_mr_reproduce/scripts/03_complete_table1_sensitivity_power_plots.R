suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ggplot2)
})

ROOT <- normalizePath(Sys.getenv("WENHUAI_MR_ROOT", "."), mustWork = FALSE)

indir <- file.path(ROOT, "results", "clean_outcome_mr")
outdir <- file.path(ROOT, "results", "complete_reproduction")
plotdir <- file.path(ROOT, "plots", "complete_reproduction")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

harm_file <- file.path(indir, "07_harmonised_mr_keep.tsv")
mr_file <- file.path(indir, "10_final5_mr_results.tsv")
ivw_file <- file.path(indir, "11_final5_ivw_only.tsv")
het_file <- file.path(indir, "12_heterogeneity.tsv")
pleio_file <- file.path(indir, "13_pleiotropy_egger_intercept.tsv")
single_file <- file.path(indir, "14_single_snp.tsv")
loo_file <- file.path(indir, "15_leave_one_out.tsv")

full1400_ivw_file <- file.path(ROOT, "results", "full1400", "mr", "10_full1400_ivw_only.tsv")

targets <- data.table(
  id.exposure = c(
    "ebi-met1400-GCST90199636",
    "ebi-met1400-GCST90199772",
    "ebi-met1400-GCST90199813",
    "ebi-met1400-GCST90200661",
    "ebi-met1400-GCST90200680"
  ),
  expected_nsnp_paper = c(5, 7, 5, 6, 4),
  paper_exposure_order = 1:5
)

stopifnot(file.exists(harm_file))
stopifnot(file.exists(mr_file))
stopifnot(file.exists(ivw_file))

dat <- fread(harm_file)
mr_res <- fread(mr_file)
ivw <- fread(ivw_file)

dat <- dat[id.exposure %in% targets$id.exposure]
mr_res <- mr_res[id.exposure %in% targets$id.exposure]
ivw <- ivw[id.exposure %in% targets$id.exposure]

# ------------------------------------------------------------
# 1. FDR: 优先使用全量 1400 的 FDR
# ------------------------------------------------------------
if (file.exists(full1400_ivw_file)) {
  full_ivw <- fread(full1400_ivw_file)
  full_ivw <- full_ivw[id.exposure %in% targets$id.exposure]

  fdr_map <- full_ivw[, .(
    id.exposure,
    IVW_Pval_FDR_full1400 = if ("p_fdr_ivw_full1400" %in% names(full_ivw)) p_fdr_ivw_full1400 else p.adjust(pval, "fdr")
  )]

  ivw <- merge(ivw, fdr_map, by = "id.exposure", all.x = TRUE)

  ivw[, IVW_Pval_FDR := IVW_Pval_FDR_full1400]
  ivw[, FDR_source := "full1400"]
} else {
  if ("p_fdr_ivw_final5" %in% names(ivw)) {
    ivw[, IVW_Pval_FDR := p_fdr_ivw_final5]
  } else {
    ivw[, IVW_Pval_FDR := p.adjust(pval, "fdr")]
  }
  ivw[, FDR_source := "final5_only"]
}

# ------------------------------------------------------------
# 2. Heterogeneity
# ------------------------------------------------------------
if (file.exists(het_file)) {
  het <- fread(het_file)
  het <- het[id.exposure %in% targets$id.exposure]
  het_ivw <- het[method == "Inverse variance weighted"]

  het_tab <- het_ivw[, .(
    id.exposure,
    Heterogeneity_Q = Q,
    Heterogeneity_pval = Q_pval
  )]
} else {
  het_tab <- data.table(
    id.exposure = targets$id.exposure,
    Heterogeneity_Q = NA_real_,
    Heterogeneity_pval = NA_real_
  )
}

# ------------------------------------------------------------
# 3. Egger pleiotropy
# ------------------------------------------------------------
if (file.exists(pleio_file)) {
  pleio <- fread(pleio_file)
  pleio <- pleio[id.exposure %in% targets$id.exposure]

  pleio_tab <- pleio[, .(
    id.exposure,
    Egger_intercept = egger_intercept,
    Pleiotropy_pval = pval
  )]
} else {
  pleio_tab <- data.table(
    id.exposure = targets$id.exposure,
    Egger_intercept = NA_real_,
    Pleiotropy_pval = NA_real_
  )
}

# ------------------------------------------------------------
# 4. Steiger test
# ------------------------------------------------------------
dat_steiger <- copy(dat)

# FinnGen R10 I9_MYOCARD: paper reports 1654 cases, 210652 controls, total 212306
dat_steiger[, samplesize.exposure := 8299]
dat_steiger[, samplesize.outcome := 212306]
dat_steiger[, ncase.outcome := 1654]
dat_steiger[, ncontrol.outcome := 210652]

steiger_res <- tryCatch({
  x <- directionality_test(as.data.frame(dat_steiger))
  as.data.table(x)
}, error = function(e) {
  message("[Steiger ERROR] ", conditionMessage(e))
  data.table()
})

if (nrow(steiger_res) > 0) {
  fwrite(steiger_res, file.path(outdir, "steiger_test_raw.tsv"), sep = "\t")

  # 常见列名：correct_causal_direction, steiger_pval
  steiger_cols <- names(steiger_res)

  dir_col <- if ("correct_causal_direction" %in% steiger_cols) {
    "correct_causal_direction"
  } else if ("steiger_dir" %in% steiger_cols) {
    "steiger_dir"
  } else {
    NA_character_
  }

  p_col <- if ("steiger_pval" %in% steiger_cols) {
    "steiger_pval"
  } else if ("pval" %in% steiger_cols) {
    "pval"
  } else {
    NA_character_
  }

  steiger_tab <- unique(steiger_res[, .(
    id.exposure,
    Steiger_direction = if (!is.na(dir_col)) as.character(get(dir_col)) else NA_character_,
    Steiger_pval = if (!is.na(p_col)) as.numeric(get(p_col)) else NA_real_
  )])
} else {
  steiger_tab <- data.table(
    id.exposure = targets$id.exposure,
    Steiger_direction = NA_character_,
    Steiger_pval = NA_real_
  )
}

# ------------------------------------------------------------
# 5. Approximate MR power
# ------------------------------------------------------------
# 说明：
# 论文使用在线 power calculator；这里用 case-control MR 的近似 NCP 公式重算。
# R2 使用 2*EAF*(1-EAF)*beta_exposure^2 近似。
approx_power_binary <- function(beta_log_or, r2, ncase, ncontrol, alpha = 0.05) {
  if (is.na(beta_log_or) || is.na(r2) || r2 <= 0) return(NA_real_)
  neff <- 4 / (1 / ncase + 1 / ncontrol)
  ncp <- neff * r2 * beta_log_or^2 / max(1e-12, 1 - r2)
  crit <- qchisq(1 - alpha, df = 1)
  power <- 1 - pchisq(crit, df = 1, ncp = ncp)
  return(power)
}

dat_power <- copy(dat)
dat_power[, eaf_for_r2 := eaf.exposure]
dat_power[is.na(eaf_for_r2) & !is.na(eaf.outcome), eaf_for_r2 := eaf.outcome]

dat_power[, snp_r2_exposure_approx := 2 * eaf_for_r2 * (1 - eaf_for_r2) * beta.exposure^2]
dat_power[is.na(snp_r2_exposure_approx) | snp_r2_exposure_approx < 0, snp_r2_exposure_approx := 0]

r2_tab <- dat_power[, .(
  exposure_r2_approx = sum(snp_r2_exposure_approx, na.rm = TRUE)
), by = .(id.exposure)]

ivw_beta <- ivw[, .(
  id.exposure,
  ivw_beta_log_or = b
)]

power_tab <- merge(r2_tab, ivw_beta, by = "id.exposure", all.x = TRUE)
power_tab[, Power_approx := mapply(
  approx_power_binary,
  beta_log_or = ivw_beta_log_or,
  r2 = exposure_r2_approx,
  MoreArgs = list(ncase = 1654, ncontrol = 210652, alpha = 0.05)
)]
power_tab[, Power_approx_percent := Power_approx * 100]

fwrite(power_tab, file.path(outdir, "power_approximation.tsv"), sep = "\t")

# ------------------------------------------------------------
# 6. MR-PRESSO
# ------------------------------------------------------------
mrpresso_tab <- data.table(
  id.exposure = targets$id.exposure,
  MR_PRESSO_global_p = NA_character_,
  MR_PRESSO_outlier_snps = NA_character_,
  MR_PRESSO_status = "not_run"
)

if (requireNamespace("MRPRESSO", quietly = TRUE)) {
  presso_list <- list()

  for (id in targets$id.exposure) {
    sub <- dat[id.exposure == id]
    sub <- as.data.frame(sub)

    if (nrow(sub) < 4) {
      mrpresso_tab[id.exposure == id, `:=`(
        MR_PRESSO_status = "skipped_nsnp_lt_4",
        MR_PRESSO_global_p = NA_character_,
        MR_PRESSO_outlier_snps = NA_character_
      )]
      next
    }

    message("[MR-PRESSO] ", id, " nsnp=", nrow(sub))

    res <- tryCatch({
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = sub,
        NbDistribution = 10000,
        SignifThreshold = 0.05
      )
    }, error = function(e) {
      message("[MR-PRESSO ERROR] ", id, " | ", conditionMessage(e))
      NULL
    })

    if (is.null(res)) {
      mrpresso_tab[id.exposure == id, MR_PRESSO_status := "error"]
      next
    }

    saveRDS(res, file.path(outdir, paste0("mrpresso_", gsub("[^A-Za-z0-9]+", "_", id), ".rds")))

    global_p <- NA_character_
    outliers <- NA_character_

    if (!is.null(res[["MR-PRESSO results"]][["Global Test"]])) {
      gt <- res[["MR-PRESSO results"]][["Global Test"]]
      if ("Pvalue" %in% names(gt)) {
        global_p <- as.character(gt$Pvalue)
      }
    }

    if (!is.null(res[["MR-PRESSO results"]][["Outlier Test"]])) {
      ot <- as.data.table(res[["MR-PRESSO results"]][["Outlier Test"]], keep.rownames = "rowid")
      fwrite(ot, file.path(outdir, paste0("mrpresso_outlier_test_", gsub("[^A-Za-z0-9]+", "_", id), ".tsv")), sep = "\t")

      # 尽量识别 outlier 行
      pcols <- intersect(c("Pvalue", "P-value", "pvalue", "pval"), names(ot))
      if (length(pcols) > 0) {
        pc <- pcols[1]
        suppressWarnings(ot[, outlier_p := as.numeric(get(pc))])
        outlier_rows <- ot[!is.na(outlier_p) & outlier_p < 0.05]
        if (nrow(outlier_rows) > 0) {
          idx <- suppressWarnings(as.integer(outlier_rows$rowid))
          idx <- idx[!is.na(idx) & idx >= 1 & idx <= nrow(sub)]
          if (length(idx) > 0) outliers <- paste(sub$SNP[idx], collapse = ",")
        } else {
          outliers <- ""
        }
      }
    }

    mrpresso_tab[id.exposure == id, `:=`(
      MR_PRESSO_status = "ok",
      MR_PRESSO_global_p = global_p,
      MR_PRESSO_outlier_snps = outliers
    )]
  }
} else {
  message("[MR-PRESSO] package not installed; skip.")
  mrpresso_tab[, MR_PRESSO_status := "MRPRESSO_not_installed"]
}

fwrite(mrpresso_tab, file.path(outdir, "mrpresso_summary.tsv"), sep = "\t")

# ------------------------------------------------------------
# 7. Final Table 1
# ------------------------------------------------------------
table1 <- ivw[, .(
  id.exposure,
  Outcome = outcome,
  Exposure = exposure,
  nSNP = nsnp,
  IVW_Pval_FDR,
  IVW_raw_pval = pval,
  OR = or,
  OR_lci95 = or_lci95,
  OR_uci95 = or_uci95,
  OR_ci95 = paste0(sprintf("%.3f", or_lci95), "-", sprintf("%.3f", or_uci95)),
  FDR_source
)]

table1 <- merge(table1, het_tab, by = "id.exposure", all.x = TRUE)
table1 <- merge(table1, pleio_tab, by = "id.exposure", all.x = TRUE)
table1 <- merge(table1, steiger_tab, by = "id.exposure", all.x = TRUE)
table1 <- merge(table1, power_tab[, .(id.exposure, exposure_r2_approx, Power_approx_percent)], by = "id.exposure", all.x = TRUE)
table1 <- merge(table1, mrpresso_tab, by = "id.exposure", all.x = TRUE)
table1 <- merge(table1, targets, by = "id.exposure", all.x = TRUE)

setorder(table1, paper_exposure_order)

final <- table1[, .(
  Outcome,
  Exposure,
  nSNP,
  expected_nsnp_paper,
  IVW_Pval_FDR,
  IVW_raw_pval,
  OR,
  OR_ci95,
  Heterogeneity_Q,
  Heterogeneity_pval,
  Egger_intercept,
  Pleiotropy_pval,
  Steiger_direction,
  Steiger_pval,
  Power_approx_percent,
  exposure_r2_approx,
  MR_PRESSO_global_p,
  MR_PRESSO_outlier_snps,
  MR_PRESSO_status,
  FDR_source,
  id.exposure
)]

fwrite(final, file.path(outdir, "table1_reproduced_complete.tsv"), sep = "\t")

# 也输出一个更接近论文展示格式的版本
paper_style <- copy(final)
paper_style[, OR := sprintf("%.3f", OR)]
paper_style[, IVW_Pval_FDR := signif(IVW_Pval_FDR, 4)]
paper_style[, Heterogeneity_Q := signif(Heterogeneity_Q, 4)]
paper_style[, Heterogeneity_pval := signif(Heterogeneity_pval, 4)]
paper_style[, Pleiotropy_pval := signif(Pleiotropy_pval, 4)]
paper_style[, Steiger_pval := signif(Steiger_pval, 4)]
paper_style[, Power_approx_percent := paste0(sprintf("%.1f", Power_approx_percent), "%")]

fwrite(paper_style, file.path(outdir, "table1_reproduced_complete_paper_style.tsv"), sep = "\t")

# ------------------------------------------------------------
# 8. Plots
# ------------------------------------------------------------
dat_df <- as.data.frame(dat)
mr_res_df <- as.data.frame(mr_res)

plot_pdf <- function(plot_list, file, width = 8, height = 6) {
  pdf(file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  if (inherits(plot_list, "ggplot")) {
    print(plot_list)
  } else if (is.list(plot_list)) {
    for (p in plot_list) {
      if (!is.null(p)) print(p)
    }
  } else {
    print(plot_list)
  }
}

tryCatch({
  p <- mr_scatter_plot(mr_res_df, dat_df)
  plot_pdf(p, file.path(plotdir, "final5_scatter.pdf"))
}, error = function(e) message("[plot scatter ERROR] ", conditionMessage(e)))

single <- tryCatch({
  mr_singlesnp(dat_df)
}, error = function(e) {
  message("[single SNP ERROR] ", conditionMessage(e))
  NULL
})

if (!is.null(single)) {
  fwrite(as.data.table(single), file.path(outdir, "single_snp_for_plots.tsv"), sep = "\t")

  tryCatch({
    p <- mr_funnel_plot(single)
    plot_pdf(p, file.path(plotdir, "final5_funnel.pdf"))
  }, error = function(e) message("[plot funnel ERROR] ", conditionMessage(e)))

  tryCatch({
    p <- mr_forest_plot(single)
    plot_pdf(p, file.path(plotdir, "final5_forest.pdf"))
  }, error = function(e) message("[plot forest ERROR] ", conditionMessage(e)))
}

loo <- tryCatch({
  mr_leaveoneout(dat_df)
}, error = function(e) {
  message("[leave-one-out ERROR] ", conditionMessage(e))
  NULL
})

if (!is.null(loo)) {
  fwrite(as.data.table(loo), file.path(outdir, "leave_one_out_for_plots.tsv"), sep = "\t")

  tryCatch({
    p <- mr_leaveoneout_plot(loo)
    plot_pdf(p, file.path(plotdir, "final5_leave_one_out.pdf"))
  }, error = function(e) message("[plot leave-one-out ERROR] ", conditionMessage(e)))
}

cat("\n[DONE]\n")
cat("[TABLE1]\n")
cat(file.path(outdir, "table1_reproduced_complete.tsv"), "\n")
cat(file.path(outdir, "table1_reproduced_complete_paper_style.tsv"), "\n")
cat("[PLOTS]\n")
cat(plotdir, "\n")
