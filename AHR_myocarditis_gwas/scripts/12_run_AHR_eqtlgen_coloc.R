
suppressPackageStartupMessages({
  library(data.table)
})

if (!requireNamespace("coloc", quietly = TRUE)) {
  stop("R package 'coloc' is not installed. Run: install.packages('coloc')")
}

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[hit + 1]
}

threads <- as.integer(get_arg("--threads", parallel::detectCores() - 1))
if (is.na(threads) || threads < 1) threads <- 1

finngen_ncase <- suppressWarnings(as.numeric(get_arg("--finngen_ncase", NA)))
finngen_ncontrol <- suppressWarnings(as.numeric(get_arg("--finngen_ncontrol", NA)))

data.table::setDTthreads(threads)

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

eqtl_file <- file.path(base, "data/formatted/exposure/eqtlgen_AHR_full_cis_eqtl.standardized.tsv.gz")

outcome_locus_dir <- file.path(base, "data/locus/outcome")
coloc_data_dir <- file.path(base, "data/coloc")
coloc_result_dir <- file.path(base, "results/coloc")
qc_dir <- file.path(base, "results/qc")
log_file <- file.path(base, "logs/12_run_AHR_eqtlgen_coloc.log")

dir.create(coloc_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(coloc_result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

read_gz_dt <- function(path) {
  fread(
    cmd = paste("zcat", shQuote(path)),
    sep = "\t",
    header = TRUE,
    nThread = threads,
    showProgress = FALSE
  )
}

is_rsid <- function(x) {
  grepl("^rs[0-9]+$", as.character(x))
}

is_palindromic <- function(a1, a2) {
  a1 <- toupper(as.character(a1))
  a2 <- toupper(as.character(a2))
  pair <- paste0(pmin(a1, a2), pmax(a1, a2))
  pair %in% c("AT", "CG")
}

prepare_eqtl <- function() {
  x <- read_gz_dt(eqtl_file)

  required <- c(
    "SNP", "chr", "pos", "effect_allele", "other_allele",
    "eaf", "beta", "se", "pval", "samplesize"
  )

  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop("Missing eQTL columns: ", paste(missing, collapse = ", "))
  }

  x <- x[, .(
    exp_SNP = as.character(SNP),
    exp_is_rsid = is_rsid(SNP),
    exp_chr = as.character(chr),
    exp_pos = as.integer(pos),
    exp_effect_allele = toupper(effect_allele),
    exp_other_allele = toupper(other_allele),
    exp_eaf = as.numeric(eaf),
    beta_eqtl = as.numeric(beta),
    se_eqtl = as.numeric(se),
    varbeta_eqtl = as.numeric(se)^2,
    pval_eqtl = as.numeric(pval),
    n_eqtl = as.numeric(samplesize),
    eqtl_z = as.numeric(z)
  )]

  x <- x[
    !is.na(beta_eqtl) &
    !is.na(se_eqtl) &
    !is.na(exp_eaf) &
    exp_eaf > 0 &
    exp_eaf < 1
  ]

  unique(x, by = "exp_SNP")
}

prepare_outcome <- function(path, dataset_name) {
  x <- read_gz_dt(path)

  required <- c(
    "dataset_key", "source_id", "open_gwas_id", "trait", "population",
    "chr", "pos", "SNP", "variant_id",
    "effect_allele", "other_allele",
    "eaf", "beta", "se", "pval",
    "n", "ncase", "ncontrol"
  )

  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop("Missing outcome columns in ", dataset_name, ": ", paste(missing, collapse = ", "))
  }

  x <- x[, .(
    outcome_dataset = dataset_key,
    outcome_source_id = source_id,
    outcome_open_gwas_id = open_gwas_id,
    outcome_trait = trait,
    outcome_population = population,
    outcome_chr = as.character(chr),
    outcome_pos = as.integer(pos),
    outcome_SNP = as.character(SNP),
    outcome_is_rsid = is_rsid(SNP),
    outcome_variant_id = as.character(variant_id),
    outcome_effect_allele = toupper(effect_allele),
    outcome_other_allele = toupper(other_allele),
    outcome_eaf = as.numeric(eaf),
    beta_gwas = as.numeric(beta),
    se_gwas = as.numeric(se),
    varbeta_gwas = as.numeric(se)^2,
    pval_gwas = as.numeric(pval),
    n_gwas = as.numeric(n),
    ncase_gwas = as.numeric(ncase),
    ncontrol_gwas = as.numeric(ncontrol)
  )]

  if (dataset_name == "FinnGen_R12_I9_MYOCARD") {
    if (!is.na(finngen_ncase) && !is.na(finngen_ncontrol)) {
      x[, ncase_gwas := finngen_ncase]
      x[, ncontrol_gwas := finngen_ncontrol]
      x[, n_gwas := finngen_ncase + finngen_ncontrol]
    }
  }

  x
}

hybrid_overlap <- function(eqtl, outcome) {
  rs_eqtl <- eqtl[exp_is_rsid == TRUE]
  rs_out <- outcome[outcome_is_rsid == TRUE]

  rs_match <- merge(
    rs_eqtl,
    rs_out,
    by.x = "exp_SNP",
    by.y = "outcome_SNP",
    allow.cartesian = TRUE
  )

  if (nrow(rs_match) > 0) {
    rs_match[, match_method := "rsID"]
    rs_match[, match_priority := 1L]
  }

  rs_matched <- unique(rs_match$exp_SNP)

  eqtl_left <- eqtl[!(exp_SNP %in% rs_matched)]

  pos_match <- merge(
    eqtl_left,
    outcome,
    by.x = c("exp_chr", "exp_pos"),
    by.y = c("outcome_chr", "outcome_pos"),
    allow.cartesian = TRUE
  )

  if (nrow(pos_match) > 0) {
    pos_match[, match_method := "chr_pos"]
    pos_match[, match_priority := 2L]
  }

  h <- rbindlist(list(rs_match, pos_match), fill = TRUE)

  if (nrow(h) == 0) return(h)

  h[, allele_match := fifelse(
    exp_effect_allele == outcome_effect_allele & exp_other_allele == outcome_other_allele,
    "aligned",
    fifelse(
      exp_effect_allele == outcome_other_allele & exp_other_allele == outcome_effect_allele,
      "swapped",
      NA_character_
    )
  )]

  h <- h[!is.na(allele_match)]

  if (nrow(h) == 0) return(h)

  h[, palindromic := is_palindromic(exp_effect_allele, exp_other_allele)]

  h[, ambiguous_palindromic := palindromic &
    (
      (!is.na(exp_eaf) & exp_eaf > 0.42 & exp_eaf < 0.58) |
      (!is.na(outcome_eaf) & outcome_eaf > 0.42 & outcome_eaf < 0.58)
    )
  ]

  h <- h[ambiguous_palindromic != TRUE]

  if (nrow(h) == 0) return(h)

  h[, beta_gwas_aligned := fifelse(allele_match == "aligned", beta_gwas, -beta_gwas)]
  h[, eaf_gwas_aligned := fifelse(allele_match == "aligned", outcome_eaf, 1 - outcome_eaf)]

  setorder(h, exp_SNP, match_priority, pval_gwas)
  h <- h[!duplicated(exp_SNP)]

  h[, snp_for_coloc := exp_SNP]
  h[!is_rsid(snp_for_coloc), snp_for_coloc := paste0(exp_chr, ":", exp_pos, ":", exp_effect_allele, ":", exp_other_allele)]

  h[, maf_for_coloc := pmin(exp_eaf, 1 - exp_eaf)]

  h <- h[
    !is.na(beta_eqtl) &
    !is.na(varbeta_eqtl) &
    !is.na(beta_gwas_aligned) &
    !is.na(varbeta_gwas) &
    !is.na(maf_for_coloc) &
    varbeta_eqtl > 0 &
    varbeta_gwas > 0 &
    maf_for_coloc > 0 &
    maf_for_coloc < 0.5
  ]

  h
}

run_coloc_for_dataset <- function(dataset_name, outcome_path, eqtl) {
  log_msg("============================================================")
  log_msg("Dataset: ", dataset_name)

  outcome <- prepare_outcome(outcome_path, dataset_name)
  overlap <- hybrid_overlap(eqtl, outcome)

  overlap_file <- file.path(coloc_data_dir, paste0("AHR_eqtlgen_", dataset_name, "_coloc_overlap.tsv.gz"))

  fwrite(
    overlap,
    overlap_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip",
    nThread = threads
  )

  log_msg("Overlap rows: ", nrow(overlap))
  log_msg("Overlap written: ", overlap_file)

  if (nrow(overlap) == 0) {
    return(list(
      summary = data.table(
        outcome_dataset = dataset_name,
        status = "skipped_no_overlap",
        n_overlap = 0
      ),
      snp = data.table()
    ))
  }

  ncase <- unique(na.omit(overlap$ncase_gwas))[1]
  ncontrol <- unique(na.omit(overlap$ncontrol_gwas))[1]

  if (is.na(ncase) || is.na(ncontrol)) {
    return(list(
      summary = data.table(
        outcome_dataset = dataset_name,
        status = "skipped_missing_ncase_ncontrol",
        n_overlap = nrow(overlap),
        overlap_file = overlap_file,
        ncase = NA_real_,
        ncontrol = NA_real_,
        s = NA_real_,
        PP.H0.abf = NA_real_,
        PP.H1.abf = NA_real_,
        PP.H2.abf = NA_real_,
        PP.H3.abf = NA_real_,
        PP.H4.abf = NA_real_
      ),
      snp = data.table()
    ))
  }

  n_total <- ncase + ncontrol
  s_case <- ncase / n_total
  n_eqtl <- round(stats::median(overlap$n_eqtl, na.rm = TRUE))

  d1 <- list(
    beta = overlap$beta_eqtl,
    varbeta = overlap$varbeta_eqtl,
    snp = overlap$snp_for_coloc,
    position = overlap$exp_pos,
    MAF = overlap$maf_for_coloc,
    N = n_eqtl,
    type = "quant"
  )

  d2 <- list(
    beta = overlap$beta_gwas_aligned,
    varbeta = overlap$varbeta_gwas,
    snp = overlap$snp_for_coloc,
    position = overlap$exp_pos,
    MAF = overlap$maf_for_coloc,
    N = n_total,
    s = s_case,
    type = "cc"
  )

  coloc_res <- coloc::coloc.abf(dataset1 = d1, dataset2 = d2)

  summary_dt <- as.data.table(as.list(coloc_res$summary))
  summary_dt[, outcome_dataset := dataset_name]
  summary_dt[, status := "ok"]
  summary_dt[, n_overlap := nrow(overlap)]
  summary_dt[, overlap_file := overlap_file]
  summary_dt[, ncase := ncase]
  summary_dt[, ncontrol := ncontrol]
  summary_dt[, s := s_case]
  setcolorder(summary_dt, c("outcome_dataset", "status", "n_overlap", "ncase", "ncontrol", "s", setdiff(names(summary_dt), c("outcome_dataset", "status", "n_overlap", "ncase", "ncontrol", "s"))))

  snp_dt <- as.data.table(coloc_res$results)
  snp_dt[, outcome_dataset := dataset_name]

  snp_anno <- overlap[, .(
    snp = snp_for_coloc,
    exp_SNP,
    exp_chr,
    exp_pos,
    match_method,
    exp_effect_allele,
    exp_other_allele,
    beta_eqtl,
    se_eqtl,
    pval_eqtl,
    beta_gwas_aligned,
    se_gwas,
    pval_gwas,
    outcome_variant_id
  )]

  snp_dt <- merge(snp_dt, snp_anno, by = "snp", all.x = TRUE, sort = FALSE)

  snp_file <- file.path(coloc_result_dir, paste0("AHR_eqtlgen_", dataset_name, "_coloc_snp_posterior.tsv"))
  fwrite(snp_dt, snp_file, sep = "\t", quote = FALSE, na = "NA", nThread = threads)

  summary_dt[, snp_posterior_file := snp_file]

  log_msg("Coloc PP.H4: ", summary_dt$PP.H4.abf)
  log_msg("SNP posterior written: ", snp_file)

  list(summary = summary_dt, snp = snp_dt)
}

main <- function() {
  log_msg("Start AHR eQTLGen coloc analysis.")
  log_msg("eQTL file: ", eqtl_file)

  eqtl <- prepare_eqtl()

  log_msg("AHR eQTL rows after basic filtering: ", nrow(eqtl))

  outcome_files <- c(
    Sakaue2021_BBJ_Myocarditis = file.path(outcome_locus_dir, "Sakaue2021_BBJ_Myocarditis.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz"),
    Sakaue2021_EUR_Myocarditis = file.path(outcome_locus_dir, "Sakaue2021_EUR_Myocarditis.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz"),
    FinnGen_R12_I9_MYOCARD = file.path(outcome_locus_dir, "FinnGen_R12_I9_MYOCARD.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz")
  )

  res <- list()

  for (nm in names(outcome_files)) {
    res[[nm]] <- run_coloc_for_dataset(nm, outcome_files[[nm]], eqtl)
  }

  summary_all <- rbindlist(lapply(res, `[[`, "summary"), fill = TRUE)

  summary_file <- file.path(coloc_result_dir, "AHR_eqtlgen_to_myocarditis_coloc_summary.tsv")
  fwrite(summary_all, summary_file, sep = "\t", quote = FALSE, na = "NA", nThread = threads)

  qc_file <- file.path(qc_dir, "AHR_eqtlgen_to_myocarditis_coloc_qc.tsv")
  fwrite(summary_all, qc_file, sep = "\t", quote = FALSE, na = "NA", nThread = threads)

  interpretation_file <- file.path(coloc_result_dir, "AHR_eqtlgen_to_myocarditis_coloc_interpretation_cn.md")

  lines <- c(
    "# AHR eQTLGen 与心肌炎 GWAS 的 coloc 共定位结果解释",
    "",
    "## 判读标准",
    "",
    "- PP.H4.abf 高：支持 AHR eQTL 和心肌炎 GWAS 共享同一因果变异。",
    "- PP.H3.abf 高：两个性状在该区域都有信号，但更可能是不同因果变异。",
    "- PP.H0/H1/H2 高：说明共定位证据不足，可能只有一个性状有信号或两个性状都缺乏明显信号。",
    "",
    "常用经验阈值：PP.H4.abf > 0.8 可视为较强共定位证据；0.5–0.8 为中等证据；低于 0.5 通常不支持强共定位。",
    "",
    "## 本次输出文件",
    "",
    paste0("- summary: ", summary_file),
    paste0("- QC: ", qc_file),
    "",
    "## 注意",
    "",
    "FinnGen 若显示 skipped_missing_ncase_ncontrol，说明其 locus overlap 已生成，但还缺少 FinnGen I9_MYOCARD 的病例数和对照数。补充 `--finngen_ncase` 与 `--finngen_ncontrol` 后可直接重跑。"
  )

  writeLines(lines, interpretation_file, useBytes = TRUE)

  log_msg("Summary written: ", summary_file)
  log_msg("QC written: ", qc_file)
  log_msg("Interpretation written: ", interpretation_file)
  log_msg("Done.")
}

main()
