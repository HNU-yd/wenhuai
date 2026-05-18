
suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (length(hit) == 0) return(default)
  if (hit == length(args)) return(default)
  args[hit + 1]
}

threads <- as.integer(get_arg("--threads", parallel::detectCores() - 1))
if (is.na(threads) || threads < 1) threads <- 1

data.table::setDTthreads(threads)

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

instrument_dir <- file.path(base, "data/instruments/exposure")
outcome_locus_dir <- file.path(base, "data/locus/outcome")
mr_dir <- file.path(base, "results/mr")
qc_dir <- file.path(base, "results/qc")
log_file <- file.path(base, "logs/10_run_AHR_eqtlgen_MR_against_myocarditis.log")

dir.create(mr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

read_dt <- function(path) {
  fread(
    path,
    sep = "\t",
    header = TRUE,
    nThread = threads,
    showProgress = FALSE
  )
}

read_gz_dt <- function(path) {
  cmd <- paste("zcat", shQuote(path))
  fread(
    cmd = cmd,
    sep = "\t",
    header = TRUE,
    nThread = threads,
    showProgress = FALSE
  )
}

is_palindromic <- function(a1, a2) {
  a1 <- toupper(as.character(a1))
  a2 <- toupper(as.character(a2))
  pair <- paste0(pmin(a1, a2), pmax(a1, a2))
  pair %in% c("AT", "CG")
}

p_from_z <- function(z) {
  2 * pnorm(abs(z), lower.tail = FALSE)
}

safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  suppressWarnings(min(x, na.rm = TRUE))
}

safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  suppressWarnings(median(x, na.rm = TRUE))
}

prepare_exposure <- function(path) {
  dt <- read_dt(path)

  setting_id <- basename(path)
  setting_id <- sub("^eqtlgen_AHR_full_cis_eqtl\\.", "", setting_id)
  setting_id <- sub("\\.clumped\\.tsv$", "", setting_id)

  required <- c("SNP", "chr", "pos", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "samplesize", "f_stat")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing exposure columns in ", path, ": ", paste(missing, collapse = ", "))
  }

  out <- dt[, .(
    instrument_setting = setting_id,
    instrument_file = path,
    exp_SNP = SNP,
    chr = as.character(chr),
    pos = as.integer(pos),
    exp_effect_allele = toupper(effect_allele),
    exp_other_allele = toupper(other_allele),
    exp_eaf = as.numeric(eaf),
    beta_exposure = as.numeric(beta),
    se_exposure = as.numeric(se),
    pval_exposure = as.numeric(pval),
    samplesize_exposure = as.numeric(samplesize),
    f_stat = as.numeric(f_stat),
    exposure = exposure,
    exposure_source = exposure_source
  )]

  out <- unique(out, by = "exp_SNP")
  out
}

prepare_outcome <- function(path) {
  dt <- read_gz_dt(path)

  required <- c("dataset_key", "source_id", "open_gwas_id", "trait", "population", "chr", "pos", "SNP", "variant_id", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "n", "ncase", "ncontrol")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing outcome columns in ", path, ": ", paste(missing, collapse = ", "))
  }

  out <- dt[, .(
    outcome_dataset = dataset_key,
    outcome_source_id = source_id,
    outcome_open_gwas_id = open_gwas_id,
    outcome_trait = trait,
    outcome_population = population,
    chr = as.character(chr),
    pos = as.integer(pos),
    outcome_SNP = SNP,
    outcome_variant_id = variant_id,
    outcome_effect_allele = toupper(effect_allele),
    outcome_other_allele = toupper(other_allele),
    outcome_eaf = as.numeric(eaf),
    beta_outcome = as.numeric(beta),
    se_outcome = as.numeric(se),
    pval_outcome = as.numeric(pval),
    n_outcome = as.numeric(n),
    ncase_outcome = as.numeric(ncase),
    ncontrol_outcome = as.numeric(ncontrol),
    outcome_locus_file = path
  )]

  out
}

harmonise_one <- function(exp_dt, out_dt) {
  n_exp <- nrow(exp_dt)

  merged <- merge(
    exp_dt,
    out_dt,
    by = c("chr", "pos"),
    allow.cartesian = TRUE
  )

  n_pos_match <- nrow(merged)

  if (nrow(merged) == 0) {
    return(list(
      harmonised = data.table(),
      qc = data.table(
        instrument_setting = unique(exp_dt$instrument_setting),
        outcome_dataset = unique(out_dt$outcome_dataset),
        n_exposure_instruments = n_exp,
        n_position_matched_rows = 0,
        n_allele_compatible_rows = 0,
        n_dropped_ambiguous_palindromic = 0,
        n_harmonised = 0,
        missing_exposure_snps = n_exp
      )
    ))
  }

  merged[, allele_match := fifelse(
    exp_effect_allele == outcome_effect_allele & exp_other_allele == outcome_other_allele,
    "aligned",
    fifelse(
      exp_effect_allele == outcome_other_allele & exp_other_allele == outcome_effect_allele,
      "swapped",
      NA_character_
    )
  )]

  allele_ok <- merged[!is.na(allele_match)]

  n_allele_ok <- nrow(allele_ok)

  if (nrow(allele_ok) == 0) {
    return(list(
      harmonised = data.table(),
      qc = data.table(
        instrument_setting = unique(exp_dt$instrument_setting),
        outcome_dataset = unique(out_dt$outcome_dataset),
        n_exposure_instruments = n_exp,
        n_position_matched_rows = n_pos_match,
        n_allele_compatible_rows = 0,
        n_dropped_ambiguous_palindromic = 0,
        n_harmonised = 0,
        missing_exposure_snps = n_exp
      )
    ))
  }

  allele_ok[, palindromic := is_palindromic(exp_effect_allele, exp_other_allele)]
  allele_ok[, ambiguous_palindromic := palindromic &
    (
      (!is.na(exp_eaf) & exp_eaf > 0.42 & exp_eaf < 0.58) |
      (!is.na(outcome_eaf) & outcome_eaf > 0.42 & outcome_eaf < 0.58)
    )
  ]

  n_drop_pal <- sum(allele_ok$ambiguous_palindromic, na.rm = TRUE)

  h <- allele_ok[ambiguous_palindromic != TRUE]

  if (nrow(h) > 0) {
    h[, beta_outcome_aligned := fifelse(allele_match == "aligned", beta_outcome, -beta_outcome)]
    h[, eaf_outcome_aligned := fifelse(allele_match == "aligned", outcome_eaf, 1 - outcome_eaf)]
    h[, effect_allele := exp_effect_allele]
    h[, other_allele := exp_other_allele]

    setorder(h, exp_SNP, pval_outcome)
    h <- h[!duplicated(exp_SNP)]

    h[, mr_keep := !is.na(beta_exposure) &
      !is.na(se_exposure) &
      !is.na(beta_outcome_aligned) &
      !is.na(se_outcome) &
      se_exposure > 0 &
      se_outcome > 0
    ]

    h <- h[mr_keep == TRUE]
  }

  qc <- data.table(
    instrument_setting = unique(exp_dt$instrument_setting),
    outcome_dataset = unique(out_dt$outcome_dataset),
    n_exposure_instruments = n_exp,
    n_position_matched_rows = n_pos_match,
    n_allele_compatible_rows = n_allele_ok,
    n_dropped_ambiguous_palindromic = n_drop_pal,
    n_harmonised = nrow(h),
    missing_exposure_snps = n_exp - uniqueN(h$exp_SNP)
  )

  list(harmonised = h, qc = qc)
}

run_mr_one <- function(h) {
  if (nrow(h) == 0) {
    return(data.table())
  }

  h[, ratio := beta_outcome_aligned / beta_exposure]
  h[, ratio_se := se_outcome / abs(beta_exposure)]
  h[, ratio_z := ratio / ratio_se]
  h[, ratio_pval := p_from_z(ratio_z)]
  h[, ratio_or := exp(ratio)]
  h[, ratio_or_lci95 := exp(ratio - 1.96 * ratio_se)]
  h[, ratio_or_uci95 := exp(ratio + 1.96 * ratio_se)]

  base_cols <- h[1, .(
    exposure = exposure,
    outcome_dataset = outcome_dataset,
    outcome_population = outcome_population,
    instrument_setting = instrument_setting
  )]

  nsnp <- nrow(h)

  res <- list()

  if (nsnp == 1) {
    beta <- h$ratio[1]
    se <- h$ratio_se[1]
    z <- beta / se

    res[[length(res) + 1]] <- cbind(
      base_cols,
      data.table(
        method = "Wald ratio",
        nsnp = nsnp,
        beta = beta,
        se = se,
        z = z,
        pval = p_from_z(z),
        OR = exp(beta),
        OR_lci95 = exp(beta - 1.96 * se),
        OR_uci95 = exp(beta + 1.96 * se),
        Q = NA_real_,
        Q_df = NA_integer_,
        Q_pval = NA_real_,
        egger_intercept = NA_real_,
        egger_intercept_se = NA_real_,
        egger_intercept_pval = NA_real_
      )
    )
  }

  if (nsnp >= 2) {
    w <- 1 / (h$ratio_se^2)
    beta <- sum(w * h$ratio) / sum(w)
    se <- sqrt(1 / sum(w))
    z <- beta / se
    Q <- sum(w * (h$ratio - beta)^2)
    Q_df <- nsnp - 1
    Q_p <- pchisq(Q, df = Q_df, lower.tail = FALSE)

    res[[length(res) + 1]] <- cbind(
      base_cols,
      data.table(
        method = "Inverse variance weighted",
        nsnp = nsnp,
        beta = beta,
        se = se,
        z = z,
        pval = p_from_z(z),
        OR = exp(beta),
        OR_lci95 = exp(beta - 1.96 * se),
        OR_uci95 = exp(beta + 1.96 * se),
        Q = Q,
        Q_df = Q_df,
        Q_pval = Q_p,
        egger_intercept = NA_real_,
        egger_intercept_se = NA_real_,
        egger_intercept_pval = NA_real_
      )
    )
  }

  if (nsnp >= 3) {
    fit <- tryCatch(
      lm(beta_outcome_aligned ~ beta_exposure, weights = 1 / se_outcome^2, data = h),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      sm <- summary(fit)
      co <- coef(sm)

      if (nrow(co) >= 2) {
        slope <- co["beta_exposure", "Estimate"]
        slope_se <- co["beta_exposure", "Std. Error"]
        slope_z <- slope / slope_se

        intercept <- co["(Intercept)", "Estimate"]
        intercept_se <- co["(Intercept)", "Std. Error"]
        intercept_p <- co["(Intercept)", "Pr(>|t|)"]

        res[[length(res) + 1]] <- cbind(
          base_cols,
          data.table(
            method = "MR Egger",
            nsnp = nsnp,
            beta = slope,
            se = slope_se,
            z = slope_z,
            pval = p_from_z(slope_z),
            OR = exp(slope),
            OR_lci95 = exp(slope - 1.96 * slope_se),
            OR_uci95 = exp(slope + 1.96 * slope_se),
            Q = NA_real_,
            Q_df = NA_integer_,
            Q_pval = NA_real_,
            egger_intercept = intercept,
            egger_intercept_se = intercept_se,
            egger_intercept_pval = intercept_p
          )
        )
      }
    }
  }

  rbindlist(res, fill = TRUE)
}

main <- function() {
  log_msg("Start harmonise and MR: AHR eQTLGen instruments vs myocarditis GWAS.")
  log_msg("Threads: ", threads)

  instrument_files <- list.files(
    instrument_dir,
    pattern = "^eqtlgen_AHR_full_cis_eqtl\\..*\\.clumped\\.tsv$",
    full.names = TRUE
  )

  if (length(instrument_files) == 0) {
    stop("No clumped instrument files found in: ", instrument_dir)
  }

  outcome_files <- c(
    Sakaue2021_BBJ_Myocarditis = file.path(outcome_locus_dir, "Sakaue2021_BBJ_Myocarditis.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz"),
    Sakaue2021_EUR_Myocarditis = file.path(outcome_locus_dir, "Sakaue2021_EUR_Myocarditis.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz"),
    FinnGen_R12_I9_MYOCARD = file.path(outcome_locus_dir, "FinnGen_R12_I9_MYOCARD.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz")
  )

  for (f in outcome_files) {
    if (!file.exists(f)) stop("Missing outcome locus file: ", f)
  }

  outcomes <- lapply(outcome_files, prepare_outcome)

  harmonised_list <- list()
  qc_list <- list()
  mr_list <- list()
  single_list <- list()

  idx <- 1

  for (inst_file in instrument_files) {
    exp_dt <- prepare_exposure(inst_file)

    log_msg("============================================================")
    log_msg("Instrument file: ", inst_file)
    log_msg("Instrument setting: ", unique(exp_dt$instrument_setting))
    log_msg("Instrument count: ", nrow(exp_dt))

    for (out_name in names(outcomes)) {
      out_dt <- outcomes[[out_name]]
      log_msg("Outcome: ", out_name)

      hx <- harmonise_one(exp_dt, out_dt)

      h <- hx$harmonised
      qc <- hx$qc

      harmonised_list[[idx]] <- h
      qc_list[[idx]] <- qc

      if (nrow(h) > 0) {
        mr <- run_mr_one(h)
        mr_list[[idx]] <- mr

        single <- h[, .(
          exposure,
          outcome_dataset,
          outcome_population,
          instrument_setting,
          exp_SNP,
          chr,
          pos,
          effect_allele,
          other_allele,
          beta_exposure,
          se_exposure,
          pval_exposure,
          f_stat,
          outcome_SNP,
          outcome_variant_id,
          beta_outcome_aligned,
          se_outcome,
          pval_outcome,
          allele_match,
          ratio,
          ratio_se,
          ratio_pval,
          ratio_or,
          ratio_or_lci95,
          ratio_or_uci95
        )]

        single_list[[idx]] <- single
      }

      log_msg("Harmonised SNPs: ", qc$n_harmonised)
      idx <- idx + 1
    }
  }

  harmonised <- rbindlist(harmonised_list, fill = TRUE)
  qc_all <- rbindlist(qc_list, fill = TRUE)
  mr_all <- rbindlist(mr_list, fill = TRUE)
  single_all <- rbindlist(single_list, fill = TRUE)

  harmonised_out <- file.path(mr_dir, "AHR_eqtlgen_to_myocarditis_harmonised.tsv")
  qc_out <- file.path(qc_dir, "AHR_eqtlgen_to_myocarditis_harmonise_qc.tsv")
  mr_out <- file.path(mr_dir, "AHR_eqtlgen_to_myocarditis_MR_results.tsv")
  single_out <- file.path(mr_dir, "AHR_eqtlgen_to_myocarditis_singleSNP.tsv")

  fwrite(harmonised, harmonised_out, sep = "\t", quote = FALSE, na = "NA")
  fwrite(qc_all, qc_out, sep = "\t", quote = FALSE, na = "NA")
  fwrite(mr_all, mr_out, sep = "\t", quote = FALSE, na = "NA")
  fwrite(single_all, single_out, sep = "\t", quote = FALSE, na = "NA")

  log_msg("Harmonised written: ", harmonised_out)
  log_msg("QC written: ", qc_out)
  log_msg("MR results written: ", mr_out)
  log_msg("Single SNP results written: ", single_out)
  log_msg("Done.")
}

main()
