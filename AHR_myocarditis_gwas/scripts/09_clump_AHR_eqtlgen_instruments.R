
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

threads <- as.integer(get_arg("--threads", "8"))
if (is.na(threads) || threads < 1) threads <- 1

plink_bin <- get_arg("--plink", "plink")
bfile <- get_arg("--bfile", "")

data.table::setDTthreads(threads)

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

exposure_file <- file.path(
  base,
  "data/formatted/exposure/eqtlgen_AHR_full_cis_eqtl.standardized.tsv.gz"
)

out_dir <- file.path(base, "data/instruments/exposure")
qc_dir <- file.path(base, "results/qc")
tmp_dir <- file.path(base, "results/instruments/tmp_clump")
log_file <- file.path(base, "logs/09_clump_AHR_eqtlgen_instruments.log")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

decompress_cmd <- function(infile) {
  pigz <- Sys.which("pigz")
  if (nzchar(pigz)) {
    paste(shQuote(pigz), "-dc -p", threads, shQuote(infile))
  } else {
    paste("zcat", shQuote(infile))
  }
}

read_exposure <- function(infile) {
  fread(
    cmd = decompress_cmd(infile),
    sep = "\t",
    header = TRUE,
    nThread = threads,
    showProgress = TRUE
  )
}

tag_num <- function(x) {
  y <- format(x, scientific = FALSE, trim = TRUE)
  y <- gsub("\\.", "p", y)
  y
}

run_one_clump <- function(dt, p_thresh, r2, kb) {
  tag <- paste0("p", tag_num(p_thresh), "_r2_", tag_num(r2), "_kb", kb)

  log_msg("============================================================")
  log_msg("Clumping setting: ", tag)

  x <- dt[
    !is.na(SNP) &
      SNP != "" &
      !is.na(pval) &
      pval < p_thresh &
      !is.na(beta) &
      !is.na(se) &
      !is.na(eaf)
  ]

  setorder(x, pval)
  x <- unique(x, by = "SNP")

  log_msg("Candidate SNPs before PLINK clumping: ", nrow(x))

  clump_input <- file.path(tmp_dir, paste0("AHR_eqtlgen_", tag, ".clump_input.tsv"))
  clump_prefix <- file.path(tmp_dir, paste0("AHR_eqtlgen_", tag))
  clump_output <- paste0(clump_prefix, ".clumped")

  out_clumped <- file.path(out_dir, paste0("eqtlgen_AHR_full_cis_eqtl.", tag, ".clumped.tsv"))
  out_removed <- file.path(out_dir, paste0("eqtlgen_AHR_full_cis_eqtl.", tag, ".removed_or_not_in_ref.tsv"))

  fwrite(
    x[, .(SNP = SNP, P = pval)],
    clump_input,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  cmd_args <- c(
    "--bfile", bfile,
    "--clump", clump_input,
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-p1", as.character(p_thresh),
    "--clump-p2", "1",
    "--clump-r2", as.character(r2),
    "--clump-kb", as.character(kb),
    "--threads", as.character(threads),
    "--out", clump_prefix
  )

  log_msg("Running PLINK command:")
  log_msg(plink_bin, " ", paste(cmd_args, collapse = " "))

  plink_out <- system2(
    command = plink_bin,
    args = cmd_args,
    stdout = TRUE,
    stderr = TRUE
  )

  cat(paste(plink_out, collapse = "\n"), "\n", file = log_file, append = TRUE)

  if (!file.exists(clump_output)) {
    log_msg("[warning] PLINK .clumped output not found: ", clump_output)
    clumped <- x[0]
  } else {
    clumped_raw <- fread(clump_output, fill = TRUE)

    if (!"SNP" %in% names(clumped_raw) || nrow(clumped_raw) == 0) {
      clumped <- x[0]
    } else {
      keep <- unique(clumped_raw$SNP)
      clumped <- x[SNP %in% keep]
      clumped[, clump_rank := match(SNP, keep)]
      setorder(clumped, clump_rank)
    }
  }

  removed <- x[!(SNP %in% clumped$SNP)]

  fwrite(clumped, out_clumped, sep = "\t", quote = FALSE, na = "NA")
  fwrite(removed, out_removed, sep = "\t", quote = FALSE, na = "NA")

  log_msg("Clumped SNPs: ", nrow(clumped))
  log_msg("Clumped output: ", out_clumped)

  data.table(
    exposure = "AHR_expression_eQTLGen_whole_blood",
    bfile = bfile,
    plink_bin = plink_bin,
    p_thresh = p_thresh,
    clump_r2 = r2,
    clump_kb = kb,
    n_candidate = nrow(x),
    n_clumped = nrow(clumped),
    n_removed_or_not_in_ref = nrow(removed),
    clumped_file = out_clumped,
    removed_file = out_removed,
    min_p_candidate = suppressWarnings(min(x$pval, na.rm = TRUE)),
    min_p_clumped = suppressWarnings(min(clumped$pval, na.rm = TRUE)),
    median_f_stat_clumped = suppressWarnings(median(clumped$f_stat, na.rm = TRUE)),
    min_f_stat_clumped = suppressWarnings(min(clumped$f_stat, na.rm = TRUE)),
    max_f_stat_clumped = suppressWarnings(max(clumped$f_stat, na.rm = TRUE)),
    top_snp = ifelse(nrow(clumped) > 0, clumped[1]$SNP, NA_character_),
    top_beta = ifelse(nrow(clumped) > 0, clumped[1]$beta, NA_real_),
    top_se = ifelse(nrow(clumped) > 0, clumped[1]$se, NA_real_),
    top_pval = ifelse(nrow(clumped) > 0, clumped[1]$pval, NA_real_)
  )
}

main <- function() {
  log_msg("Start AHR eQTLGen LD clumping.")
  log_msg("Exposure file: ", exposure_file)
  log_msg("PLINK: ", plink_bin)
  log_msg("BFILE: ", bfile)
  log_msg("Threads: ", threads)

  if (bfile == "") {
    stop("Missing --bfile")
  }

  if (!file.exists(paste0(bfile, ".bed")) ||
      !file.exists(paste0(bfile, ".bim")) ||
      !file.exists(paste0(bfile, ".fam"))) {
    stop("Missing PLINK reference files for bfile prefix: ", bfile)
  }

  if (!file.exists(exposure_file)) {
    stop("Missing exposure file: ", exposure_file)
  }

  dt <- read_exposure(exposure_file)

  log_msg("Loaded exposure rows: ", nrow(dt))
  log_msg("p < 5e-8: ", sum(dt$pval < 5e-8, na.rm = TRUE))
  log_msg("p < 1e-6: ", sum(dt$pval < 1e-6, na.rm = TRUE))
  log_msg("p < 1e-5: ", sum(dt$pval < 1e-5, na.rm = TRUE))

  settings <- rbindlist(list(
    data.table(p_thresh = 5e-8, r2 = 0.001, kb = 10000),
    data.table(p_thresh = 1e-6, r2 = 0.001, kb = 10000),
    data.table(p_thresh = 1e-5, r2 = 0.001, kb = 10000),
    data.table(p_thresh = 5e-8, r2 = 0.01, kb = 10000),
    data.table(p_thresh = 1e-6, r2 = 0.01, kb = 10000),
    data.table(p_thresh = 1e-5, r2 = 0.01, kb = 10000)
  ))

  qc_list <- vector("list", nrow(settings))

  for (i in seq_len(nrow(settings))) {
    qc_list[[i]] <- run_one_clump(
      dt = dt,
      p_thresh = settings$p_thresh[i],
      r2 = settings$r2[i],
      kb = settings$kb[i]
    )
  }

  qc <- rbindlist(qc_list, fill = TRUE)

  qc_out <- file.path(qc_dir, "eqtlgen_AHR_ld_clump_qc.tsv")
  fwrite(qc, qc_out, sep = "\t", quote = FALSE, na = "NA")

  log_msg("QC written: ", qc_out)
  log_msg("Done.")
}

main()
