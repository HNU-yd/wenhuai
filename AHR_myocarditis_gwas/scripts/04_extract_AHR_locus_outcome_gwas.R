
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

chr_target <- get_arg("--chr", "7")
region_start <- as.integer(get_arg("--start", "16200000"))
region_end <- as.integer(get_arg("--end", "18500000"))

data.table::setDTthreads(threads)

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

formatted_dir <- file.path(base, "data/formatted/outcome")
locus_dir <- file.path(base, "data/locus/outcome")
qc_dir <- file.path(base, "results/qc")
locus_result_dir <- file.path(base, "results/locus")
log_file <- file.path(base, "logs/04_extract_AHR_locus_outcome_gwas.log")

dir.create(locus_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(locus_result_dir, recursive = TRUE, showWarnings = FALSE)
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

extract_locus_one <- function(dataset_key, infile) {
  log_msg("============================================================")
  log_msg("Dataset: ", dataset_key)
  log_msg("Input: ", infile)

  if (!file.exists(infile)) {
    stop("Missing input file: ", infile)
  }

  region_tag <- paste0("chr", chr_target, "_", region_start, "_", region_end)

  raw_out <- file.path(
    locus_dir,
    paste0(dataset_key, ".AHR_locus.", region_tag, ".raw.tsv.gz")
  )

  dedup_out <- file.path(
    locus_dir,
    paste0(dataset_key, ".AHR_locus.", region_tag, ".dedup_by_variant.tsv.gz")
  )

  top_out <- file.path(
    locus_result_dir,
    paste0(dataset_key, ".AHR_locus.", region_tag, ".top50.tsv")
  )

  awk_script <- sprintf(
    'NR==1 || ($8=="%s" && $9 >= %d && $9 <= %d)',
    chr_target,
    region_start,
    region_end
  )

  cmd <- paste(decompress_cmd(infile), "| awk", shQuote(paste0('BEGIN{FS=OFS="\\t"} ', awk_script)))

  log_msg("Extract command prepared.")
  log_msg("Reading locus rows with threads: ", threads)

  dt <- fread(
    cmd = cmd,
    sep = "\t",
    header = TRUE,
    nThread = threads,
    showProgress = TRUE
  )

  log_msg("Raw locus rows: ", nrow(dt))

  if (nrow(dt) == 0) {
    warning("No rows extracted for ", dataset_key)
    return(data.table(
      dataset_key = dataset_key,
      chr = chr_target,
      region_start = region_start,
      region_end = region_end,
      raw_locus_file = raw_out,
      dedup_locus_file = dedup_out,
      top50_file = top_out,
      n_raw = 0,
      n_dedup_by_variant = 0,
      n_duplicate_variant_id = 0,
      min_pval = NA_real_,
      top_variant_id = NA_character_,
      top_snp = NA_character_,
      top_pos = NA_integer_,
      top_effect_allele = NA_character_,
      top_other_allele = NA_character_,
      top_beta = NA_real_,
      top_se = NA_real_,
      top_pval = NA_real_
    ))
  }

  dt[, analysis_region := paste0("AHR_locus_", region_tag)]

  fwrite(
    dt,
    raw_out,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip",
    nThread = threads
  )

  setorder(dt, variant_id, pval)
  dt_dedup <- dt[!duplicated(variant_id)]

  fwrite(
    dt_dedup,
    dedup_out,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip",
    nThread = threads
  )

  setorder(dt_dedup, pval)

  top_n <- min(50, nrow(dt_dedup))
  top_dt <- dt_dedup[seq_len(top_n)]

  fwrite(
    top_dt,
    top_out,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    nThread = threads
  )

  top1 <- dt_dedup[1]

  qc <- data.table(
    dataset_key = dataset_key,
    chr = chr_target,
    region_start = region_start,
    region_end = region_end,
    raw_locus_file = raw_out,
    dedup_locus_file = dedup_out,
    top50_file = top_out,
    n_raw = nrow(dt),
    n_dedup_by_variant = nrow(dt_dedup),
    n_duplicate_variant_id = nrow(dt) - uniqueN(dt$variant_id),
    min_pval = suppressWarnings(min(dt_dedup$pval, na.rm = TRUE)),
    top_variant_id = top1$variant_id,
    top_snp = top1$SNP,
    top_pos = top1$pos,
    top_effect_allele = top1$effect_allele,
    top_other_allele = top1$other_allele,
    top_beta = top1$beta,
    top_se = top1$se,
    top_pval = top1$pval
  )

  log_msg("Raw locus written: ", raw_out)
  log_msg("Dedup locus written: ", dedup_out)
  log_msg("Top50 written: ", top_out)

  rm(dt, dt_dedup, top_dt)
  gc()

  qc
}

main <- function() {
  log_msg("Start extracting AHR locus from formatted outcome GWAS.")
  log_msg("Region: chr", chr_target, ":", region_start, "-", region_end)
  log_msg("data.table threads: ", data.table::getDTthreads())

  files <- list(
    Sakaue2021_BBJ_Myocarditis = file.path(formatted_dir, "Sakaue2021_BBJ_Myocarditis.formatted.tsv.gz"),
    Sakaue2021_EUR_Myocarditis = file.path(formatted_dir, "Sakaue2021_EUR_Myocarditis.formatted.tsv.gz"),
    FinnGen_R12_I9_MYOCARD = file.path(formatted_dir, "FinnGen_R12_I9_MYOCARD.formatted.tsv.gz")
  )

  qc_list <- list()

  for (nm in names(files)) {
    qc_list[[nm]] <- extract_locus_one(nm, files[[nm]])
  }

  qc <- rbindlist(qc_list, fill = TRUE)

  qc_out <- file.path(
    qc_dir,
    paste0("AHR_locus_outcome_gwas_qc.chr", chr_target, "_", region_start, "_", region_end, ".tsv")
  )

  fwrite(qc, qc_out, sep = "\t", quote = FALSE, na = "NA", nThread = threads)

  log_msg("QC written: ", qc_out)
  log_msg("Done.")
}

main()
