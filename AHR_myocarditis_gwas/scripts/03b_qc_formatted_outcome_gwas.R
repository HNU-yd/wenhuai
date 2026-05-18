
# cd /home/data1/wenhuai/AHR_myocarditis_gwas
# Rscript scripts/03b_qc_formatted_outcome_gwas.R --threads 24

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
out_dir <- file.path(base, "data/formatted/outcome")
qc_dir <- file.path(base, "results/qc")
log_file <- file.path(base, "logs/03b_qc_formatted_outcome_gwas.log")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

summarise_one <- function(dataset_key, infile) {
  log_msg("Reading formatted GWAS for QC: ", dataset_key)
  log_msg("Input: ", infile)

  if (!file.exists(infile)) {
    stop("Missing formatted file: ", infile)
  }

  needed_cols <- c(
    "chr",
    "SNP",
    "variant_id",
    "effect_allele",
    "other_allele",
    "eaf",
    "beta",
    "se",
    "pval"
  )

  cmd <- paste("zcat", shQuote(infile))

  dt <- fread(
    cmd = cmd,
    sep = "\t",
    header = TRUE,
    select = needed_cols,
    nThread = threads,
    showProgress = TRUE
  )

  log_msg("Rows loaded: ", nrow(dt))

  qc <- data.table(
    dataset_key = dataset_key,
    formatted_file = infile,
    n_rows = nrow(dt),
    n_chr = uniqueN(dt$chr),
    min_pval = suppressWarnings(min(dt$pval, na.rm = TRUE)),
    n_missing_snp = sum(is.na(dt$SNP) | dt$SNP == ""),
    n_missing_beta = sum(is.na(dt$beta)),
    n_missing_se = sum(is.na(dt$se)),
    n_missing_pval = sum(is.na(dt$pval)),
    n_missing_eaf = sum(is.na(dt$eaf)),
    n_duplicate_variant_id = nrow(dt) - uniqueN(dt$variant_id),
    n_duplicate_snp = nrow(dt) - uniqueN(dt$SNP),
    n_snv = sum(
      nchar(dt$effect_allele) == 1 &
      nchar(dt$other_allele) == 1,
      na.rm = TRUE
    ),
    n_indel = sum(
      nchar(dt$effect_allele) != 1 |
      nchar(dt$other_allele) != 1,
      na.rm = TRUE
    )
  )

  rm(dt)
  gc()

  qc
}

log_msg("Start QC for formatted outcome GWAS.")
log_msg("data.table threads: ", data.table::getDTthreads())

files <- list(
  Sakaue2021_BBJ_Myocarditis = file.path(out_dir, "Sakaue2021_BBJ_Myocarditis.formatted.tsv.gz"),
  Sakaue2021_EUR_Myocarditis = file.path(out_dir, "Sakaue2021_EUR_Myocarditis.formatted.tsv.gz"),
  FinnGen_R12_I9_MYOCARD = file.path(out_dir, "FinnGen_R12_I9_MYOCARD.formatted.tsv.gz")
)

qc_list <- list()

for (nm in names(files)) {
  qc_list[[nm]] <- summarise_one(nm, files[[nm]])
}

qc <- rbindlist(qc_list, fill = TRUE)

qc_out <- file.path(qc_dir, "formatted_outcome_gwas_qc.tsv")
fwrite(qc, qc_out, sep = "\t", quote = FALSE, na = "NA", nThread = threads)

log_msg("QC written: ", qc_out)
log_msg("Done.")
