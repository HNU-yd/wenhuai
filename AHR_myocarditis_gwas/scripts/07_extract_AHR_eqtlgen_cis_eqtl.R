
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

gene_target <- get_arg("--gene", "AHR")

data.table::setDTthreads(threads)

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

raw_dir <- file.path(base, "data/raw/exposure/eqtlgen")
out_dir <- file.path(base, "data/formatted/exposure")
qc_dir <- file.path(base, "results/qc")
res_dir <- file.path(base, "results/exposure")
log_file <- file.path(base, "logs/07_extract_AHR_eqtlgen_cis_eqtl.log")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
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

read_header <- function(infile) {
  cmd <- paste(decompress_cmd(infile), "| head -n 1")
  h <- fread(cmd = cmd, sep = "\t", header = TRUE, nrows = 0)
  names(h)
}

detect_col <- function(cols, candidates, required = FALSE, label = "") {
  for (x in candidates) {
    hit <- which(tolower(cols) == tolower(x))
    if (length(hit) > 0) return(cols[hit[1]])
  }

  if (required) {
    stop(
      "Cannot detect required column: ", label,
      "\nAvailable columns:\n",
      paste(seq_along(cols), cols, sep = ": ", collapse = "\n")
    )
  }

  NA_character_
}

extract_one <- function(file_key, infile, output_prefix) {
  log_msg("============================================================")
  log_msg("File key: ", file_key)
  log_msg("Input: ", infile)

  if (!file.exists(infile)) {
    stop("Missing input file: ", infile)
  }

  cols <- read_header(infile)

  gene_col <- detect_col(
    cols,
    c("GeneSymbol", "GeneName", "HGNCName", "Symbol", "Gene", "gene", "gene_name", "hgnc_symbol"),
    required = TRUE,
    label = "gene symbol"
  )

  p_col <- detect_col(
    cols,
    c("Pvalue", "PValue", "pvalue", "p.value", "P", "pval"),
    required = FALSE,
    label = "P value"
  )

  snp_col <- detect_col(
    cols,
    c("SNPName", "SNP", "rsid", "rsids", "variant", "variant_id"),
    required = FALSE,
    label = "SNP"
  )

  z_col <- detect_col(
    cols,
    c("OverallZScore", "ZScore", "Zscore", "Z", "zscore"),
    required = FALSE,
    label = "Z score"
  )

  allele_col <- detect_col(
    cols,
    c("AlleleAssessed", "AssessedAllele", "effect_allele", "EA"),
    required = FALSE,
    label = "assessed allele"
  )

  fdr_col <- detect_col(
    cols,
    c("FDR", "fdr"),
    required = FALSE,
    label = "FDR"
  )

  bonf_col <- detect_col(
    cols,
    c("BonferroniP", "Bonferroni", "bonferroni_p"),
    required = FALSE,
    label = "Bonferroni"
  )

  gene_idx <- which(cols == gene_col)

  log_msg("Detected gene column: ", gene_col, " ; column index = ", gene_idx)
  log_msg("Detected p column: ", p_col)
  log_msg("Detected SNP column: ", snp_col)
  log_msg("Detected Z column: ", z_col)
  log_msg("Detected allele column: ", allele_col)

  awk_body <- sprintf(
    'BEGIN{FS=OFS="\\t"} NR==1{print; next} $%d=="%s"{print}',
    gene_idx,
    gene_target
  )

  cmd <- paste(decompress_cmd(infile), "| awk", shQuote(awk_body))

  log_msg("Extracting rows for gene: ", gene_target)

  dt <- fread(
    cmd = cmd,
    sep = "\t",
    header = TRUE,
    nThread = threads,
    showProgress = TRUE
  )

  log_msg("Extracted rows: ", nrow(dt))

  if (nrow(dt) > 0) {
    dt[, source_file_key := file_key]
    dt[, target_gene := gene_target]
    dt[, gene_column_used := gene_col]
  }

  out_file <- file.path(out_dir, paste0(output_prefix, ".extracted.tsv.gz"))
  top_file <- file.path(res_dir, paste0(output_prefix, ".top50.tsv"))

  fwrite(
    dt,
    out_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip",
    nThread = threads
  )

  if (!is.na(p_col) && p_col %in% names(dt) && nrow(dt) > 0) {
    dt[, p_for_sort := suppressWarnings(as.numeric(get(p_col)))]
    setorder(dt, p_for_sort)
    top_dt <- dt[seq_len(min(50, .N))]
    top_dt[, p_for_sort := NULL]
  } else {
    top_dt <- dt[seq_len(min(50, nrow(dt)))]
  }

  fwrite(
    top_dt,
    top_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    nThread = threads
  )

  get_min_numeric <- function(col) {
    if (is.na(col) || !(col %in% names(dt)) || nrow(dt) == 0) return(NA_real_)
    suppressWarnings(min(as.numeric(dt[[col]]), na.rm = TRUE))
  }

  count_p_lt <- function(col, cutoff) {
    if (is.na(col) || !(col %in% names(dt)) || nrow(dt) == 0) return(NA_integer_)
    sum(suppressWarnings(as.numeric(dt[[col]])) < cutoff, na.rm = TRUE)
  }

  n_unique_snp <- if (!is.na(snp_col) && snp_col %in% names(dt)) uniqueN(dt[[snp_col]]) else NA_integer_

  qc <- data.table(
    file_key = file_key,
    input_file = infile,
    output_file = out_file,
    top50_file = top_file,
    target_gene = gene_target,
    gene_col = gene_col,
    p_col = p_col,
    snp_col = snp_col,
    z_col = z_col,
    allele_col = allele_col,
    fdr_col = fdr_col,
    bonferroni_col = bonf_col,
    n_rows = nrow(dt),
    n_unique_snp = n_unique_snp,
    min_pval = get_min_numeric(p_col),
    n_p_lt_5e_8 = count_p_lt(p_col, 5e-8),
    n_p_lt_1e_6 = count_p_lt(p_col, 1e-6),
    n_p_lt_1e_5 = count_p_lt(p_col, 1e-5),
    n_p_lt_1e_4 = count_p_lt(p_col, 1e-4)
  )

  log_msg("Written extracted file: ", out_file)
  log_msg("Written top50 file: ", top_file)

  rm(dt, top_dt)
  gc()

  qc
}

main <- function() {
  log_msg("Start extracting eQTLGen cis-eQTL for gene: ", gene_target)
  log_msg("data.table threads: ", data.table::getDTthreads())

  files <- list(
    eqtlgen_significant_FDR0.05 = file.path(
      raw_dir,
      "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
    ),
    eqtlgen_full_cis_eQTL = file.path(
      raw_dir,
      "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
    )
  )

  qc_list <- list(
    extract_one(
      "eqtlgen_significant_FDR0.05",
      files$eqtlgen_significant_FDR0.05,
      paste0("eqtlgen_", gene_target, "_FDR0.05_cis_eqtl")
    ),
    extract_one(
      "eqtlgen_full_cis_eQTL",
      files$eqtlgen_full_cis_eQTL,
      paste0("eqtlgen_", gene_target, "_full_cis_eqtl")
    )
  )

  qc <- rbindlist(qc_list, fill = TRUE)

  qc_out <- file.path(qc_dir, paste0("eqtlgen_", gene_target, "_cis_eqtl_extraction_qc.tsv"))

  fwrite(
    qc,
    qc_out,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    nThread = threads
  )

  log_msg("QC written: ", qc_out)
  log_msg("Done.")
}

main()
