
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

in_dir <- file.path(base, "data/formatted/exposure")
raw_eqtlgen_dir <- file.path(base, "data/raw/exposure/eqtlgen")
out_dir <- file.path(base, "data/formatted/exposure")
qc_dir <- file.path(base, "results/qc")
res_dir <- file.path(base, "results/exposure")
log_file <- file.path(base, "logs/08_standardize_AHR_eqtlgen_exposure.log")

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

read_gz_dt <- function(infile, select = NULL) {
  cmd <- decompress_cmd(infile)
  fread(
    cmd = cmd,
    sep = "\t",
    header = TRUE,
    select = select,
    nThread = threads,
    showProgress = TRUE
  )
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

read_header_cols <- function(infile) {
  h <- fread(
    cmd = paste(decompress_cmd(infile), "| head -n 1"),
    sep = "\t",
    header = TRUE,
    nrows = 0
  )
  names(h)
}

standardize_one <- function(file_key, eqtl_file, af_file, output_prefix) {
  log_msg("============================================================")
  log_msg("Standardizing: ", file_key)
  log_msg("eQTL file: ", eqtl_file)
  log_msg("AF file: ", af_file)

  eqtl <- read_gz_dt(eqtl_file)

  required_eqtl <- c(
    "Pvalue",
    "SNP",
    "SNPChr",
    "SNPPos",
    "AssessedAllele",
    "OtherAllele",
    "Zscore",
    "Gene",
    "GeneSymbol",
    "GeneChr",
    "GenePos",
    "NrSamples"
  )

  missing_eqtl <- setdiff(required_eqtl, names(eqtl))
  if (length(missing_eqtl) > 0) {
    stop("Missing eQTL columns: ", paste(missing_eqtl, collapse = ", "))
  }

  af_cols <- read_header_cols(af_file)

  snp_col <- detect_col(
    af_cols,
    c("SNP", "SNPName", "rsid", "rsID", "rsids", "MarkerName"),
    required = TRUE,
    label = "SNP in allele frequency file"
  )

  allele_a_col <- detect_col(
    af_cols,
    c("AlleleA", "Allele_A", "A1", "Allele1", "OtherAllele"),
    required = FALSE,
    label = "AlleleA"
  )

  allele_b_col <- detect_col(
    af_cols,
    c("AlleleB", "Allele_B", "A2", "Allele2", "AssessedAllele"),
    required = FALSE,
    label = "AlleleB"
  )

  af_b_col <- detect_col(
    af_cols,
    c(
      "AlleleB_all",
      "AlleleB_All",
      "AF_AlleleB",
      "AlleleBFrequency",
      "AlleleB_frequency",
      "Freq_AlleleB",
      "FreqAlleleB",
      "AF",
      "MAF"
    ),
    required = TRUE,
    label = "AlleleB frequency"
  )

  log_msg("AF SNP column: ", snp_col)
  log_msg("AF AlleleA column: ", allele_a_col)
  log_msg("AF AlleleB column: ", allele_b_col)
  log_msg("AF AlleleB frequency column: ", af_b_col)

  af_select <- unique(na.omit(c(snp_col, allele_a_col, allele_b_col, af_b_col)))

  af <- read_gz_dt(af_file, select = af_select)

  setnames(af, snp_col, "SNP")
  if (!is.na(allele_a_col) && allele_a_col %in% names(af)) {
    setnames(af, allele_a_col, "AF_AlleleA")
  } else {
    af[, AF_AlleleA := NA_character_]
  }

  if (!is.na(allele_b_col) && allele_b_col %in% names(af)) {
    setnames(af, allele_b_col, "AF_AlleleB")
  } else {
    af[, AF_AlleleB := NA_character_]
  }

  setnames(af, af_b_col, "AF_for_AlleleB")
  af[, AF_for_AlleleB := suppressWarnings(as.numeric(AF_for_AlleleB))]

  af <- unique(af, by = "SNP")

  log_msg("eQTL rows: ", nrow(eqtl))
  log_msg("AF rows: ", nrow(af))

  dt <- merge(eqtl, af, by = "SNP", all.x = TRUE, sort = FALSE)

  dt[, Pvalue := suppressWarnings(as.numeric(Pvalue))]
  dt[, Zscore := suppressWarnings(as.numeric(Zscore))]
  dt[, NrSamples := suppressWarnings(as.numeric(NrSamples))]
  dt[, SNPChr := as.character(SNPChr)]
  dt[, SNPPos := as.integer(SNPPos)]

  dt[, eaf := NA_real_]

  has_allele_cols <- all(c("AF_AlleleA", "AF_AlleleB") %in% names(dt))

  if (has_allele_cols) {
    dt[AssessedAllele == AF_AlleleB, eaf := AF_for_AlleleB]
    dt[AssessedAllele == AF_AlleleA, eaf := 1 - AF_for_AlleleB]
  }

  dt[is.na(eaf), eaf := AF_for_AlleleB]

  dt[eaf < 0 | eaf > 1, eaf := NA_real_]

  dt[, beta := NA_real_]
  dt[, se := NA_real_]

  dt[, valid_beta_se_calc := !is.na(Zscore) &
    !is.na(NrSamples) &
    !is.na(eaf) &
    eaf > 0 &
    eaf < 1]

  dt[valid_beta_se_calc == TRUE, se := 1 / sqrt(2 * eaf * (1 - eaf) * (NrSamples + Zscore^2))]
  dt[valid_beta_se_calc == TRUE, beta := Zscore * se]
  dt[, valid_beta_se_calc := NULL]

  dt[, exposure := paste0(gene_target, "_expression_eQTLGen_whole_blood")]
  dt[, exposure_source := "eQTLGen whole blood cis-eQTL"]
  dt[, effect_allele := AssessedAllele]
  dt[, other_allele := OtherAllele]
  dt[, pval := Pvalue]
  dt[, samplesize := NrSamples]
  dt[, chr := SNPChr]
  dt[, pos := SNPPos]
  dt[, gene_id := Gene]
  dt[, gene_symbol := GeneSymbol]
  dt[, z := Zscore]
  dt[, f_stat := beta^2 / se^2]
  dt[, beta_se_approx_note := "beta and se approximated from Zscore, sample size and effect allele frequency"]

  standard_cols <- c(
    "exposure",
    "exposure_source",
    "SNP",
    "chr",
    "pos",
    "effect_allele",
    "other_allele",
    "eaf",
    "beta",
    "se",
    "pval",
    "samplesize",
    "z",
    "f_stat",
    "gene_id",
    "gene_symbol",
    "GeneChr",
    "GenePos",
    "NrCohorts",
    "FDR",
    "BonferroniP",
    "AF_AlleleA",
    "AF_AlleleB",
    "AF_for_AlleleB",
    "beta_se_approx_note"
  )

  for (cc in standard_cols) {
    if (!cc %in% names(dt)) dt[, (cc) := NA]
  }

  std <- dt[, ..standard_cols]

  out_file <- file.path(out_dir, paste0(output_prefix, ".standardized.tsv.gz"))
  top_file <- file.path(res_dir, paste0(output_prefix, ".standardized.top50.tsv"))
  mr_candidate_file <- file.path(res_dir, paste0(output_prefix, ".mr_candidates_p5e-8.tsv"))

  fwrite(
    std,
    out_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip",
    nThread = threads
  )

  setorder(std, pval)
  fwrite(
    std[seq_len(min(50, .N))],
    top_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    nThread = threads
  )

  mr_candidates <- std[
    !is.na(beta) &
    !is.na(se) &
    !is.na(eaf) &
    !is.na(pval) &
    pval < 5e-8
  ]

  fwrite(
    mr_candidates,
    mr_candidate_file,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    nThread = threads
  )

  qc <- data.table(
    file_key = file_key,
    input_eqtl_file = eqtl_file,
    input_af_file = af_file,
    output_file = out_file,
    top50_file = top_file,
    mr_candidate_file = mr_candidate_file,
    n_rows = nrow(std),
    n_unique_snp = uniqueN(std$SNP),
    n_with_af = sum(!is.na(std$eaf)),
    n_with_beta_se = sum(!is.na(std$beta) & !is.na(std$se)),
    n_p_lt_5e_8 = sum(std$pval < 5e-8, na.rm = TRUE),
    n_p_lt_1e_6 = sum(std$pval < 1e-6, na.rm = TRUE),
    n_p_lt_1e_5 = sum(std$pval < 1e-5, na.rm = TRUE),
    min_pval = suppressWarnings(min(std$pval, na.rm = TRUE)),
    top_snp = std[1]$SNP,
    top_beta = std[1]$beta,
    top_se = std[1]$se,
    top_z = std[1]$z,
    top_pval = std[1]$pval,
    median_f_stat_p5e_8 = suppressWarnings(median(mr_candidates$f_stat, na.rm = TRUE)),
    min_f_stat_p5e_8 = suppressWarnings(min(mr_candidates$f_stat, na.rm = TRUE)),
    max_f_stat_p5e_8 = suppressWarnings(max(mr_candidates$f_stat, na.rm = TRUE))
  )

  log_msg("Standardized file written: ", out_file)
  log_msg("MR candidates written: ", mr_candidate_file)

  rm(eqtl, af, dt, std, mr_candidates)
  gc()

  qc
}

main <- function() {
  log_msg("Start standardizing eQTLGen exposure for gene: ", gene_target)
  log_msg("data.table threads: ", data.table::getDTthreads())

  af_file <- file.path(
    raw_eqtlgen_dir,
    "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"
  )

  files <- list(
    FDR0.05 = file.path(in_dir, "eqtlgen_AHR_FDR0.05_cis_eqtl.extracted.tsv.gz"),
    full = file.path(in_dir, "eqtlgen_AHR_full_cis_eqtl.extracted.tsv.gz")
  )

  qc_list <- list(
    standardize_one(
      "eqtlgen_AHR_FDR0.05_cis_eqtl",
      files$FDR0.05,
      af_file,
      "eqtlgen_AHR_FDR0.05_cis_eqtl"
    ),
    standardize_one(
      "eqtlgen_AHR_full_cis_eqtl",
      files$full,
      af_file,
      "eqtlgen_AHR_full_cis_eqtl"
    )
  )

  qc <- rbindlist(qc_list, fill = TRUE)

  qc_out <- file.path(qc_dir, "eqtlgen_AHR_exposure_standardization_qc.tsv")
  fwrite(qc, qc_out, sep = "\t", quote = FALSE, na = "NA", nThread = threads)

  log_msg("QC written: ", qc_out)
  log_msg("Done.")
}

main()
