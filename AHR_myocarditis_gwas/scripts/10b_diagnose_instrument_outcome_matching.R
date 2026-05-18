
suppressPackageStartupMessages({
  library(data.table)
})

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

instrument_dir <- file.path(base, "data/instruments/exposure")
outcome_locus_dir <- file.path(base, "data/locus/outcome")
qc_dir <- file.path(base, "results/qc")
log_file <- file.path(base, "logs/10b_diagnose_instrument_outcome_matching.log")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

read_gz <- function(path) {
  fread(cmd = paste("zcat", shQuote(path)), sep = "\t", header = TRUE)
}

instrument_files <- list.files(
  instrument_dir,
  pattern = "^eqtlgen_AHR_full_cis_eqtl\\..*\\.clumped\\.tsv$",
  full.names = TRUE
)

outcome_files <- c(
  Sakaue2021_BBJ_Myocarditis = file.path(outcome_locus_dir, "Sakaue2021_BBJ_Myocarditis.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz"),
  Sakaue2021_EUR_Myocarditis = file.path(outcome_locus_dir, "Sakaue2021_EUR_Myocarditis.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz"),
  FinnGen_R12_I9_MYOCARD = file.path(outcome_locus_dir, "FinnGen_R12_I9_MYOCARD.AHR_locus.chr7_16200000_18500000.dedup_by_variant.tsv.gz")
)

outcomes <- lapply(outcome_files, function(f) {
  x <- read_gz(f)
  x[, outcome_SNP := as.character(SNP)]
  x[, chr := as.character(chr)]
  x[, pos := as.integer(pos)]
  x
})

res <- list()
k <- 1

for (inst_file in instrument_files) {
  inst <- fread(inst_file, sep = "\t")
  setting <- basename(inst_file)
  setting <- sub("^eqtlgen_AHR_full_cis_eqtl\\.", "", setting)
  setting <- sub("\\.clumped\\.tsv$", "", setting)

  inst[, exp_SNP := as.character(SNP)]
  inst[, chr := as.character(chr)]
  inst[, pos := as.integer(pos)]

  for (oname in names(outcomes)) {
    out <- outcomes[[oname]]

    rs_match <- merge(
      inst[, .(exp_SNP)],
      out[, .(outcome_SNP)],
      by.x = "exp_SNP",
      by.y = "outcome_SNP"
    )

    pos_match <- merge(
      inst[, .(chr, pos, exp_SNP)],
      out[, .(chr, pos, outcome_SNP)],
      by = c("chr", "pos"),
      allow.cartesian = TRUE
    )

    res[[k]] <- data.table(
      instrument_setting = setting,
      outcome_dataset = oname,
      n_instruments = nrow(inst),
      n_rsid_matched = uniqueN(rs_match$exp_SNP),
      n_position_matched = uniqueN(pos_match$exp_SNP),
      rsid_matched_snps = paste(sort(unique(rs_match$exp_SNP)), collapse = ";"),
      position_matched_snps = paste(sort(unique(pos_match$exp_SNP)), collapse = ";")
    )
    k <- k + 1
  }
}

ans <- rbindlist(res, fill = TRUE)
out_file <- file.path(qc_dir, "AHR_instrument_outcome_matching_diagnosis.tsv")
fwrite(ans, out_file, sep = "\t", quote = FALSE, na = "NA")

log_msg("Written: ", out_file)
