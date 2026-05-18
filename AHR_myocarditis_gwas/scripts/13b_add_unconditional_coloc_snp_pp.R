
suppressPackageStartupMessages({
  library(data.table)
})

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

summary_file <- file.path(base, "results/coloc/AHR_eqtlgen_coloc_main_result.tsv")
top_file <- file.path(base, "results/coloc/AHR_eqtlgen_coloc_top_snp_posterior.tsv")
out_file <- file.path(base, "results/coloc/AHR_eqtlgen_coloc_top_snp_posterior.with_unconditional.tsv")
log_file <- file.path(base, "logs/13b_add_unconditional_coloc_snp_pp.log")

if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

main <- function() {
  log_msg("Start adding unconditional SNP H4 posterior.")

  s <- fread(summary_file, sep = "\t")
  x <- fread(top_file, sep = "\t")

  s <- s[, .(
    outcome_dataset,
    PP.H4.abf = as.numeric(PP.H4.abf),
    PP.H3.abf = as.numeric(PP.H3.abf),
    PP.H1.abf = as.numeric(PP.H1.abf)
  )]

  x <- merge(x, s, by = "outcome_dataset", all.x = TRUE, sort = FALSE)

  if (!"SNP.PP.H4" %in% names(x)) {
    stop("SNP.PP.H4 column not found in top SNP posterior table.")
  }

  x[, SNP.PP.H4 := as.numeric(SNP.PP.H4)]
  x[, unconditional_SNP_H4 := PP.H4.abf * SNP.PP.H4]

  setorder(x, outcome_dataset, -unconditional_SNP_H4)

  fwrite(x, out_file, sep = "\t", quote = FALSE, na = "NA")

  log_msg("Written: ", out_file)
  log_msg("Done.")
}

main()
