suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

ROOT <- normalizePath(Sys.getenv("WENHUAI_MR_ROOT", "."), mustWork = FALSE)

exp_file <- file.path(ROOT, "results", "final5_exposure_instruments.tsv")
outcome_file <- file.path(ROOT, "outcome", "finngen_R10_I9_MYOCARD.gz")
outdir <- file.path(ROOT, "results", "clean_outcome_mr")
plotdir <- file.path(ROOT, "plots", "clean_outcome_mr")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

targets <- data.table(
  id.exposure = c(
    "ebi-met1400-GCST90199636",
    "ebi-met1400-GCST90199772",
    "ebi-met1400-GCST90199813",
    "ebi-met1400-GCST90200661",
    "ebi-met1400-GCST90200680"
  ),
  expected_nsnp_paper = c(5, 7, 5, 6, 4)
)

message("========== READ EXPOSURE INSTRUMENTS ==========")

if (!file.exists(exp_file)) {
  stop("Missing exposure instrument file: ", exp_file)
}

exp_dat <- fread(exp_file)
exp_dat <- as.data.table(exp_dat)

required_cols <- c(
  "SNP",
  "id.exposure",
  "exposure",
  "beta.exposure",
  "se.exposure",
  "effect_allele.exposure",
  "other_allele.exposure",
  "pval.exposure",
  "chr.exposure",
  "pos.exposure"
)

miss_cols <- setdiff(required_cols, names(exp_dat))
if (length(miss_cols) > 0) {
  stop("Exposure file missing columns: ", paste(miss_cols, collapse = ", "))
}

exp_dat[, chr_exp := as.character(chr.exposure)]
exp_dat[, pos_exp := as.integer(pos.exposure)]
exp_dat[, ea_exp := toupper(effect_allele.exposure)]
exp_dat[, oa_exp := toupper(other_allele.exposure)]
exp_dat[, locus := paste(chr_exp, pos_exp, sep = ":")]

target_snps <- unique(exp_dat$SNP)
target_loci <- unique(exp_dat[!is.na(chr_exp) & !is.na(pos_exp), locus])

message("[TARGET SNP N] ", length(target_snps))
message("[TARGET LOCUS N] ", length(target_loci))

exp_counts <- exp_dat[, .(exposure_nsnp = uniqueN(SNP)), by = .(id.exposure, exposure)]
exp_counts <- merge(exp_counts, targets, by = "id.exposure", all.x = TRUE)
fwrite(exp_counts, file.path(outdir, "01_exposure_nsnp_counts.tsv"), sep = "\t")
print(exp_counts)

message("========== READ FINNGEN OUTCOME ==========")

out_raw <- fread(
  cmd = sprintf("zcat %s", shQuote(outcome_file)),
  select = c("#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"),
  showProgress = TRUE
)

setnames(out_raw, old = "#chrom", new = "chrom", skip_absent = TRUE)

out_raw[, chrom := as.character(chrom)]
out_raw[, pos := as.integer(pos)]
out_raw[, ref := toupper(ref)]
out_raw[, alt := toupper(alt)]
out_raw[, locus := paste(chrom, pos, sep = ":")]

message("========== MATCH OUTCOME BY LOCUS + ALLELE ==========")

exp_locus <- unique(exp_dat[, .(
  SNP,
  id.exposure,
  exposure,
  chrom = chr_exp,
  pos = pos_exp,
  ea_exp,
  oa_exp,
  locus
)])

out_locus_hit <- out_raw[locus %in% target_loci]

locus_merged <- merge(
  exp_locus,
  out_locus_hit,
  by = c("chrom", "pos"),
  allow.cartesian = TRUE
)

locus_allele_hit <- locus_merged[
  (ea_exp == alt & oa_exp == ref) |
  (ea_exp == ref & oa_exp == alt)
]

if (nrow(locus_allele_hit) > 0) {
  locus_allele_hit <- locus_allele_hit[, .(
    SNP,
    chrom,
    pos,
    ref,
    alt,
    pval,
    beta,
    sebeta,
    af_alt,
    rsids,
    locus = paste(chrom, pos, sep = ":"),
    match_source = "locus_allele"
  )]
} else {
  locus_allele_hit <- data.table()
}

message("[LOCUS + ALLELE HIT SNP N] ", uniqueN(locus_allele_hit$SNP))

message("========== MATCH OUTCOME BY RSID FALLBACK ==========")

rs_pattern <- paste0("(^|,)(", paste(target_snps, collapse = "|"), ")(,|$)")
out_rsid_hit <- out_raw[!is.na(rsids) & grepl(rs_pattern, rsids)]

if (nrow(out_rsid_hit) > 0) {
  rsid_hit <- out_rsid_hit[, {
    ids <- unlist(strsplit(rsids, ",", fixed = TRUE))
    ids <- intersect(ids, target_snps)
    if (length(ids) == 0) {
      NULL
    } else {
      .(SNP = ids)
    }
  }, by = .(chrom, pos, ref, alt, pval, beta, sebeta, af_alt, rsids, locus)]

  if (nrow(rsid_hit) > 0) {
    rsid_hit[, match_source := "rsid"]
  }
} else {
  rsid_hit <- data.table()
}

message("[RSID HIT SNP N] ", uniqueN(rsid_hit$SNP))

message("========== COMBINE OUTCOME HITS ==========")

outcome_hits <- rbindlist(
  list(
    locus_allele_hit,
    rsid_hit[, .(SNP, chrom, pos, ref, alt, pval, beta, sebeta, af_alt, rsids, locus, match_source)]
  ),
  fill = TRUE
)

if (nrow(outcome_hits) == 0) {
  stop("No outcome hits found.")
}

# 关键：优先选择 allele-compatible 的 locus_allele，而不是 rsid
outcome_hits[, match_rank := fifelse(match_source == "locus_allele", 1L, 2L)]
setorder(outcome_hits, SNP, match_rank)
outcome_hits <- unique(outcome_hits, by = "SNP")

match_status <- data.table(SNP = target_snps)
match_status[, matched := SNP %in% outcome_hits$SNP]
match_status <- merge(
  match_status,
  unique(outcome_hits[, .(SNP, match_source, chrom, pos, ref, alt, rsids)]),
  by = "SNP",
  all.x = TRUE
)

fwrite(match_status, file.path(outdir, "02_outcome_match_status.tsv"), sep = "\t")
fwrite(outcome_hits, file.path(outdir, "03_finngen_outcome_hits.tsv"), sep = "\t")

missing <- match_status[matched == FALSE, SNP]
if (length(missing) > 0) {
  fwrite(data.table(SNP = missing), file.path(outdir, "04_finngen_missing_outcome_snps.tsv"), sep = "\t")
  message("[MISSING OUTCOME SNPs] ", paste(missing, collapse = ", "))
} else {
  message("[OUTCOME] all exposure SNPs matched.")
}

message("========== FORMAT OUTCOME ==========")

out_dat <- format_data(
  as.data.frame(outcome_hits),
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval",
  chr_col = "chrom",
  pos_col = "pos"
)

out_dat$outcome <- "Myocarditis"
out_dat$id.outcome <- "finngen_R10_I9_MYOCARD"

fwrite(as.data.table(out_dat), file.path(outdir, "05_finngen_outcome_formatted.tsv"), sep = "\t")

message("========== HARMONISE ==========")

exp_df <- as.data.frame(exp_dat)
dat <- harmonise_data(exp_df, out_dat, action = 2)
dat_dt <- as.data.table(dat)

fwrite(dat_dt, file.path(outdir, "06_harmonised_all.tsv"), sep = "\t")

if ("mr_keep" %in% names(dat_dt)) {
  dat_keep <- dat_dt[mr_keep == TRUE]
} else {
  dat_keep <- dat_dt
}

fwrite(dat_keep, file.path(outdir, "07_harmonised_mr_keep.tsv"), sep = "\t")

harm_all_counts <- dat_dt[, .(harmonised_all_nsnp = uniqueN(SNP)), by = .(id.exposure, exposure)]
harm_keep_counts <- dat_keep[, .(harmonised_mr_keep_nsnp = uniqueN(SNP)), by = .(id.exposure, exposure)]

counts <- merge(exp_counts, harm_all_counts, by = c("id.exposure", "exposure"), all.x = TRUE)
counts <- merge(counts, harm_keep_counts, by = c("id.exposure", "exposure"), all.x = TRUE)
counts <- merge(counts, targets, by = "id.exposure", all.x = TRUE, suffixes = c("", ".target"))

if ("expected_nsnp_paper.target" %in% names(counts)) {
  counts[, expected_nsnp_paper := fifelse(
    is.na(expected_nsnp_paper),
    expected_nsnp_paper.target,
    expected_nsnp_paper
  )]
  counts[, expected_nsnp_paper.target := NULL]
}

fwrite(counts, file.path(outdir, "08_nsnp_counts_all_steps.tsv"), sep = "\t")
print(counts)

message("========== GCST90200680 DETAIL ==========")
print(dat_dt[id.exposure == "ebi-met1400-GCST90200680", .(
  SNP,
  exposure,
  effect_allele.exposure,
  other_allele.exposure,
  effect_allele.outcome,
  other_allele.outcome,
  beta.exposure,
  beta.outcome,
  eaf.exposure,
  eaf.outcome,
  mr_keep,
  remove,
  palindromic,
  ambiguous
)])

fwrite(
  dat_dt[id.exposure == "ebi-met1400-GCST90200680"],
  file.path(outdir, "09_GCST90200680_harmonised_detail.tsv"),
  sep = "\t"
)

message("========== MR ==========")

dat_keep_df <- as.data.frame(dat_keep)

mr_res <- mr(
  dat_keep_df,
  method_list = c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median",
    "mr_weighted_mode"
  )
)

mr_or <- generate_odds_ratios(mr_res)

ivw_idx <- mr_or$method == "Inverse variance weighted"
mr_or$p_fdr_ivw_final5 <- NA_real_
mr_or$p_fdr_ivw_final5[ivw_idx] <- p.adjust(mr_or$pval[ivw_idx], method = "fdr")

fwrite(as.data.table(mr_or), file.path(outdir, "10_final5_mr_results.tsv"), sep = "\t")

ivw <- as.data.table(mr_or)[method == "Inverse variance weighted"]
ivw <- merge(ivw, targets, by = "id.exposure", all.x = TRUE)
fwrite(ivw, file.path(outdir, "11_final5_ivw_only.tsv"), sep = "\t")

print(ivw[, .(
  id.exposure,
  exposure,
  nsnp,
  expected_nsnp_paper,
  b,
  se,
  pval,
  or,
  or_lci95,
  or_uci95,
  p_fdr_ivw_final5
)])

message("========== SENSITIVITY ==========")

het <- tryCatch(mr_heterogeneity(dat_keep_df), error = function(e) {
  message("[mr_heterogeneity ERROR] ", conditionMessage(e))
  NULL
})
if (!is.null(het)) fwrite(as.data.table(het), file.path(outdir, "12_heterogeneity.tsv"), sep = "\t")

pleio <- tryCatch(mr_pleiotropy_test(dat_keep_df), error = function(e) {
  message("[mr_pleiotropy ERROR] ", conditionMessage(e))
  NULL
})
if (!is.null(pleio)) fwrite(as.data.table(pleio), file.path(outdir, "13_pleiotropy_egger_intercept.tsv"), sep = "\t")

single <- tryCatch(mr_singlesnp(dat_keep_df), error = function(e) {
  message("[mr_singlesnp ERROR] ", conditionMessage(e))
  NULL
})
if (!is.null(single)) fwrite(as.data.table(single), file.path(outdir, "14_single_snp.tsv"), sep = "\t")

loo <- tryCatch(mr_leaveoneout(dat_keep_df), error = function(e) {
  message("[mr_leaveoneout ERROR] ", conditionMessage(e))
  NULL
})
if (!is.null(loo)) fwrite(as.data.table(loo), file.path(outdir, "15_leave_one_out.tsv"), sep = "\t")

message("========== DONE ==========")
message("[COUNTS] ", file.path(outdir, "08_nsnp_counts_all_steps.tsv"))
message("[IVW]    ", file.path(outdir, "11_final5_ivw_only.tsv"))
message("[DETAIL] ", file.path(outdir, "09_GCST90200680_harmonised_detail.tsv"))
