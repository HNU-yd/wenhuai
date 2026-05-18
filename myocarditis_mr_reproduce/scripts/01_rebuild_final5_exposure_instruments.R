suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

ROOT <- normalizePath(Sys.getenv("WENHUAI_MR_ROOT", "."), mustWork = FALSE)

exposure_dir <- file.path(ROOT, "exposure")
ld_bfile <- file.path(ROOT, "ld", "EUR")
outdir <- file.path(ROOT, "results", "clean_rebuild_exposure_v2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

targets <- data.table(
  id.exposure = c(
    "ebi-met1400-GCST90199636",
    "ebi-met1400-GCST90199772",
    "ebi-met1400-GCST90199813",
    "ebi-met1400-GCST90200661",
    "ebi-met1400-GCST90200680"
  ),
  gcst = c(
    "GCST90199636",
    "GCST90199772",
    "GCST90199813",
    "GCST90200661",
    "GCST90200680"
  ),
  exposure = c(
    "Kynurenine levels",
    "1-stearoyl-GPE (18:0) levels",
    "Deoxycarnitine levels",
    "X-25422 levels",
    "5-acetylamino-6-formylamino-3-methyluracil levels"
  ),
  expected_nsnp_paper = c(5, 7, 5, 6, 4)
)

find_plink <- function() {
  p <- Sys.which("plink")
  if (!is.na(p) && p != "") return(p)

  if (requireNamespace("genetics.binaRies", quietly = TRUE)) {
    p <- genetics.binaRies::get_plink_binary()
    if (!is.na(p) && file.exists(p)) return(p)
  }

  stop("Cannot find plink.")
}

plink_bin <- find_plink()
cat("[PLINK]", plink_bin, "\n")

pick_col <- function(nms, candidates, required = TRUE) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[1])

  hit2 <- nms[tolower(nms) %in% tolower(candidates)]
  if (length(hit2) > 0) return(hit2[1])

  if (required) {
    stop("Cannot find column among: ", paste(candidates, collapse = ", "),
         "\nAvailable: ", paste(nms, collapse = ", "))
  }

  NA_character_
}

find_exposure_file <- function(gcst) {
  d <- file.path(exposure_dir, gcst)
  files <- list.files(d, recursive = TRUE, full.names = TRUE)

  files <- files[grepl(paste0(gcst, ".*\\.tsv(\\.gz)?$"), basename(files))]
  files <- files[!grepl("\\.tbi$", files)]
  files <- files[!grepl("robots|index|html|README|md5", basename(files), ignore.case = TRUE)]
  files <- files[order(!grepl("buildGRCh38", basename(files)), basename(files))]

  if (length(files) == 0) {
    stop("No exposure file found for ", gcst,
         "\nCurrent files: ", paste(list.files(d, recursive = TRUE), collapse = ", "))
  }

  files[1]
}

read_exposure_raw <- function(gcst, id.exposure, exposure_name) {
  f <- find_exposure_file(gcst)
  cat("\n========== READ", gcst, exposure_name, "==========\n")
  cat("[FILE]", f, "\n")

  dt <- fread(f, showProgress = TRUE)
  nms <- names(dt)

  snp_col <- pick_col(nms, c("hm_rsid", "rsid", "SNP", "variant_id", "rs_id"))
  beta_col <- pick_col(nms, c("hm_beta", "beta", "BETA", "effect"))
  se_col <- pick_col(nms, c("standard_error", "se", "SE", "sebeta"))
  p_col <- pick_col(nms, c("p_value", "pval", "p", "P"))
  ea_col <- pick_col(nms, c("hm_effect_allele", "effect_allele", "EA", "ALT", "alt"))
  oa_col <- pick_col(nms, c("hm_other_allele", "other_allele", "NEA", "REF", "ref"))

  eaf_col <- pick_col(nms, c("effect_allele_frequency", "hm_effect_allele_frequency", "eaf", "EAF", "af"), required = FALSE)
  chr_col <- pick_col(nms, c("hm_chrom", "chromosome", "chr", "CHR", "chrom"), required = FALSE)
  pos_col <- pick_col(nms, c("hm_pos", "base_pair_location", "pos", "POS", "position"), required = FALSE)

  x <- data.table(
    SNP = as.character(dt[[snp_col]]),
    beta = as.numeric(dt[[beta_col]]),
    se = as.numeric(dt[[se_col]]),
    pval = as.numeric(dt[[p_col]]),
    effect_allele = toupper(as.character(dt[[ea_col]])),
    other_allele = toupper(as.character(dt[[oa_col]]))
  )

  if (!is.na(eaf_col)) x[, eaf := as.numeric(dt[[eaf_col]])]
  if (!is.na(chr_col)) x[, chr_name := as.character(dt[[chr_col]])]
  if (!is.na(pos_col)) x[, chrom_start := as.integer(dt[[pos_col]])]

  x <- x[!is.na(SNP) & SNP != "" & grepl("^rs", SNP)]
  x <- x[!is.na(beta) & !is.na(se) & !is.na(pval)]
  x <- x[pval < 5e-8]
  x[, F_stat_raw := (beta / se)^2]
  setorder(x, pval)
  x <- unique(x, by = "SNP")

  x[, `:=`(
    id.exposure = id.exposure,
    exposure = exposure_name,
    gcst = gcst
  )]

  cat("[RAW P < 5e-8]", nrow(x), "\n")
  fwrite(x, file.path(outdir, paste0(gcst, "_raw_p5e8.tsv")), sep = "\t")

  x
}

plink_clump_raw <- function(raw, gcst) {
  d <- file.path(outdir, "plink_clump", gcst)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

  input <- raw[, .(SNP, P = pval)]
  input <- input[!is.na(SNP) & !is.na(P)]
  input <- unique(input, by = "SNP")
  setorder(input, P)

  input_file <- file.path(d, paste0(gcst, "_rawP_clump_input.tsv"))
  fwrite(input, input_file, sep = "\t")

  plink_out <- file.path(d, paste0(gcst, "_rawP_clump"))

  cmd <- sprintf(
    "%s --bfile %s --clump %s --clump-kb 10000 --clump-r2 0.001 --clump-p1 5e-8 --clump-p2 5e-8 --out %s",
    shQuote(plink_bin),
    shQuote(ld_bfile),
    shQuote(input_file),
    shQuote(plink_out)
  )

  cat("[CMD]", cmd, "\n")
  status <- system(cmd)

  if (status != 0) {
    stop("PLINK failed for ", gcst, ". Check ", paste0(plink_out, ".log"))
  }

  clumped_file <- paste0(plink_out, ".clumped")
  if (!file.exists(clumped_file)) {
    stop("No clumped file: ", clumped_file)
  }

  clumped <- fread(clumped_file, fill = TRUE)
  index_snps <- unique(clumped$SNP)

  selected_raw <- raw[SNP %in% index_snps]
  selected_raw[, F_stat := (beta / se)^2]
  selected_raw <- selected_raw[F_stat >= 10]
  setorder(selected_raw, pval)

  fwrite(selected_raw, file.path(d, paste0(gcst, "_selected_rawP_F_ge_10.tsv")), sep = "\t")

  cat("[SELECTED RAW]", gcst, nrow(selected_raw), paste(selected_raw$SNP, collapse = ", "), "\n")

  selected_raw
}

format_selected_for_mr <- function(selected_raw) {
  fmt <- format_data(
    as.data.frame(selected_raw),
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = if ("eaf" %in% names(selected_raw)) "eaf" else NULL,
    pval_col = "pval",
    chr_col = if ("chr_name" %in% names(selected_raw)) "chr_name" else NULL,
    pos_col = if ("chrom_start" %in% names(selected_raw)) "chrom_start" else NULL,
    min_pval = 1e-320
  )

  fmt$id.exposure <- unique(selected_raw$id.exposure)
  fmt$exposure <- unique(selected_raw$exposure)
  fmt$gcst <- unique(selected_raw$gcst)
  fmt$F_stat <- (fmt$beta.exposure / fmt$se.exposure)^2

  # 保留原始 p 值，方便审计。TwoSampleMR 的 pval.exposure 可能会受 min_pval 影响。
  raw_p_map <- selected_raw[, .(SNP, pval.exposure.raw = pval)]
  fmt_dt <- as.data.table(fmt)
  fmt_dt <- merge(fmt_dt, raw_p_map, by = "SNP", all.x = TRUE, sort = FALSE)

  fmt_dt
}

all <- list()

for (i in seq_len(nrow(targets))) {
  raw <- read_exposure_raw(
    gcst = targets$gcst[i],
    id.exposure = targets$id.exposure[i],
    exposure_name = targets$exposure[i]
  )

  selected_raw <- plink_clump_raw(raw, targets$gcst[i])
  selected_fmt <- format_selected_for_mr(selected_raw)

  all[[i]] <- selected_fmt
}

exp_dat <- rbindlist(all, fill = TRUE)

counts <- exp_dat[, .(
  exposure_nsnp = uniqueN(SNP),
  SNPs = paste(SNP, collapse = ",")
), by = .(id.exposure, exposure)]

counts <- merge(
  counts,
  targets[, .(id.exposure, expected_nsnp_paper)],
  by = "id.exposure",
  all.x = TRUE
)

cat("\n========== FINAL COUNTS ==========\n")
print(counts)

fwrite(exp_dat, file.path(outdir, "final5_exposure_instruments.clean_v2.tsv"), sep = "\t")
fwrite(counts, file.path(outdir, "final5_exposure_counts.clean_v2.tsv"), sep = "\t")

# 覆盖主 results 文件，给 02 脚本继续使用
fwrite(exp_dat, file.path(ROOT, "results", "final5_exposure_instruments.tsv"), sep = "\t")
fwrite(counts, file.path(ROOT, "results", "final5_exposure_nsnp_counts.tsv"), sep = "\t")

cat("\n[DONE]\n")
cat("[OUT]", file.path(ROOT, "results", "final5_exposure_instruments.tsv"), "\n")
