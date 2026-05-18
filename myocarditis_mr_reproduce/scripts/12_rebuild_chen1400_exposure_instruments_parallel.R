suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(parallel)
})

ROOT <- normalizePath(Sys.getenv("WENHUAI_MR_ROOT", "."), mustWork = FALSE)

args <- commandArgs(trailingOnly = TRUE)
workers <- ifelse(length(args) >= 1, as.integer(args[1]), 24)
workers <- max(1L, workers)

exposure_dir <- file.path(ROOT, "exposure_full")
acc_list_file <- file.path(ROOT, "results/full1400/accessions/chen1400_accessions.txt")
ld_bfile <- file.path(ROOT, "ld", "EUR")

outdir <- file.path(ROOT, "results/full1400")
raw_dir <- file.path(outdir, "exposure_raw")
inst_dir <- file.path(outdir, "exposure_instruments")
clump_dir <- file.path(outdir, "plink_clump")
log_dir <- file.path(outdir, "logs/exposure")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(inst_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clump_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

find_plink <- function() {
  p <- Sys.which("plink")
  if (!is.na(p) && p != "") return(p)

  if (requireNamespace("genetics.binaRies", quietly = TRUE)) {
    p <- genetics.binaRies::get_plink_binary()
    if (!is.na(p) && file.exists(p)) return(p)
  }

  stop("Cannot find plink.")
}

PLINK_BIN <- find_plink()

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
    stop("No exposure file found for ", gcst)
  }

  files[1]
}

read_exposure_raw <- function(gcst) {
  f <- find_exposure_file(gcst)
  dt <- fread(f, showProgress = FALSE)
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
    id.exposure = paste0("ebi-met1400-", gcst),
    exposure = paste0(gcst, " levels"),
    gcst = gcst
  )]

  x
}

plink_clump_raw <- function(raw, gcst) {
  d <- file.path(clump_dir, gcst)
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
    shQuote(PLINK_BIN),
    shQuote(ld_bfile),
    shQuote(input_file),
    shQuote(plink_out)
  )

  status <- system(cmd, ignore.stdout = TRUE, ignore.stderr = FALSE)

  if (status != 0) {
    stop("PLINK failed for ", gcst, ". Check ", paste0(plink_out, ".log"))
  }

  clumped_file <- paste0(plink_out, ".clumped")
  if (!file.exists(clumped_file)) {
    stop("No clumped file: ", clumped_file)
  }

  clumped <- fread(clumped_file, fill = TRUE, showProgress = FALSE)
  index_snps <- unique(clumped$SNP)

  selected_raw <- raw[SNP %in% index_snps]
  selected_raw[, F_stat := (beta / se)^2]
  selected_raw <- selected_raw[F_stat >= 10]
  setorder(selected_raw, pval)

  selected_raw
}

format_selected_for_mr <- function(selected_raw) {
  if (nrow(selected_raw) == 0) {
    return(data.table())
  }

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

  raw_p_map <- selected_raw[, .(SNP, pval.exposure.raw = pval)]
  fmt_dt <- as.data.table(fmt)
  fmt_dt <- merge(fmt_dt, raw_p_map, by = "SNP", all.x = TRUE, sort = FALSE)

  fmt_dt
}

process_one <- function(gcst) {
  out_file <- file.path(inst_dir, paste0(gcst, "_instruments.tsv"))
  log_file <- file.path(log_dir, paste0(gcst, ".log"))

  if (file.exists(out_file) && file.size(out_file) > 0) {
    return(data.table(gcst = gcst, status = "skip_exists", nsnp = NA_integer_, error = NA_character_))
  }

  zz <- file(log_file, open = "wt")
  sink(zz, type = "output")
  sink(zz, type = "message")

  res <- tryCatch({
    cat("[START]", gcst, "\n")

    raw <- read_exposure_raw(gcst)
    fwrite(raw, file.path(raw_dir, paste0(gcst, "_raw_p5e8.tsv")), sep = "\t")

    cat("[RAW P<5e-8]", nrow(raw), "\n")

    if (nrow(raw) == 0) {
      empty <- data.table()
      fwrite(empty, out_file, sep = "\t")
      list(status = "no_p5e8", nsnp = 0L, error = NA_character_)
    } else {
      selected_raw <- plink_clump_raw(raw, gcst)
      cat("[CLUMPED F>=10]", nrow(selected_raw), "\n")
      cat("[SNPs]", paste(selected_raw$SNP, collapse = ","), "\n")

      fmt <- format_selected_for_mr(selected_raw)
      fwrite(fmt, out_file, sep = "\t")

      list(status = "ok", nsnp = nrow(fmt), error = NA_character_)
    }
  }, error = function(e) {
    cat("[ERROR]", conditionMessage(e), "\n")
    list(status = "error", nsnp = NA_integer_, error = conditionMessage(e))
  })

  sink(type = "message")
  sink(type = "output")
  close(zz)

  data.table(gcst = gcst, status = res$status, nsnp = res$nsnp, error = res$error)
}

acc <- readLines(acc_list_file)
acc <- acc[nzchar(acc)]

message("[WORKERS] ", workers)
message("[N ACC] ", length(acc))
message("[PLINK] ", PLINK_BIN)

status_list <- mclapply(acc, process_one, mc.cores = workers)
status <- rbindlist(status_list, fill = TRUE)

fwrite(status, file.path(outdir, "full1400_exposure_build_status.tsv"), sep = "\t")

ok_files <- list.files(inst_dir, pattern = "_instruments\\.tsv$", full.names = TRUE)

all <- rbindlist(lapply(ok_files, function(f) {
  x <- tryCatch(fread(f, showProgress = FALSE), error = function(e) data.table())
  if (nrow(x) == 0) return(data.table())
  x
}), fill = TRUE)

fwrite(all, file.path(outdir, "full1400_exposure_instruments.tsv"), sep = "\t")

counts <- all[, .(
  exposure_nsnp = uniqueN(SNP),
  SNPs = paste(unique(SNP), collapse = ",")
), by = .(id.exposure, exposure, gcst)]

fwrite(counts, file.path(outdir, "full1400_exposure_nsnp_counts.tsv"), sep = "\t")

message("[DONE]")
message("[STATUS] ", file.path(outdir, "full1400_exposure_build_status.tsv"))
message("[INSTRUMENTS] ", file.path(outdir, "full1400_exposure_instruments.tsv"))
message("[COUNTS] ", file.path(outdir, "full1400_exposure_nsnp_counts.tsv"))

print(status[, .N, by = status])
