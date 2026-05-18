
suppressPackageStartupMessages({
  library(data.table)
})

base <- normalizePath(Sys.getenv("AHR_GWAS_ROOT", "."), mustWork = FALSE)

raw_dir <- file.path(base, "data/raw/outcome")
out_dir <- file.path(base, "data/formatted/outcome")
qc_dir  <- file.path(base, "results/qc")
log_file <- file.path(base, "logs/03_format_outcome_gwas.log")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) {
  msg <- paste0(..., collapse = "")
  line <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

num <- function(x) suppressWarnings(as.numeric(x))

clean_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- sub("^0+([0-9]+)$", "\\1", x)
  x
}

first_rsid <- function(x) {
  x <- as.character(x)
  x[x == ""] <- NA_character_
  x <- sub(",.*$", "", x)
  x
}

write_or_append <- function(dt, outfile, append = FALSE) {
  fwrite(
    dt,
    file = outfile,
    sep = "\t",
    quote = FALSE,
    na = "NA",
    compress = "gzip",
    append = append,
    col.names = !append
  )
}

standard_cols <- c(
  "dataset_key",
  "source_id",
  "open_gwas_id",
  "trait",
  "population",
  "analysis_role",
  "genome_build",
  "chr",
  "pos",
  "SNP",
  "variant_id",
  "variant_id_raw",
  "effect_allele",
  "other_allele",
  "eaf",
  "beta",
  "se",
  "pval",
  "n",
  "ncase",
  "ncontrol",
  "info",
  "direction",
  "af_cases",
  "af_controls",
  "af_allele2_ukb",
  "af_allele2_fg",
  "source_file",
  "inner_file",
  "effect_allele_note"
)

ensure_standard_cols <- function(dt) {
  for (cc in standard_cols) {
    if (!cc %in% names(dt)) {
      dt[, (cc) := NA]
    }
  }
  dt[, ..standard_cols]
}

format_sakaue_bbj_one <- function(zip_path, inner_file) {
  log_msg("Reading BBJ inner file: ", inner_file)

  cmd <- paste("unzip -p", shQuote(zip_path), shQuote(inner_file), "| gzip -dc")
  dt <- fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = TRUE, nThread = max(1, parallel::detectCores() - 1))

  out <- data.table(
    dataset_key = "Sakaue2021_BBJ_Myocarditis",
    source_id = "GCST90018662",
    open_gwas_id = "ebi-a-GCST90018662",
    trait = "Myocarditis",
    population = "East Asian / BBJ",
    analysis_role = "primary_parallel",
    genome_build = "source_original_not_lifted",
    chr = clean_chr(dt$CHR),
    pos = as.integer(dt$POS),
    SNP = fifelse(!is.na(dt$SNPID) & dt$SNPID != "", as.character(dt$SNPID), as.character(dt$v)),
    variant_id = paste0(clean_chr(dt$CHR), ":", dt$POS, ":", dt$Allele1, ":", dt$Allele2),
    variant_id_raw = as.character(dt$v),
    effect_allele = as.character(dt$Allele2),
    other_allele = as.character(dt$Allele1),
    eaf = num(dt$AF_Allele2),
    beta = num(dt$BETA),
    se = num(dt$SE),
    pval = num(dt$p.value),
    n = num(dt$N),
    ncase = 102,
    ncontrol = 177745,
    info = num(dt$imputationInfo),
    direction = NA_character_,
    af_cases = num(dt$AF.Cases),
    af_controls = num(dt$AF.Controls),
    af_allele2_ukb = NA_real_,
    af_allele2_fg = NA_real_,
    source_file = zip_path,
    inner_file = inner_file,
    effect_allele_note = "BETA, SE, p.value and AF are aligned to Allele2"
  )

  ensure_standard_cols(out)
}

format_sakaue_eur_one <- function(zip_path, inner_file) {
  log_msg("Reading EUR inner file: ", inner_file)

  cmd <- paste("unzip -p", shQuote(zip_path), shQuote(inner_file), "| gzip -dc")
  dt <- fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = TRUE, nThread = max(1, parallel::detectCores() - 1))

  af_ukb <- num(dt$AF_Allele2_UKB)
  af_fg <- num(dt$AF_Allele2_FG)
  eaf_used <- fifelse(!is.na(af_fg), af_fg, af_ukb)

  out <- data.table(
    dataset_key = "Sakaue2021_EUR_Myocarditis",
    source_id = "GCST90018882",
    open_gwas_id = "ebi-a-GCST90018882",
    trait = "Myocarditis",
    population = "European",
    analysis_role = "primary_parallel",
    genome_build = "source_original_not_lifted",
    chr = clean_chr(dt$CHR),
    pos = as.integer(dt$POS),
    SNP = paste0(clean_chr(dt$CHR), ":", dt$POS, ":", dt$Allele1, ":", dt$Allele2),
    variant_id = paste0(clean_chr(dt$CHR), ":", dt$POS, ":", dt$Allele1, ":", dt$Allele2),
    variant_id_raw = as.character(dt$v),
    effect_allele = as.character(dt$Allele2),
    other_allele = as.character(dt$Allele1),
    eaf = eaf_used,
    beta = num(dt$BETA),
    se = num(dt$SE),
    pval = num(dt$p.value),
    n = 427911,
    ncase = 633,
    ncontrol = 427278,
    info = NA_real_,
    direction = as.character(dt$Direction),
    af_cases = NA_real_,
    af_controls = NA_real_,
    af_allele2_ukb = af_ukb,
    af_allele2_fg = af_fg,
    source_file = zip_path,
    inner_file = inner_file,
    effect_allele_note = "BETA, SE, p.value and AF are aligned to Allele2; eaf prefers AF_Allele2_FG then AF_Allele2_UKB"
  )

  ensure_standard_cols(out)
}

format_finngen <- function(gz_path) {
  log_msg("Reading FinnGen file: ", gz_path)

  cmd <- paste("zcat", shQuote(gz_path))
  dt <- fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = TRUE, nThread = max(1, parallel::detectCores() - 1))

  if ("#chrom" %in% names(dt)) {
    setnames(dt, "#chrom", "chrom")
  }

  rsid <- first_rsid(dt$rsids)
  variant <- paste0(clean_chr(dt$chrom), ":", dt$pos, ":", dt$ref, ":", dt$alt)

  out <- data.table(
    dataset_key = "FinnGen_R12_I9_MYOCARD",
    source_id = "FINNGEN_R12_I9_MYOCARD",
    open_gwas_id = "finn-b-I9_MYOCARD",
    trait = "Myocarditis",
    population = "Finnish / European",
    analysis_role = "primary_parallel",
    genome_build = "GRCh38_reported_by_FinnGen",
    chr = clean_chr(dt$chrom),
    pos = as.integer(dt$pos),
    SNP = fifelse(!is.na(rsid) & rsid != "", rsid, variant),
    variant_id = variant,
    variant_id_raw = variant,
    effect_allele = as.character(dt$alt),
    other_allele = as.character(dt$ref),
    eaf = num(dt$af_alt),
    beta = num(dt$beta),
    se = num(dt$sebeta),
    pval = num(dt$pval),
    n = NA_real_,
    ncase = NA_real_,
    ncontrol = NA_real_,
    info = NA_real_,
    direction = NA_character_,
    af_cases = num(dt$af_alt_cases),
    af_controls = num(dt$af_alt_controls),
    af_allele2_ukb = NA_real_,
    af_allele2_fg = NA_real_,
    source_file = gz_path,
    inner_file = NA_character_,
    effect_allele_note = "BETA, SE, pval and AF are aligned to alt allele"
  )

  ensure_standard_cols(out)
}

summarise_formatted <- function(infile, dataset_key) {
  log_msg("QC summary for: ", dataset_key)

  dt <- fread(infile, sep = "\t", nThread = max(1, parallel::detectCores() - 1))

  data.table(
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
    n_indel = sum(nchar(dt$effect_allele) != 1 | nchar(dt$other_allele) != 1, na.rm = TRUE),
    n_snv = sum(nchar(dt$effect_allele) == 1 & nchar(dt$other_allele) == 1, na.rm = TRUE)
  )
}

write_column_map <- function() {
  m <- rbindlist(list(
    data.table(dataset_key="Sakaue2021_BBJ_Myocarditis", raw_column="Allele2", standard_column="effect_allele", note="BETA/SE/p.value/AF_Allele2 correspond to Allele2"),
    data.table(dataset_key="Sakaue2021_BBJ_Myocarditis", raw_column="Allele1", standard_column="other_allele", note="Allele1 is the non-effect allele in the standardized file"),
    data.table(dataset_key="Sakaue2021_BBJ_Myocarditis", raw_column="AF_Allele2", standard_column="eaf", note="Effect allele frequency"),
    data.table(dataset_key="Sakaue2021_BBJ_Myocarditis", raw_column="BETA", standard_column="beta", note="Effect of Allele2"),
    data.table(dataset_key="Sakaue2021_BBJ_Myocarditis", raw_column="SE", standard_column="se", note="Standard error of BETA"),
    data.table(dataset_key="Sakaue2021_BBJ_Myocarditis", raw_column="p.value", standard_column="pval", note="Association P value"),

    data.table(dataset_key="Sakaue2021_EUR_Myocarditis", raw_column="Allele2", standard_column="effect_allele", note="BETA/SE/p.value are aligned to Allele2"),
    data.table(dataset_key="Sakaue2021_EUR_Myocarditis", raw_column="Allele1", standard_column="other_allele", note="Allele1 is the non-effect allele in the standardized file"),
    data.table(dataset_key="Sakaue2021_EUR_Myocarditis", raw_column="AF_Allele2_FG / AF_Allele2_UKB", standard_column="eaf", note="Prefer FinnGen AF, otherwise UKB AF"),
    data.table(dataset_key="Sakaue2021_EUR_Myocarditis", raw_column="BETA", standard_column="beta", note="Effect of Allele2"),
    data.table(dataset_key="Sakaue2021_EUR_Myocarditis", raw_column="SE", standard_column="se", note="Standard error of BETA"),
    data.table(dataset_key="Sakaue2021_EUR_Myocarditis", raw_column="p.value", standard_column="pval", note="Association P value"),

    data.table(dataset_key="FinnGen_R12_I9_MYOCARD", raw_column="alt", standard_column="effect_allele", note="FinnGen beta/sebeta/pval are aligned to alt allele"),
    data.table(dataset_key="FinnGen_R12_I9_MYOCARD", raw_column="ref", standard_column="other_allele", note="Reference allele"),
    data.table(dataset_key="FinnGen_R12_I9_MYOCARD", raw_column="af_alt", standard_column="eaf", note="Effect allele frequency"),
    data.table(dataset_key="FinnGen_R12_I9_MYOCARD", raw_column="beta", standard_column="beta", note="Effect of alt allele"),
    data.table(dataset_key="FinnGen_R12_I9_MYOCARD", raw_column="sebeta", standard_column="se", note="Standard error of beta"),
    data.table(dataset_key="FinnGen_R12_I9_MYOCARD", raw_column="pval", standard_column="pval", note="Association P value")
  ))

  outfile <- file.path(qc_dir, "formatted_outcome_gwas_column_map.tsv")
  fwrite(m, outfile, sep = "\t", quote = FALSE, na = "NA")
  log_msg("Column map written: ", outfile)
}

main <- function() {
  if (file.exists(log_file)) file.remove(log_file)

  log_msg("Start formatting outcome GWAS.")

  bbj_zip <- file.path(raw_dir, "sakaue_2021_BBJ/hum0197.v3.BBJ.Myo.v1.zip")
  eur_zip <- file.path(raw_dir, "sakaue_2021_EUR/hum0197.v3.EUR.Myo.v1.zip")
  fg_gz   <- file.path(raw_dir, "finngen_R12/finngen_R12_I9_MYOCARD.gz")

  bbj_out <- file.path(out_dir, "Sakaue2021_BBJ_Myocarditis.formatted.tsv.gz")
  eur_out <- file.path(out_dir, "Sakaue2021_EUR_Myocarditis.formatted.tsv.gz")
  fg_out  <- file.path(out_dir, "FinnGen_R12_I9_MYOCARD.formatted.tsv.gz")

  unlink(c(bbj_out, eur_out, fg_out))

  bbj_inner <- system2("unzip", c("-Z1", shQuote(bbj_zip)), stdout = TRUE)
  bbj_inner <- bbj_inner[grepl("\\.txt\\.gz$", bbj_inner)]

  for (i in seq_along(bbj_inner)) {
    x <- format_sakaue_bbj_one(bbj_zip, bbj_inner[i])
    write_or_append(x, bbj_out, append = i > 1)
    rm(x)
    gc()
  }
  log_msg("Written: ", bbj_out)

  eur_inner <- system2("unzip", c("-Z1", shQuote(eur_zip)), stdout = TRUE)
  eur_inner <- eur_inner[grepl("\\.txt\\.gz$", eur_inner)]

  for (i in seq_along(eur_inner)) {
    x <- format_sakaue_eur_one(eur_zip, eur_inner[i])
    write_or_append(x, eur_out, append = i > 1)
    rm(x)
    gc()
  }
  log_msg("Written: ", eur_out)

  fg <- format_finngen(fg_gz)
  fwrite(fg, fg_out, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
  rm(fg)
  gc()
  log_msg("Written: ", fg_out)

  qc <- rbindlist(list(
    summarise_formatted(bbj_out, "Sakaue2021_BBJ_Myocarditis"),
    summarise_formatted(eur_out, "Sakaue2021_EUR_Myocarditis"),
    summarise_formatted(fg_out, "FinnGen_R12_I9_MYOCARD")
  ), fill = TRUE)

  qc_out <- file.path(qc_dir, "formatted_outcome_gwas_qc.tsv")
  fwrite(qc, qc_out, sep = "\t", quote = FALSE, na = "NA")
  log_msg("QC written: ", qc_out)

  write_column_map()

  log_msg("Done.")
}

main()
