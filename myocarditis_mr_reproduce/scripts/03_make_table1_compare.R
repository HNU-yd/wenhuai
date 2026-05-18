suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- normalizePath(Sys.getenv("WENHUAI_MR_ROOT", "."), mustWork = FALSE)

ivw <- fread(file.path(ROOT, "results/clean_outcome_mr/11_final5_ivw_only.tsv"))
het <- fread(file.path(ROOT, "results/clean_outcome_mr/12_heterogeneity.tsv"))
pleio <- fread(file.path(ROOT, "results/clean_outcome_mr/13_pleiotropy_egger_intercept.tsv"))

paper <- data.table(
  id.exposure = c(
    "ebi-met1400-GCST90199636",
    "ebi-met1400-GCST90199772",
    "ebi-met1400-GCST90199813",
    "ebi-met1400-GCST90200661",
    "ebi-met1400-GCST90200680"
  ),
  paper_exposure = c(
    "Kynurenine",
    "1-stearoyl-GPE (18:0)",
    "Deoxycarnitine",
    "X-25,422",
    "5-acetylamino-6-formylamino-3-methyluracil"
  ),
  paper_nsnp = c(5, 7, 5, 6, 4),
  paper_or = c(1.441, 1.263, 0.813, 0.721, 0.864),
  paper_p = c(0.018, 0.029, 0.029, 0.018, 0.009)
)

x <- merge(
  ivw[, .(
    id.exposure,
    exposure,
    nsnp,
    b,
    se,
    pval,
    or,
    or_lci95,
    or_uci95,
    p_fdr_ivw_final5
  )],
  paper,
  by = "id.exposure",
  all.x = TRUE
)

x[, or_diff := or - paper_or]
x[, nsnp_match := nsnp == paper_nsnp]

fwrite(
  x,
  file.path(ROOT, "results/final_release/table1_compare_finngen.tsv"),
  sep = "\t"
)

print(x)
