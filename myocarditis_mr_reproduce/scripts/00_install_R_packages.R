options(timeout = 1000)

# MRC IEU 官方 r-universe + CRAN
options(repos = c(
  mrcieu = "https://mrcieu.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
))

base_pkgs <- c(
  "data.table", "dplyr", "stringr", "tibble", "readr",
  "ggplot2", "remotes", "devtools", "meta", "metafor",
  "ieugwasr", "genetics.binaRies"
)

for (p in base_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("[INSTALL] ", p)
    install.packages(p)
  } else {
    message("[OK] ", p)
  }
}

# TwoSampleMR 不在 CRAN，直接从 MRC IEU r-universe 装
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  message("[INSTALL] TwoSampleMR from MRC IEU r-universe")
  install.packages("TwoSampleMR", repos = c(
    "https://mrcieu.r-universe.dev",
    "https://cloud.r-project.org"
  ))
} else {
  message("[OK] TwoSampleMR")
}

# MendelianRandomization 在 CRAN，有版本要求但先装通
if (!requireNamespace("MendelianRandomization", quietly = TRUE)) {
  message("[INSTALL] MendelianRandomization")
  install.packages("MendelianRandomization")
} else {
  message("[OK] MendelianRandomization")
}

# coloc 可以先尝试指定版本；失败就装最新版
if (!requireNamespace("coloc", quietly = TRUE)) {
  message("[INSTALL] coloc")
  tryCatch({
    remotes::install_version("coloc", version = "5.2.3", upgrade = "never")
  }, error = function(e) {
    message("[WARN] coloc 5.2.3 failed, install latest coloc")
    install.packages("coloc")
  })
} else {
  message("[OK] coloc")
}

# MRPRESSO 从 GitHub 装
if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  message("[INSTALL] MRPRESSO from GitHub")
  devtools::install_github("rondolab/MR-PRESSO", upgrade = "never")
} else {
  message("[OK] MRPRESSO")
}

cat("\n[DONE] package versions:\n")
for (p in c("TwoSampleMR", "MendelianRandomization", "MRPRESSO", "coloc", "ieugwasr", "genetics.binaRies")) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(p, as.character(utils::packageVersion(p)), "\n")
  } else {
    cat(p, "NOT_INSTALLED\n")
  }
}
