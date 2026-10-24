#!/usr/bin/env Rscript
# install.R — install the full environment for the submission figure reports.
# Usage:  Rscript install.R
#
# Installs every package the four reports need (CRAN + Bioconductor) in one go.
# Idempotent: already-installed packages are skipped. Uses BiocManager, which
# resolves BOTH CRAN and Bioconductor packages against a single Bioconductor
# release, so versions stay mutually compatible.
#
# For a byte-for-byte reproducible lockfile, see the renv note at the bottom.

# --- R version gate -------------------------------------------------------- #
if (getRversion() < "4.5.0")
  stop("R >= 4.5 required (Seurat 5.x target); found ", getRversion())

# --- package set (mirrors DESCRIPTION Imports) ----------------------------- #
pkgs <- c(
  # CRAN — table / plotting / reporting stack
  "data.table", "dplyr", "tidyr", "stringr", "ggplot2",
  "rmarkdown", "knitr",
  # CRAN — single-cell (01 QC) + interactive 3D (04 MERFISH)
  "Seurat", "SeuratObject", "plotly", "htmlwidgets",
  # Bioconductor — heatmaps (02 DE / 03 Permutation)
  "ComplexHeatmap"
)

# --- bootstrap BiocManager, then install the lot --------------------------- #
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) {
  message("Installing: ", paste(need, collapse = ", "))
  BiocManager::install(need, update = FALSE, ask = FALSE)
} else {
  message("All packages already installed.")
}

# --- verify ---------------------------------------------------------------- #
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
if (all(ok)) {
  message("\nEnvironment ready. Seurat ",
          as.character(utils::packageVersion("Seurat")),
          " on R ", getRversion(), ".")
} else {
  stop("Failed to install: ", paste(pkgs[!ok], collapse = ", "))
}

# --- optional: freeze an exact, reproducible lockfile ---------------------- #
# Run these in THIS environment to capture exact versions into renv.lock:
#   install.packages("renv")
#   renv::init(bare = TRUE)          # adopt the repo as an renv project
#   renv::snapshot()                 # writes renv.lock from installed versions
# Thereafter anyone can reproduce with:  renv::restore()
