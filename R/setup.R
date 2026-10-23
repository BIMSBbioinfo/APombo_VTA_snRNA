# ============================================================================ #
# setup.R — shared harness for the APombo_Cocaine submission figure repo.
#
# Sourced at the top of every reports/*.Rmd. Provides:
#   - common libraries (light; heavy ones loaded by the Rmd that needs them)
#   - repo root + input resolution (resolve_input / inputs)
#   - the manuscript colour palettes and ggplot2 theme
#
# Self-contained: no dependency on the analysis project, Config.txt, or Mylib.
# Target environment: R 4.5 + Seurat 5.x (QC + MERFISH); the DE and permutation
# reports need only the light libraries below.
# ============================================================================ #

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# ---------------------------------------------------------------------------- #
# Repo root: the parent of this R/ directory. Works whether setup.R is sourced
# from reports/ (knitting) or the repo root (render_all.R).
.repo_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  file_arg <- file_arg[file_arg != "-" & file.exists(file_arg)]   # ignore stdin / non-files
  ofile    <- if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else NULL
  src      <- if (length(file_arg)) file_arg[1] else if (!is.null(ofile)) ofile else NA
  if (!is.na(src)) dirname(dirname(normalizePath(src))) else normalizePath(getwd())
})
# When knitting, the chunk's working dir is reports/; fall back to its parent.
if (!dir.exists(file.path(.repo_root, "data"))) {
  cand <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  if (dir.exists(file.path(cand, "data"))) .repo_root <- cand
}
DATA_DIR    <- file.path(.repo_root, "data")
FIGURE_DIR  <- file.path(.repo_root, "figures")

# ---------------------------------------------------------------------------- #
# resolve_input(pattern): find an input file in data/ by its canonical basename.
# Files may be dropped verbatim with their YYMMDD. date prefix; the newest match
# wins. `pattern` is matched EXACTLY against the basename OR against
# ^<6digits>\.<pattern>$, so e.g. "table-differential_expression_ss_dd.tsv" does
# NOT match the "..._pval0.01_lfc0.25.tsv" threshold variants.
resolve_input <- function(pattern, required = TRUE) {
  files <- list.files(DATA_DIR, full.names = FALSE)
  esc   <- gsub("([.^$*+?()\\[\\]{}|\\\\])", "\\\\\\1", pattern, perl = TRUE)
  hit   <- files[files == pattern | grepl(paste0("^[0-9]{6}\\.", esc, "$"), files)]
  if (!length(hit)) {
    if (required)
      stop(sprintf("Missing required input in data/: '%s'\n  Place the exported file (dated or undated) in: %s",
                   pattern, DATA_DIR), call. = FALSE)
    return(NA_character_)
  }
  file.path(DATA_DIR, sort(hit, decreasing = TRUE)[1])   # newest YYMMDD. wins
}

# Manifest of every input the four reports can consume. Patterns are canonical
# basenames (date prefix optional). See data/README.md for descriptions.
INPUT_PATTERNS <- c(
  # objects (QC, MERFISH)
  all_cells_rds      = "All_Filtered_Cells_Combined.rds",
  dopa_rds           = "Dopaminergic_Cells_Combined.rds",
  # pseudobulk (QC, DE)
  pb_counts_all      = "table-pseudobulk_counts_per_replicate_all.tsv",
  pb_counts_dopa     = "table-pseudobulk_counts_per_replicate_dopa.tsv",
  pb_sheet_all       = "table-pseudobulk_sample_sheet_all.tsv",
  pb_sheet_dopa      = "table-pseudobulk_sample_sheet_dopa.tsv",
  # differential expression (DE)
  de_table           = "table-differential_expression_ss_dd.tsv",
  # gene-set permutation summaries (Permutation)
  perm_all           = "table-gene_set_permutation-all_dopa.tsv",
  perm_subcluster    = "table-gene_set_permutation-per_subcluster.tsv",
  perm_metacluster   = "table-gene_set_permutation-per_metacluster.tsv",
  # gene-set permutation raw draws (Permutation density plots)
  perm_draws_all     = "table-gene_set_permutation_draws-all_dopa.tsv",
  perm_draws_sub     = "table-gene_set_permutation_draws-per_subcluster.tsv",
  perm_draws_meta    = "table-gene_set_permutation_draws-per_metacluster.tsv",
  # metacluster markers (DE)
  markers            = "table-metacluster_markers_wilcox.tsv",
  # MERFISH
  merfish_cells      = "merfish_vta_cells.tsv",
  merfish_mapping    = "merfish_scrna_cluster_mapping.tsv",
  merfish_objects    = "merfish_3d_objects.rds"
)

# inputs(key, required): resolve one manifest entry by key.
inputs <- function(key, required = TRUE) {
  stopifnot(key %in% names(INPUT_PATTERNS))
  resolve_input(INPUT_PATTERNS[[key]], required = required)
}

# Report which inputs are present / missing (used by render_all.R and each Rmd).
input_status <- function(keys = names(INPUT_PATTERNS)) {
  data.frame(
    key     = keys,
    pattern = unname(INPUT_PATTERNS[keys]),
    path    = vapply(keys, function(k) {
      p <- tryCatch(inputs(k, required = FALSE), error = function(e) NA_character_)
      if (is.na(p)) "MISSING" else basename(p)
    }, character(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------- #
# Palettes — copied verbatim from the analysis project (Config.txt
# globals$params$colors) so figures match the manuscript.

# Condition (treatment_time) colours.
CONDITION_COLS <- c(
  sal_d01 = "#269a37", sal_d14 = "#265a37",
  coc_d01 = "#D96E1A", coc_d14 = "#B1529F",
  # reinjection (dd); not in the original ss palette — repo additions:
  sal_24h = "#6BAF74", coc_24h = "#2E7D32"
)

# Contrast colours (cocaine vs saline per timepoint). Reinjection = dark green,
# matching the report's permutation_colors = c(coc_d01, coc_d14, "darkgreen").
CONTRAST_COLS <- c(
  "coc_d01-sal_d01" = "#D96E1A",
  "coc_d14-sal_d14" = "#B1529F",
  "coc_24h-sal_24h" = "darkgreen"
)
# Friendly contrast labels used across reports.
CONTRAST_LABELS <- c(
  "coc_d01-sal_d01" = "Coc d01",
  "coc_d14-sal_d14" = "Coc d14",
  "coc_24h-sal_24h" = "Coc d14 - reinjection"
)

# Up/Down/No colours for volcano/MA.
DIFF_COLS <- c(Up = "#D7191C", Down = "#2C7BB6", No = "grey70")

# Subcluster colours: ggplot2 default hue palette with C10 overridden to
# firebrick (Config.txt globals$params$make_subcluster_colors / overrides).
make_subcluster_colors <- function(cluster_names,
                                    overrides = c(C10 = "firebrick")) {
  cluster_names <- sort(unique(as.character(cluster_names)))
  n <- length(cluster_names)
  hues <- seq(15, 375, length.out = n + 1)[seq_len(n)]
  cols <- grDevices::hcl(h = hues, c = 100, l = 65)
  names(cols) <- cluster_names
  ov <- overrides[names(overrides) %in% cluster_names]
  if (length(ov)) cols[names(ov)] <- unname(ov)
  cols
}
# Canonical 16 DN subclusters C00..C15.
SUBCLUSTER_LEVELS <- sprintf("C%02d", 0:15)
SUBCLUSTER_COLS   <- make_subcluster_colors(SUBCLUSTER_LEVELS)

# Metacluster A/B/C/D and their member subclusters (Config merfish$combinations).
METACLUSTER_MEMBERS <- list(
  A = c("C08","C10","C11","C14","C15"),
  B = c("C09","C13","C06","C07","C02"),
  C = c("C01","C00"),
  D = c("C03","C04")
)
# Qualitative A-D palette (repo choice — the analysis report has no fixed A-D map).
METACLUSTER_COLS <- c(A = "#1B9E77", B = "#D95F02", C = "#7570B3", D = "#E7298A")

# Manuscript x-axis order for the permutation subcluster heatmaps (from the report).
SUBCLUSTER_ORDER <- c("C13","C09","C02","C06","C07","C01","C00","C03",
                      "C04","C08","C15","C10","C14","C11","C05","C12")

# ---------------------------------------------------------------------------- #
# Shared ggplot2 theme (matches the analysis project's globals$theme spirit).
theme_submission <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
      axis.line  = element_line(colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = NA),
      legend.key = element_blank()
    )
}
theme_set(theme_submission())

# Helper: save a ggplot to figures/<topic>/<name>.pdf (+ create the dir).
save_fig <- function(plot, topic, name, width = 7, height = 5, dev = "pdf") {
  d <- file.path(FIGURE_DIR, topic)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(d, paste0(name, ".", dev)), plot, width = width, height = height)
}

message("setup.R loaded — repo root: ", .repo_root)
