# ============================================================================ #
# plot_helpers.R — self-contained figure functions for the submission reports.
#
# Trimmed, dependency-free equivalents of the analysis project's Plot_* helpers
# (no source.pro / cacheFile / Mylib). Sourced after R/setup.R, so palettes,
# theme, CONTRAST_* and *_COLS are already in scope.
#
# Needs: ggplot2, data.table, dplyr, tidyr, stringr (loaded by setup.R).
# ComplexHeatmap is loaded lazily inside the heatmap helper.
# ============================================================================ #

# ---- Differential expression --------------------------------------------- #

# Re-derive Up/Down/No from explicit thresholds (transparent; does not trust the
# table's own `diff`). Raw p_val, matching the analysis convention.
de_call <- function(de, p_thr = 0.05, lfc_thr = 0.5) {
  de %>% dplyr::mutate(diff = dplyr::case_when(
    p_val < p_thr & logFC >  lfc_thr ~ "Up",
    p_val < p_thr & logFC < -lfc_thr ~ "Down",
    TRUE ~ "No"))
}

# Volcano, one panel per contrast.
plot_volcano <- function(de, p_thr = 0.05, lfc_thr = 0.5, labels = CONTRAST_LABELS) {
  d <- de_call(de, p_thr, lfc_thr) %>%
    dplyr::mutate(contrast = factor(labels[contrast], levels = unname(labels)))
  ggplot(d, aes(logFC, -log10(p_val), colour = diff)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2, colour = "grey50") +
    geom_hline(yintercept = -log10(p_thr), linetype = 2, colour = "grey50") +
    scale_colour_manual(values = DIFF_COLS, name = NULL) +
    facet_wrap(~ contrast, nrow = 1) +
    labs(x = expression(log[2]~fold~change), y = expression(-log[10]~italic(P))) +
    theme_submission()
}

# MA plot (abundance vs log2FC), one panel per contrast. x = logCPM column.
plot_ma <- function(de, p_thr = 0.05, lfc_thr = 0.5, labels = CONTRAST_LABELS) {
  stopifnot("logCPM" %in% colnames(de))
  d <- de_call(de, p_thr, lfc_thr) %>%
    dplyr::mutate(contrast = factor(labels[contrast], levels = unname(labels)))
  ggplot(d, aes(logCPM, logFC, colour = diff)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_hline(yintercept = c(-lfc_thr, 0, lfc_thr), linetype = c(2,1,2), colour = "grey50") +
    scale_colour_manual(values = DIFF_COLS, name = NULL) +
    facet_wrap(~ contrast, nrow = 1) +
    labs(x = expression(mean~log[2]~CPM), y = expression(log[2]~fold~change)) +
    theme_submission()
}

# DE-gene counts (Up/Down per contrast).
plot_de_counts <- function(de, p_thr = 0.05, lfc_thr = 0.5, labels = CONTRAST_LABELS) {
  d <- de_call(de, p_thr, lfc_thr) %>%
    dplyr::filter(diff != "No") %>%
    dplyr::count(contrast, diff) %>%
    dplyr::mutate(n = ifelse(diff == "Down", -n, n),
                  contrast = factor(labels[contrast], levels = unname(labels)))
  ggplot(d, aes(contrast, n, fill = diff)) +
    geom_col() +
    geom_hline(yintercept = 0, colour = "black") +
    scale_fill_manual(values = DIFF_COLS[c("Up","Down")], name = NULL) +
    labs(x = NULL, y = "DE genes (Up above / Down below)",
         subtitle = sprintf("p < %g & |log2FC| > %g", p_thr, lfc_thr)) +
    theme_submission()
}

# Metacluster marker volcano (A/B/C/D vs rest). x = avg_log2FC, y = -log10 p_val_adj.
plot_marker_volcano <- function(markers, lfc_thr = 1, padj_thr = 0.01) {
  d <- markers %>% dplyr::mutate(
    sig = dplyr::case_when(
      p_val_adj < padj_thr & avg_log2FC >  lfc_thr ~ "Up",
      p_val_adj < padj_thr & avg_log2FC < -lfc_thr ~ "Down",
      TRUE ~ "No"),
    metacluster = factor(metacluster, levels = c("A","B","C","D")))
  ggplot(d, aes(avg_log2FC, -log10(p_val_adj + 1e-300), colour = sig)) +
    geom_point(size = 0.6, alpha = 0.6) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2, colour = "grey50") +
    geom_hline(yintercept = -log10(padj_thr), linetype = 2, colour = "grey50") +
    scale_colour_manual(values = DIFF_COLS, name = NULL) +
    facet_wrap(~ metacluster, nrow = 2) +
    labs(x = expression(log[2]~(metacluster/rest)), y = expression(-log[10]~adj.~italic(P))) +
    theme_submission()
}

# Pseudobulk heatmap of selected genes (rows z-scored). counts_dt has gene_name,
# gene_id + one column per library; sheet orders/annotates the columns.
plot_pseudobulk_heatmap <- function(counts_dt, genes, sheet,
                                     value = c("gene_name","gene_id")[1],
                                     title = NULL) {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  cnt  <- as.data.frame(counts_dt)
  keep <- cnt[[value]] %in% genes
  if (!any(keep)) stop("None of the requested genes are in the counts table.")
  samp_cols <- setdiff(colnames(cnt), c("gene_name","gene_id"))
  m <- as.matrix(cnt[keep, samp_cols, drop = FALSE])
  rownames(m) <- make.unique(cnt[[value]][keep])
  m <- log2(m + 1)
  m <- t(scale(t(m)))                                  # z-score per gene
  m <- m[, intersect(sheet$sample_id, colnames(m)), drop = FALSE]  # ss -> dd order
  ann <- sheet[match(colnames(m), sheet$sample_id), ]
  ha <- ComplexHeatmap::HeatmapAnnotation(
    experiment = ann$experiment, group = ann$group_id,
    col = list(group = CONDITION_COLS[intersect(names(CONDITION_COLS), unique(ann$group_id))]))
  ComplexHeatmap::Heatmap(
    m, name = "z(log2 pb)", top_annotation = ha,
    cluster_columns = FALSE, show_row_names = nrow(m) <= 60,
    column_title = title)
}

# ---- Gene-set permutation -------------------------------------------------- #

# Density of the permuted medians with observed median overlaid, per gene set x
# contrast. draws = dperm (columns gene_set_name, contrast, cluster, med, mean);
# dstat = summary (for data_median / perm_median). Optionally restrict to one
# cluster (default the single all-dopa grouping).
plot_perm_density <- function(draws, dstat, gene_set, cluster = NULL,
                              labels = CONTRAST_LABELS, xlim = 0.25) {
  dp <- draws %>% dplyr::filter(gene_set_name == gene_set)
  ds <- dstat %>% dplyr::filter(gene_set_name == gene_set)
  if (!is.null(cluster)) { dp <- dp %>% dplyr::filter(cluster == !!cluster)
                           ds <- ds %>% dplyr::filter(cluster == !!cluster) }
  lab <- function(x) factor(labels[x], levels = unname(labels))
  dp$contrast <- lab(dp$contrast); ds$contrast <- lab(ds$contrast)
  ggplot(dp, aes(med)) +
    geom_density(aes(fill = contrast), alpha = 0.2, colour = "white") +
    geom_vline(data = ds, aes(xintercept = perm_median), colour = "grey40") +
    geom_vline(data = ds, aes(xintercept = data_median, colour = contrast), linetype = 2) +
    scale_fill_manual(values = unname(CONTRAST_COLS)) +
    scale_colour_manual(values = unname(CONTRAST_COLS)) +
    facet_wrap(~ contrast) + coord_cartesian(xlim = c(-xlim, xlim)) +
    labs(x = "Median permuted log2FC", y = "Density", title = gene_set) +
    theme_submission() + theme(legend.position = "none")
}

# Median z-score heatmap: x = cluster, y = gene set, fill = median_z_score,
# stars for median_p_value < 0.05 (*) / < 0.001 (**). Facet by contrast.
plot_perm_zscore_heatmap <- function(dstat, cluster_levels = NULL,
                                     labels = CONTRAST_LABELS) {
  d <- dstat
  if (!is.null(cluster_levels))
    d <- d %>% dplyr::filter(cluster %in% cluster_levels) %>%
         dplyr::mutate(cluster = factor(cluster, levels = cluster_levels))
  d <- d %>% dplyr::mutate(
    contrast = factor(labels[contrast], levels = unname(labels)),
    star = dplyr::case_when(median_p_value < 0.001 ~ "**",
                            median_p_value < 0.05  ~ "*", TRUE ~ ""))
  lim <- max(abs(d$median_z_score), na.rm = TRUE)
  ggplot(d, aes(cluster, gene_set_name, fill = median_z_score)) +
    geom_tile(colour = "grey90") +
    geom_text(aes(label = star), size = 3, vjust = 0.75) +
    scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C",
                         midpoint = 0, limits = c(-lim, lim), name = "median Z") +
    facet_wrap(~ contrast, nrow = 1) +
    labs(x = NULL, y = NULL) +
    theme_submission() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          panel.grid = element_blank())
}
