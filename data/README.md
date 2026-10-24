# `data/` — required input files

Drop the exported submission files here. Files may keep their `YYMMDD.` date
prefix (the newest match wins) or be undated. Nothing in `data/` is committed
(see `.gitignore`), except this README.

| file (canonical name) | used by | description |
|---|---|---|
| `All_Filtered_Cells_Combined.rds` | 01 QC | All filtered nuclei (ss+dd), Seurat v5. |
| `Dopaminergic_Cells_Combined.rds` | 01 QC, (04) | Dopaminergic nuclei (ss+dd), Seurat v5. |
| `table-pseudobulk_counts_per_replicate_all.tsv` | 01 QC | Per-replicate pseudobulk counts, all nuclei (gene_name, gene_id, library columns). |
| `table-pseudobulk_counts_per_replicate_dopa.tsv` | 01 QC, 02 DE | Per-replicate pseudobulk counts, dopaminergic nuclei. |
| `table-pseudobulk_sample_sheet_all.tsv` | 01 QC | Library → condition map + n_nuclei (all). |
| `table-pseudobulk_sample_sheet_dopa.tsv` | 01 QC, 02 DE | Library → condition map + n_nuclei (dopa). |
| `table-differential_expression_ss_dd.tsv` | 02 DE | Combined ss+dd gene-level DE (edgeR TREAT); gene_id, gene_name, logFC, logCPM, p_val, contrast. |
| `table-metacluster_markers_wilcox.tsv` | 02 DE | A/B/C/D marker genes (Wilcoxon, metacluster vs rest). |
| `table-gene_set_permutation-all_dopa.tsv` | 03 Perm | Gene-set permutation summary, whole dopa. |
| `table-gene_set_permutation-per_subcluster.tsv` | 03 Perm | Permutation summary per DN subcluster (C00–C15). |
| `table-gene_set_permutation-per_metacluster.tsv` | 03 Perm | Permutation summary per metacluster (A–D). |
| `table-gene_set_permutation_draws-all_dopa.tsv` | 03 Perm | Raw per-permutation draws (null), whole dopa — needed for density plots. |
| `table-gene_set_permutation_draws-per_subcluster.tsv` | 03 Perm | Raw draws per subcluster. |
| `table-gene_set_permutation_draws-per_metacluster.tsv` | 03 Perm | Raw draws per metacluster. |
| `merfish_vta_cells.tsv` | 04 MERFISH | MERFISH VTA cells: brain_section_label, x, y, z, subclass, cluster_id, scrna_subcluster, metacluster. |
| `merfish_scrna_cluster_mapping.tsv` | 04 MERFISH | MERFISH cluster → best scRNA subcluster (correlation, metacluster). |
| `merfish_3d_objects.rds` | 04 MERFISH | Bundle: merfish_3d, meta_groups, subcluster_cols, cluster_mapping, make_3d_plot(). |

Column-level descriptions for the tables are in the exported data dictionaries
(`*column*description*.md`).

Check what's resolved from a running R session:
```r
source("R/setup.R"); input_status()
```
