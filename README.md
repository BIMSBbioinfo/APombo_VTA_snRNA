# APombo cocaine VTA — submission figures

Standalone, reproducible regeneration of the manuscript figures from the
**exported** data only. No dependency on the analysis project, its `Config.txt`,
or `Mylib` — clone, drop the exported files into `data/`, and knit.

Study: single-nucleus RNA-seq of mouse VTA dopaminergic neurons after cocaine.
Two cohorts — **ss** = single-injection (timepoints d01 = 24 h, d14) and
**dd** = re-injection (24 h after re-exposure).

## Layout

```
submission_repo/
├── R/
│   ├── setup.R          # libraries, input resolver, palettes, theme
│   └── plot_helpers.R   # self-contained figure functions
├── reports/
│   ├── 01_QC.Rmd                    # quality control (needs Seurat)
│   ├── 02_DifferentialExpression.Rmd# DE volcano / MA / heatmap / markers
│   ├── 03_Permutation.Rmd           # gene-set permutation density + z-heatmaps
│   └── 04_MERFISH.Rmd               # 3D projections + spatial slices (needs Seurat/plotly)
├── data/                # <- drop exported inputs here (gitignored); see data/README.md
├── figures/             # <- rendered HTML + panels (gitignored)
├── render_all.R
└── DESCRIPTION          # pinned environment
```

## Requirements

One environment for all four reports: **R 4.5 + Seurat 5.x**. Package list in
`DESCRIPTION`. 01/04 exercise Seurat/plotly/htmlwidgets; 02/03 need only
`data.table`, `dplyr`, `tidyr`, `ggplot2`, `ComplexHeatmap`.

## Use

1. Place the exported files in `data/` (see `data/README.md` for the list and
   the bimsbstatic URLs). Dated filenames are fine — the newest match wins.
2. Render:
   ```sh
   Rscript render_all.R          # all reports
   Rscript render_all.R 02 03    # only matching reports
   ```
   or knit a single `reports/*.Rmd` in RStudio.
3. Outputs (HTML + PDF/PNG panels + MERFISH 3D widgets) land in `figures/`.

Check which inputs are resolved:
```r
source("R/setup.R"); input_status()
```

## Report → inputs

| report | inputs |
|---|---|
| 01 QC | `*Combined.rds` (both) + pseudobulk sample sheets |
| 02 DifferentialExpression | `table-differential_expression_ss_dd.tsv` + pseudobulk counts/sheet (dopa) + `table-metacluster_markers_wilcox.tsv` |
| 03 Permutation | the three `gene_set_permutation-*` summaries + the three `gene_set_permutation_draws-*` draw tables |
| 04 MERFISH | `merfish_3d_objects.rds` + `merfish_vta_cells.tsv` + `merfish_scrna_cluster_mapping.tsv` |

## Notes

- `03_Permutation` density plots need the `*_draws-*` tables (the raw permutation
  draws). If only the summary tables are present, the density chunk errors with a
  clear "missing input" message; the z-score heatmaps still render.
- `git remote` is intentionally unset — add `git remote add origin <private repo>`
  when ready to push.
