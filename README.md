# APombo cocaine VTA — submission figures

Single-nucleus RNA-seq of mouse VTA dopaminergic neurons after cocaine.
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

One environment for all four reports: **R ≥ 4.5** with **Seurat 5.x**. The full
declared list (with the pinned floor versions) is in `DESCRIPTION`.

Dependencies by source:

| source | packages |
|---|---|
| CRAN | `data.table`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, `plotly`, `htmlwidgets`, `rmarkdown`, `knitr`, `Seurat` (≥ 5.0.0), `SeuratObject` (≥ 5.0.0) |
| Bioconductor | `ComplexHeatmap` |

Which reports need what:

| report | needs |
|---|---|
| 02 DifferentialExpression, 03 Permutation | CRAN table stack + `ComplexHeatmap` only (no Seurat) |
| 01 QC, 04 MERFISH | the above **plus** `Seurat`/`SeuratObject` (01) and `plotly`/`htmlwidgets` (04) |

### Install dependencies

One command installs everything (CRAN + Bioconductor, via `BiocManager`):

```sh
Rscript install.R
```

`install.R` is idempotent (skips already-installed packages), gates on R ≥ 4.5,
and verifies the install at the end. Equivalent manual steps:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "data.table","dplyr","tidyr","stringr","ggplot2","rmarkdown","knitr",
  "Seurat","SeuratObject","plotly","htmlwidgets","ComplexHeatmap"))
```

### Reproducible lockfile (optional, recommended)

For a byte-for-byte reproducible environment, snapshot a real `renv.lock` **in
the environment you render in** (versions + hashes can only come from installed
packages, so this must be run on your machine, not committed blind):

```r
install.packages("renv")
renv::init(bare = TRUE)   # adopt this repo as an renv project
renv::snapshot()          # writes renv.lock from the installed versions
```

Commit the resulting `renv.lock` (+ `renv/activate.R`, `renv/settings.json`);
`.gitignore` already excludes the `renv/library/`. Anyone can then reproduce with
`renv::restore()`. Exact versions are also recorded by the `sessionInfo()` chunk
at the end of every rendered HTML.

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
