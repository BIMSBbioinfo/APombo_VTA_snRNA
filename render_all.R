#!/usr/bin/env Rscript
# Render all submission figure reports to self-contained HTML in figures/.
# Usage:  Rscript render_all.R            # render all available reports
#         Rscript render_all.R 02 03      # render only matching reports
#
# Target environment: R 4.5 + Seurat 5.x (see DESCRIPTION / README).
# Missing input files cause the corresponding report to fail with a clear
# message from resolve_input(); other reports still render.

suppressPackageStartupMessages(library(rmarkdown))
root     <- normalizePath(dirname(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
reports  <- sort(list.files(file.path(root, "reports"), pattern = "\\.Rmd$", full.names = TRUE))
sel      <- commandArgs(trailingOnly = TRUE)
if (length(sel)) reports <- reports[grepl(paste(sel, collapse = "|"), basename(reports))]

# Report which inputs are present before rendering.
source(file.path(root, "R", "setup.R"))
message("\nInput status:"); print(input_status()); message("")

dir.create(file.path(root, "figures"), showWarnings = FALSE)
for (rmd in reports) {
  message(">>> rendering ", basename(rmd))
  ok <- tryCatch({
    rmarkdown::render(rmd, output_dir = file.path(root, "figures"),
                      envir = new.env(), quiet = TRUE); TRUE
  }, error = function(e) { message("  FAILED: ", conditionMessage(e)); FALSE })
  message(if (ok) "  OK" else "  skipped (see error above)")
}
message("\nDone. HTML + panels in: ", file.path(root, "figures"))
