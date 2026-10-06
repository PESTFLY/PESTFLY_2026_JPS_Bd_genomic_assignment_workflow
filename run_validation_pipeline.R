#!/usr/bin/env Rscript

# Run the seven supplementary validation and sensitivity analyses.
# Execute from the repository root after restoring all required inputs.
# Analyses 04 and 07 require regenerated Step 00 FASTA alignments.
# Default output directories can replace archived supplementary outputs.

steps <- c(
  "steps/validation/01_out_of_scope_challenge/run.R",
  "steps/validation/02_grouped_geographic_cv/run.R",
  "steps/validation/03_reference_set_audit/run.R",
  "steps/validation/04_mascarene_coding_sensitivity/run.R",
  "steps/validation/05_reference_downsampling/run.R",
  "steps/validation/06_marker_resampling/run.R",
  "steps/validation/07_reference_exclusion_sensitivity/run.R"
)

missing <- steps[!file.exists(steps)]
if (length(missing) > 0L) {
  stop(
    "Run from the repository root. Missing scripts: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

rscript <- file.path(R.home("bin"), "Rscript")
if (.Platform$OS.type == "windows") rscript <- paste0(rscript, ".exe")

for (s in steps) {
  message("Running supplementary analysis: ", s)
  status <- system2(rscript, shQuote(s))
  if (!identical(status, 0L)) {
    stop("Supplementary analyses stopped because this script failed: ", s, call. = FALSE)
  }
}

message("Supplementary validation analyses completed.")
