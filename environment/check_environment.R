#!/usr/bin/env Rscript

# Read-only check of the main dependencies recorded for the benchmark.
# Run from the repository root; --strict also rejects version differences.

args <- commandArgs(trailingOnly = TRUE)
strict <- "--strict" %in% args
versions_file <- "environment/benchmark_versions.tsv"
if (!file.exists(versions_file)) stop("Run this check from the repository root.")
expected <- read.delim(versions_file, stringsAsFactors = FALSE)

observed <- vapply(expected$component, function(package) {
  if (package == "R") return(as.character(getRversion()))
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(package))
}, character(1))

status <- ifelse(is.na(observed), "MISSING",
                 ifelse(observed == expected$benchmark_version, "MATCH", "DIFFERENT"))
result <- data.frame(component = expected$component,
                     benchmark_version = expected$benchmark_version,
                     installed_version = unname(observed), status = unname(status))
print(result, row.names = FALSE)

required <- c("data.table", "optparse", "readxl", "Biostrings", "openxlsx")
missing_required <- result$component[result$status == "MISSING" & result$component %in% required]
rf_available <- any(result$component %in% c("ranger", "randomForest") & result$status != "MISSING")
if (length(missing_required) > 0L || !rf_available) {
  stop("Required dependencies are missing. See docs/INSTALL.md.")
}
if (strict && any(result$status != "MATCH")) {
  stop("Installed main versions differ from the recorded benchmark environment.")
}
if (any(result$status != "MATCH")) {
  message("Differences are reported above. Exact benchmark reproduction requires the recorded environment and RF backend.")
} else {
  message("All recorded main component versions match.")
}
