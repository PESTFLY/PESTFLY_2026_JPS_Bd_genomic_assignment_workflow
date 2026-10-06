#!/usr/bin/env Rscript
# Base R checks for report eligibility. Run from the repository root.
# These checks do not rerun any genomic analysis.
source("steps/07_authority_facing_report/reporting_criteria.R")
p <- read_reporting_parameters("results/04_origin_assignment/step4_params.tsv")
saved <- read.delim("results/04_origin_assignment/final_assignment.tsv",
                   stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
names(saved) <- tolower(names(saved))
benchmark <- assess_reporting_criteria(saved, p)
stopifnot(sum(benchmark$summary$resolution == "macroregion" & benchmark$summary$status == "PASS") == 19L,
          sum(benchmark$summary$resolution == "macroregion" & benchmark$summary$status == "FAIL") == 3L,
          sum(benchmark$summary$resolution == "subregion" & benchmark$summary$status == "PASS") == 12L,
          sum(benchmark$summary$resolution == "subregion" & benchmark$summary$status == "FAIL") == 7L,
          sum(benchmark$summary$resolution == "subregion" & benchmark$summary$status == "NOT_EVALUATED") == 3L,
          nrow(benchmark$detail) == 264L)

rules <- reporting_threshold_table(p)
base <- saved[which(saved$macroregion_confidence == "High" & saved$subregion_confidence == "High")[1], , drop = FALSE]
base$sample_id <- "threshold_boundary"
for (i in seq_len(nrow(rules))) base[[rules$observed_field[i]]] <- rules$threshold[i]
base$macroregion_confidence <- "Moderate"
base$subregion_confidence <- "Moderate"
stopifnot(all(assess_reporting_criteria(base, p)$summary$status == "PASS"))

expect_error <- function(expr) {
  failed <- tryCatch({ force(expr); FALSE }, error = function(e) TRUE)
  stopifnot(failed)
}

for (resolution in c("macroregion", "subregion")) {
  rr <- rules[rules$resolution == resolution, , drop = FALSE]
  for (i in seq_len(nrow(rr))) {
    x <- base
    delta <- if (rr$criterion[i] %in% c("usable_snps", "usable_K_count")) 1 else 1e-8
    x[[rr$observed_field[i]]] <- rr$threshold[i] - delta
    x[[paste0(resolution, "_confidence")]] <- "Uncertain"
    if (resolution == "macroregion") {
      subcols <- grep("^subregion_", names(x), value = TRUE)
      for (key in subcols) x[[key]] <- NA
      x$subregion_confidence <- "Uncertain"
    }
    report <- assess_reporting_criteria(x, p)
    failed <- report$summary[report$summary$resolution == resolution, , drop = FALSE]
    stopifnot(failed$status == "FAIL", failed$criteria_met == 5L,
              failed$criteria_unmet == 1L, failed$unmet_criteria == rr$criterion[i])
    if (resolution == "macroregion") {
      stopifnot(report$summary$status[report$summary$resolution == "subregion"] == "NOT_EVALUATED")
    }
  }
}

# Missing observed values fail individually. Invalid schemas and contradictory saved
# categories stop generation, rather than appearing as a valid PASS or silent downgrade.
x <- base
x$subregion_agreement <- NA_real_
x$subregion_confidence <- "Uncertain"
stopifnot(assess_reporting_criteria(x, p)$summary$unmet_criteria[2] == "global_K_agreement")
expect_error(assess_reporting_criteria(base[, names(base) != "macroregion_maxk_nsnps"], p))
expect_error(assess_reporting_criteria(rbind(base, base), p))
x <- base; x$macroregion_confidence <- "High"
expect_error(assess_reporting_criteria(x, p))
x <- base; x$subregion_tail_fraction <- 0.25
expect_error(assess_reporting_criteria(x, p))
x <- base; x$subregion_maxk_post <- 1.01
expect_error(assess_reporting_criteria(x, p))
x <- base; x$subregion_nk_usable_for_stability <- 6.5
expect_error(assess_reporting_criteria(x, p))
x <- base; x$macroregion_maxk_post <- p$high_posterior; x$macroregion_confidence <- "High"
stopifnot(assess_reporting_criteria(x, p)$summary$status[1] == "PASS")
x$macroregion_maxk_post <- p$high_posterior - 1e-8; x$macroregion_confidence <- "Moderate"
stopifnot(assess_reporting_criteria(x, p)$summary$status[1] == "PASS")

cat("Reporting checks passed: benchmark, inclusive thresholds, twelve single criterion failures,\n",
    "conditional branches, missing metrics and invalid inputs.\n", sep = "")
