# PESTFLY: final native R reporting verification.
# Run this file from the extracted repository root.
# In RStudio, set that directory as the working directory, then run:
# source("verify_release.R")
# This runs isolated function checks and the Step 07 reporting formatter.
# Genomic inference, marker discovery and statistical models are not rerun.

local({
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  required <- c(
    "tools/check_implementation.R",
    "steps/07_authority_facing_report/check_reporting.R",
    "steps/07_authority_facing_report/run.R",
    "results/04_origin_assignment/final_assignment.tsv",
    "results/04_origin_assignment/step4_params.tsv",
    "results/07_authority_report/FINAL_origin_tracing_authority_criteria.tsv",
    "results/07_authority_report/reporting_thresholds.tsv",
    "results/07_authority_report/FINAL_origin_tracing_authority_simplified.tsv",
    "results/07_authority_report/FINAL_origin_tracing_authority_extended.tsv"
  )
  missing_files <- required[!file.exists(required)]
  if (length(missing_files)) {
    stop("Set the working directory to the extracted repository root. Missing: ",
         paste(missing_files, collapse = ", "), call. = FALSE)
  }
  packages <- c("optparse", "data.table", "openxlsx")
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) {
    stop("Install the missing reporting packages, then run this checker again: ",
         paste(missing_packages, collapse = ", "), call. = FALSE)
  }
  executable <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  candidates <- c(file.path(R.home("bin"), executable),
                  file.path(R.home(), "bin", executable),
                  file.path(R.home(), "bin", "x64", executable))
  candidates <- unique(candidates[file.exists(candidates)])
  if (!length(candidates)) stop("Rscript was not found in the running R installation.", call. = FALSE)
  rscript <- normalizePath(candidates[1], winslash = "/", mustWork = TRUE)

  stamp <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  out <- file.path(root, "release_check", stamp)
  if (dir.exists(out)) stop("Verification output directory already exists: ", out)
  dir.create(out, recursive = TRUE)
  report_dir <- file.path(out, "authority_report")
  log_file <- file.path(out, "verification_log.txt")
  lines <- c("PESTFLY native R reporting verification", paste("Started:", Sys.time()),
             paste("R:", R.version.string), paste("Repository:", root),
             "No genomic or statistical model rerun is requested.")
  log <- function(text) {
    lines <<- c(lines, text)
    writeLines(lines, log_file, useBytes = TRUE)
    cat(paste(text, collapse = "\n"), "\n", sep = "")
  }
  writeLines(capture.output(sessionInfo()), file.path(out, "session_info.txt"))
  versions <- data.frame(package = packages,
                         version = vapply(packages, function(x) as.character(packageVersion(x)), character(1)))
  write.table(versions, file.path(out, "reporting_package_versions.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  archived <- sort(list.files(file.path(root, "results"), recursive = TRUE, full.names = TRUE))
  before <- tools::md5sum(archived)
  saveRDS(before, file.path(out, "archived_result_md5_before.rds"))
  log("Outputs are written to a new release_check directory.")

  run <- function(label, script, args = character()) {
    log(paste("RUN:", label))
    output <- system2(rscript, args = c("--vanilla", shQuote(script), shQuote(args)),
                      stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    log(output)
    if (status != 0L) stop(label, " failed. See ", log_file, call. = FALSE)
    log(paste("PASS:", label))
  }

  read_tsv <- function(path) {
    if (!file.exists(path)) stop("Missing output: ", path)
    x <- read.delim(path, colClasses = "character", stringsAsFactors = FALSE,
                    check.names = FALSE, na.strings = c("", "NA"), fileEncoding = "UTF-8")
    names(x) <- tolower(names(x))
    if (anyDuplicated(names(x))) stop("Duplicate table columns: ", path)
    x
  }
  align <- function(x, keys, label) {
    if (!all(keys %in% names(x))) stop("Missing comparison key in ", label)
    if (anyNA(x[keys])) stop("Missing comparison key value in ", label)
    key <- do.call(paste, c(unname(x[keys]), sep = "\r"))
    if (anyDuplicated(key)) stop("Duplicate comparison key in ", label)
    list(data = x[order(key), , drop = FALSE], key = sort(key))
  }
  same_values <- function(a, b) {
    a <- as.character(a); b <- as.character(b)
    same <- is.na(a) & is.na(b)
    present <- !is.na(a) & !is.na(b)
    same[present] <- a[present] == b[present]
    # Ignore numeric text formatting and logical letter case, while preserving
    # missing values and exact scientific category labels.
    an <- suppressWarnings(as.numeric(a)); bn <- suppressWarnings(as.numeric(b))
    numeric <- present & is.finite(an) & is.finite(bn)
    same[numeric] <- abs(an[numeric] - bn[numeric]) <= 1e-12 * pmax(1, abs(an[numeric]), abs(bn[numeric]))
    logical <- present & toupper(a) %in% c("TRUE", "FALSE") & toupper(b) %in% c("TRUE", "FALSE")
    same[logical] <- toupper(a[logical]) == toupper(b[logical])
    same
  }
  comparisons <- list()
  # These two descriptive fields used "Step 4" in the independently assembled
  # checkpoint and "Step 04" in native R. Normalise this exact label spelling
  # only. Changed scientific values or other descriptive text still fail.
  label_fields <- c("rf_corroboration_detail", "commodity_origin_mismatch_detail")
  same_field_values <- function(a, b, field) {
    if (field %in% label_fields) {
      a <- gsub("Step 4 ", "Step 04 ", as.character(a), fixed = TRUE)
      b <- gsub("Step 4 ", "Step 04 ", as.character(b), fixed = TRUE)
    }
    same_values(a, b)
  }
  compare_tables <- function(label, expected, observed, keys, ignore = character()) {
    a <- align(expected, keys, paste(label, "checkpoint"))
    b <- align(observed, keys, paste(label, "native R"))
    if (!identical(a$key, b$key)) stop("Row keys differ: ", label)
    fields <- setdiff(names(a$data), ignore)
    if (!all(fields %in% names(b$data))) stop("Native output lacks columns in ", label, ": ",
                                              paste(setdiff(fields, names(b$data)), collapse = ", "))
    failures <- lapply(fields, function(field) {
      bad <- which(!same_field_values(a$data[[field]], b$data[[field]], field))
      if (!length(bad)) return(NULL)
      data.frame(table = label, key = a$key[bad], field = field,
                 checkpoint = a$data[[field]][bad], native_R = b$data[[field]][bad],
                 stringsAsFactors = FALSE)
    })
    failures <- Filter(Negate(is.null), failures)
    comparisons[[length(comparisons) + 1L]] <<- data.frame(
      table = label, rows = nrow(a$data), fields = length(fields),
      result = if (length(failures)) "FAIL" else "PASS", stringsAsFactors = FALSE)
    if (length(failures)) {
      write.table(do.call(rbind, failures), file.path(out, "comparison_failures.tsv"),
                  sep = "\t", row.names = FALSE, quote = TRUE)
      stop("Saved values differ from native R output in ", label, ". See comparison_failures.tsv.")
    }
    log(paste("PASS:", label, "comparison;", nrow(a$data), "rows,", length(fields), "fields"))
  }

  success <- tryCatch({
    run("Hudson FST and MI implementation checks", "tools/check_implementation.R")
    run("Six criterion reporting checks", "steps/07_authority_facing_report/check_reporting.R")
    run("Step 07 formatter", "steps/07_authority_facing_report/run.R",
        c("--out_dir", report_dir))
    reference <- file.path(root, "results/07_authority_report")
    compare_file <- function(file, keys, ignore = character()) {
      compare_tables(file, read_tsv(file.path(reference, file)),
                     read_tsv(file.path(report_dir, file)), keys, ignore)
    }
    compare_file("FINAL_origin_tracing_authority_criteria.tsv", c("sample_id", "resolution", "criterion"))
    compare_file("reporting_thresholds.tsv", c("resolution", "criterion"))
    prose <- c("reporting_basis", "subregion_withheld_explanation", "action_note", "report_statement")
    compare_file("FINAL_origin_tracing_authority_extended.tsv", "sample_id", prose)
    compare_file("FINAL_origin_tracing_authority_simplified.tsv", "sample_id", prose)

    workbook <- file.path(report_dir, "FINAL_origin_tracing_authority_REPORT.xlsx")
    expected_sheets <- c("Simplified", "Extended", "Criteria", "Reporting_Thresholds", "Legend", "Column_Key",
                         "Intercept_Origin_Map", "Country_Region_Map", "Vanbergen_Reference", "Step4b_Validation",
                         "Step4b_PerClass", "Step5_RF_Validation", "Step5_RF_PerClass")
    if (!setequal(openxlsx::getSheetNames(workbook), expected_sheets)) stop("Workbook sheet names differ.")
    for (item in list(c("Simplified", "FINAL_origin_tracing_authority_simplified.tsv", "sample_id"),
                      c("Extended", "FINAL_origin_tracing_authority_extended.tsv", "sample_id"),
                      c("Criteria", "FINAL_origin_tracing_authority_criteria.tsv", "sample_id", "resolution", "criterion"),
                      c("Reporting_Thresholds", "reporting_thresholds.tsv", "resolution", "criterion"))) {
      sheet <- openxlsx::read.xlsx(workbook, sheet = item[1], check.names = FALSE)
      names(sheet) <- tolower(names(sheet))
      compare_tables(paste("Workbook", item[1]), read_tsv(file.path(report_dir, item[2])), sheet, item[-c(1, 2)])
    }
    for (kind in c("simplified", "extended")) {
      stem <- paste0("FINAL_origin_tracing_authority_", kind)
      obj <- as.data.frame(readRDS(file.path(report_dir, paste0(stem, ".rds"))))
      names(obj) <- tolower(names(obj))
      compare_tables(paste("RDS", kind), read_tsv(file.path(report_dir, paste0(stem, ".tsv"))), obj, "sample_id")
    }
    run_info <- readRDS(file.path(report_dir, "step7_run_info.rds"))
    if (!isTRUE(run_info$decision_consistency_with_step4) || run_info$n_samples != 22L || is.null(run_info$session_info)) {
      stop("Native R run record is incomplete or inconsistent.")
    }
    TRUE
  }, error = function(error) {
    log(paste("FAIL:", conditionMessage(error)))
    FALSE
  })
  after <- tools::md5sum(archived)
  unchanged <- identical(before, after) && identical(archived, sort(list.files(
    file.path(root, "results"), recursive = TRUE, full.names = TRUE)))
  log(if (unchanged) "PASS: archived results remain unchanged." else "FAIL: archived result files changed.")
  success <- success && unchanged
  if (length(comparisons)) write.table(do.call(rbind, comparisons), file.path(out, "comparison_summary.tsv"),
                                       sep = "\t", row.names = FALSE, quote = FALSE)
  summary <- list(status = if (success) "PASS" else "FAIL", completed = Sys.time(),
                  R_version = R.version.string, reporting_package_versions = versions,
                  archived_results_unchanged = unchanged,
                  output_directory = out, genomic_models_recomputed = FALSE)
  saveRDS(summary, file.path(out, "release_check_summary.rds"))
  writeLines(c(paste("Overall result:", summary$status), paste("R:", R.version.string),
                paste("Archived results unchanged:", unchanged), "Genomic models recomputed: FALSE",
                paste("Completed:", Sys.time())), file.path(out, "release_check_summary.txt"))
  log(paste("Overall result:", summary$status))
  log(paste("Results folder:", normalizePath(out, winslash = "/", mustWork = TRUE)))
  if (!success) stop("Verification did not pass. Return this results folder for review.", call. = FALSE)
  cat("\nVerification passed. Return this results folder for review.\n")
})
