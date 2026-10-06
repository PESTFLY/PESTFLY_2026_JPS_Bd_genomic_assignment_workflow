#!/usr/bin/env Rscript

# PESTFLY supplementary analysis: Reference exclusion sensitivity
#
# Purpose
# Compare the retained primary benchmark with exclusions of twelve published caution references,
# two lower completeness references, and their combined set of fourteen.
#
# Interpretation
# These are predefined sensitivity scenarios, not post hoc removals chosen to improve query
# assignments. Changes quantify dependence on the retained references; unchanged calls do not
# prove that every retained specimen is taxonomically or biologically unproblematic.
#
# Technical notes
# Scenario runs repeat Steps 02, 04 and 05 and require FASTA alignments. Results use the short
# 07_exclusion directory to limit Windows path length. Preexisting validated scenario outputs
# may be reused and are identified by run records.
#
# Run from the repository root:
#   Rscript steps/validation/07_reference_exclusion_sensitivity/run.R
# See the adjacent README.md for inputs, outputs and complete CLI defaults.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# =============================================================================
# General helpers
# =============================================================================

find_repo_root <- function(max_up = 10L) {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  candidate <- wd

  for (i in seq_len(max_up + 1L)) {
    if (dir.exists(file.path(candidate, "data")) &&
        dir.exists(file.path(candidate, "steps")) &&
        dir.exists(file.path(candidate, "results"))) {
      return(candidate)
    }

    parent <- normalizePath(
      file.path(candidate, ".."),
      winslash = "/",
      mustWork = FALSE
    )

    if (identical(parent, candidate)) break
    candidate <- parent
  }

  stop(
    "Cannot identify the repository root. Run this script from within the ",
    "pipeline repository."
  )
}

repo_root <- find_repo_root()

is_abs_path <- function(path) {
  grepl("^(?:[A-Za-z]:[/\\\\]|/)", path)
}

resolve_path <- function(path, root = repo_root) {
  if (is.null(path) || !nzchar(path)) return(path)

  if (is_abs_path(path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  normalizePath(file.path(root, path), winslash = "/", mustWork = FALSE)
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("[\u200B-\u200D\uFEFF]", "", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x
}

as_clean_logical <- function(x) {
  if (is.logical(x)) return(x)

  value <- tolower(clean_text(x))
  out <- rep(NA, length(value))

  out[value %in% c("true", "t", "1", "yes", "y", "reference", "ref")] <- TRUE
  out[value %in% c("false", "f", "0", "no", "n", "query", "intercept")] <- FALSE

  out
}

same_with_na <- function(x, y) {
  (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_cor <- function(x, y, method = "spearman") {
  x <- safe_numeric(x)
  y <- safe_numeric(y)
  keep <- is.finite(x) & is.finite(y)

  if (sum(keep) < 3L) return(NA_real_)
  if (length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) {
    return(NA_real_)
  }

  suppressWarnings(cor(x[keep], y[keep], method = method))
}

write_table_pair <- function(x, stem, out_dir) {
  x <- as.data.table(x)

  fwrite(
    x,
    file.path(out_dir, paste0(stem, ".tsv")),
    sep = "\t",
    na = "NA",
    quote = FALSE
  )

  saveRDS(x, file.path(out_dir, paste0(stem, ".rds")))
  invisible(x)
}

nonempty_files_exist <- function(paths) {
  if (length(paths) == 0L || !all(file.exists(paths))) return(FALSE)
  sizes <- file.info(paths)$size
  all(is.finite(sizes) & sizes > 0)
}

cli_value <- function(flag, value) {
  c(flag, shQuote(as.character(value)))
}

# =============================================================================
# Command line options
# =============================================================================

option_list <- list(
  make_option(
    "--qc_dir",
    type = "character",
    default = "results/01_qc",
    help = "Primary Step 01 QC directory [default %default]"
  ),
  make_option(
    "--align_dir",
    type = "character",
    default = "results/00_fasta",
    help = "Step 00 FASTA directory [default %default]"
  ),
  make_option(
    "--baseline_step2_dir",
    type = "character",
    default = "results/02_snp_panels",
    help = "Corrected baseline Step 02 directory [default %default]"
  ),
  make_option(
    "--baseline_step4_dir",
    type = "character",
    default = "results/04_origin_assignment",
    help = "Corrected baseline Step 04 directory [default %default]"
  ),
  make_option(
    "--baseline_step5_dir",
    type = "character",
    default = "results/05_loo_validation",
    help = "Corrected baseline Step 05 directory [default %default]"
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/validation/07_exclusion",
    help = "Supplementary validation output directory [default %default]"
  ),
  make_option(
    "--cores",
    type = "integer",
    default = 0L,
    help = "Step 02 workers; 0 uses all detected cores minus one [default %default]"
  ),
  make_option(
    "--chunk_size",
    type = "integer",
    default = 300L,
    help = "Step 02 FASTA files per processing chunk [default %default]"
  ),
  make_option(
    "--resume",
    type = "logical",
    default = TRUE,
    help = "Reuse complete scenario steps after an interrupted run [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$qc_dir <- resolve_path(opt$qc_dir)
opt$align_dir <- resolve_path(opt$align_dir)
opt$baseline_step2_dir <- resolve_path(opt$baseline_step2_dir)
opt$baseline_step4_dir <- resolve_path(opt$baseline_step4_dir)
opt$baseline_step5_dir <- resolve_path(opt$baseline_step5_dir)
opt$out_dir <- resolve_path(opt$out_dir)
opt$cores <- as.integer(opt$cores)
opt$chunk_size <- max(1L, as.integer(opt$chunk_size))

return_dir <- file.path(opt$out_dir, "return_bundle")
scenario_root <- file.path(opt$out_dir, "scenarios")

if (.Platform$OS.type == "windows") {
  longest_expected_path <- file.path(
    opt$out_dir,
    "scenarios",
    "exclude_published_caution",
    "02_snp_panels",
    "panels",
    "P1_macroregion_africa_vs_asia",
    "reference_group_counts.tsv"
  )

  if (nchar(longest_expected_path, type = "chars") > 245L) {
    stop(
      "The validation output path is too long for reliable Windows ",
      "file creation (", nchar(longest_expected_path, type = "chars"),
      " characters). Use a shorter repository location or run with ",
      "--out_dir results/v07."
    )
  }
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(return_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scenario_root, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Scripts, expected outputs, and baseline preflight
# =============================================================================

step2_script <- file.path(
  repo_root,
  "steps",
  "02_snp_panel_discovery_and_ranking",
  "run.R"
)

step4_script <- file.path(
  repo_root,
  "steps",
  "04_hierarchical_origin_assignment",
  "run.R"
)

step5_script <- file.path(
  repo_root,
  "steps",
  "05_full_panel_loo_validation",
  "run.R"
)

metadata_path <- file.path(opt$qc_dir, "metadata_clean.tsv")
ogs_pass_path <- file.path(opt$qc_dir, "ogs_pass_qc.txt")

panel_ids <- c(
  "P1_macroregion_africa_vs_asia",
  "P2_subregion_within_africa",
  "P3_subregion_within_asia"
)

loo_summary_files <- c(
  P1_macroregion_africa_vs_asia = "LOO_P1_macroregion_loo_summary.tsv",
  P2_subregion_within_africa = "LOO_P2_africa_subregion_loo_summary.tsv",
  P3_subregion_within_asia = "LOO_P3_asia_subregion_loo_summary.tsv"
)

loo_perclass_files <- c(
  P1_macroregion_africa_vs_asia = "LOO_P1_macroregion_loo_perclass.tsv",
  P2_subregion_within_africa = "LOO_P2_africa_subregion_loo_perclass.tsv",
  P3_subregion_within_asia = "LOO_P3_asia_subregion_loo_perclass.tsv"
)

required_inputs <- c(
  step2_script,
  step4_script,
  step5_script,
  metadata_path,
  ogs_pass_path,
  opt$align_dir,
  file.path(opt$baseline_step2_dir, "panels_index.tsv"),
  file.path(
    opt$baseline_step2_dir,
    "panels",
    panel_ids,
    "snp_map.tsv"
  ),
  file.path(
    opt$baseline_step2_dir,
    "panels",
    panel_ids,
    "snp_matrix.rds"
  ),
  file.path(opt$baseline_step4_dir, "final_assignment.tsv"),
  file.path(opt$baseline_step5_dir, unname(loo_summary_files)),
  file.path(opt$baseline_step5_dir, unname(loo_perclass_files))
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s):\n",
    paste0("  ", missing_inputs, collapse = "\n")
  )
}

step2_code <- paste(
  readLines(step2_script, warn = FALSE, encoding = "UTF-8"),
  collapse = "\n"
)

if (!grepl("hudson_fst", step2_code, fixed = TRUE)) {
  stop("The current Step 02 script does not contain hudson_fst().")
}

legacy_hudson_patterns <- c(
  "h1\\s*<-\\s*2\\s*\\*\\s*p1\\s*\\*\\s*\\(1\\s*-\\s*p1\\)",
  "h2\\s*<-\\s*2\\s*\\*\\s*p2\\s*\\*\\s*\\(1\\s*-\\s*p2\\)"
)

if (any(vapply(
  legacy_hudson_patterns,
  function(pattern) grepl(pattern, step2_code, perl = TRUE),
  logical(1)
))) {
  stop(
    "The current Step 02 script still contains the superseded Hudson FST ",
    "correction. Use the validated implementation before running this test."
  )
}

meta <- fread(metadata_path)

required_meta_cols <- c(
  "sample_id",
  "is_reference",
  "macroregion_3",
  "subregion",
  "country"
)

missing_meta_cols <- setdiff(required_meta_cols, names(meta))

if (length(missing_meta_cols) > 0L) {
  stop(
    "metadata_clean.tsv is missing required column(s): ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

meta[, sample_id := clean_text(sample_id)]
meta[, is_reference_clean := as_clean_logical(is_reference)]

if (anyNA(meta$is_reference_clean)) {
  stop("metadata_clean.tsv contains unrecognised is_reference values.")
}

if (anyNA(meta$sample_id) || anyDuplicated(meta$sample_id)) {
  stop("metadata_clean.tsv must contain unique, nonmissing sample_id values.")
}

n_references <- meta[is_reference_clean %in% TRUE, .N]
n_queries <- meta[is_reference_clean %in% FALSE, .N]

if (n_references != 330L || n_queries != 22L) {
  stop(
    "This script expects the frozen 330 reference and 22 query snapshot. ",
    "Observed references=", n_references, "; queries=", n_queries, "."
  )
}

# =============================================================================
# Evidence based exclusion sets
# =============================================================================

published_evidence <- data.table(
  sample_id = c(
    "CDKB02",
    "IDJI04",
    "LKAD01",
    "LKAD04",
    "MMYG01",
    "MMYG02",
    "MMYG04",
    "MMYG05",
    "YNJH01",
    "YNJH02",
    "YNJH03",
    "YNJH05"
  ),
  evidence_category = c(
    "COI_dorsalis_kandiensis_ambiguity",
    "nuclear_PCA_outlier; published_hybrid_candidate",
    "COI_dorsalis_kandiensis_ambiguity",
    "COI_dorsalis_kandiensis_ambiguity",
    "nuclear_PCA_outlier; published_hybrid_candidate",
    "nuclear_PCA_outlier; published_hybrid_candidate",
    "nuclear_PCA_outlier; published_hybrid_candidate",
    "nuclear_PCA_outlier; published_hybrid_candidate",
    "nuclear_PCA_outlier",
    "nuclear_PCA_outlier",
    "nuclear_PCA_outlier",
    "nuclear_PCA_outlier"
  ),
  evidence_source = c(
    "Vanbergen et al. 2025 Supplementary Table S15",
    "Vanbergen et al. 2025 Supplementary Table S9; Zhang et al. 2023",
    "Vanbergen et al. 2025 Supplementary Table S15",
    "Vanbergen et al. 2025 Supplementary Table S15",
    "Vanbergen et al. 2025 Supplementary Table S9; Zhang et al. 2023",
    "Vanbergen et al. 2025 Supplementary Table S9; Zhang et al. 2023",
    "Vanbergen et al. 2025 Supplementary Table S9; Zhang et al. 2023",
    "Vanbergen et al. 2025 Supplementary Table S9; Zhang et al. 2023",
    "Vanbergen et al. 2025 Supplementary Table S9",
    "Vanbergen et al. 2025 Supplementary Table S9",
    "Vanbergen et al. 2025 Supplementary Table S9",
    "Vanbergen et al. 2025 Supplementary Table S9"
  ),
  evidence_detail = c(
    "COI similarity included both B. dorsalis and B. kandiensis",
    "Outside the main nuclear PCA groups and reported with B. dorsalis/B. carambolae hybrid ancestry",
    "COI similarity included both B. dorsalis and B. kandiensis",
    "COI similarity included both B. dorsalis and B. kandiensis",
    "Outside the main nuclear PCA groups and reported with B. dorsalis/B. carambolae hybrid ancestry",
    "Outside the main nuclear PCA groups and reported with B. dorsalis/B. carambolae hybrid ancestry",
    "Outside the main nuclear PCA groups and reported with B. dorsalis/B. carambolae hybrid ancestry",
    "Outside the main nuclear PCA groups and reported with B. dorsalis/B. carambolae hybrid ancestry",
    "Outside the main B. dorsalis nuclear PCA groups",
    "Outside the main B. dorsalis nuclear PCA groups",
    "Outside the main B. dorsalis nuclear PCA groups",
    "Outside the main B. dorsalis nuclear PCA groups"
  )
)

low_coverage_ids <- c("ReuPie1", "ReuPie4")
all_expected_ids <- union(published_evidence$sample_id, low_coverage_ids)

missing_expected_ids <- setdiff(all_expected_ids, meta$sample_id)

if (length(missing_expected_ids) > 0L) {
  stop(
    "Expected sensitivity sample(s) are absent from metadata_clean.tsv: ",
    paste(missing_expected_ids, collapse = ", ")
  )
}

if (meta[sample_id %in% all_expected_ids, any(!is_reference_clean)]) {
  stop("Every sensitivity sample must remain marked as a reference.")
}

p1_matrix <- readRDS(file.path(
  opt$baseline_step2_dir,
  "panels",
  "P1_macroregion_africa_vs_asia",
  "snp_matrix.rds"
))

if (is.null(colnames(p1_matrix)) || is.null(rownames(p1_matrix))) {
  stop("The corrected baseline P1 snp_matrix.rds lacks row or column names.")
}

missing_matrix_samples <- setdiff(all_expected_ids, colnames(p1_matrix))

if (length(missing_matrix_samples) > 0L) {
  stop(
    "Expected sensitivity sample(s) are absent from the P1 matrix: ",
    paste(missing_matrix_samples, collapse = ", ")
  )
}

coverage_all_samples <- data.table(
  sample_id = colnames(p1_matrix),
  p1_markers = nrow(p1_matrix),
  p1_genotypes_called = colSums(!is.na(p1_matrix)),
  p1_call_fraction = colMeans(!is.na(p1_matrix))
)

coverage_evidence <- coverage_all_samples[
  sample_id %in% meta[is_reference_clean %in% TRUE, sample_id]
]

setorder(coverage_evidence, p1_call_fraction, sample_id)
coverage_evidence[, p1_coverage_rank_low_to_high := .I]

low_coverage_check <- coverage_evidence[sample_id %in% low_coverage_ids]

if (nrow(low_coverage_check) != 2L ||
    !setequal(low_coverage_check$p1_coverage_rank_low_to_high, 1:2)) {
  stop(
    "ReuPie1 and ReuPie4 are no longer the two lowest P1 call fraction ",
    "references. Review the frozen snapshot before proceeding."
  )
}

low_coverage_evidence <- low_coverage_check[, .(
  sample_id,
  evidence_category = "extreme_low_P1_genotype_completeness",
  evidence_source = "Current corrected P1 SNP matrix",
  evidence_detail = paste0(
    "P1 call fraction=",
    signif(p1_call_fraction, 6),
    "; rank ",
    p1_coverage_rank_low_to_high,
    " of ",
    nrow(coverage_evidence),
    " from lowest to highest"
  )
)]

reference_evidence <- rbindlist(
  list(published_evidence, low_coverage_evidence),
  use.names = TRUE,
  fill = TRUE
)

reference_evidence <- merge(
  reference_evidence,
  coverage_evidence[, .(
    sample_id,
    p1_markers,
    p1_genotypes_called,
    p1_call_fraction,
    p1_coverage_rank_low_to_high
  )],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

reference_evidence <- merge(
  reference_evidence,
  meta[, .(
    sample_id,
    original_macroregion_3 = macroregion_3,
    original_subregion = subregion,
    country,
    site = if ("site" %in% names(meta)) site else NA_character_,
    collection = if ("collection" %in% names(meta)) collection else NA_character_,
    sample_accession = if ("sample_accession" %in% names(meta)) sample_accession else NA_character_
  )],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(reference_evidence, evidence_category, sample_id)
write_table_pair(reference_evidence, "reference_exclusion_evidence", return_dir)

# =============================================================================
# Scenario definitions and isolated metadata snapshots
# =============================================================================

scenario_definitions <- data.table(
  scenario_id = c(
    "primary_all_references",
    "exclude_published_caution",
    "exclude_low_coverage",
    "exclude_combined"
  ),
  scenario_order = 1:4,
  scenario_label = c(
    "Corrected primary analysis with all references",
    "Exclude 12 published caution references",
    "Exclude ReuPie1 and ReuPie4",
    "Exclude all 14 caution or low coverage references"
  ),
  n_excluded = c(0L, 12L, 2L, 14L),
  run_pipeline = c(FALSE, TRUE, TRUE, TRUE)
)

scenario_exclusions <- list(
  primary_all_references = character(),
  exclude_published_caution = published_evidence$sample_id,
  exclude_low_coverage = low_coverage_ids,
  exclude_combined = union(published_evidence$sample_id, low_coverage_ids)
)

observed_exclusion_counts <- vapply(
  scenario_definitions$scenario_id,
  function(scenario_id) length(unique(scenario_exclusions[[scenario_id]])),
  integer(1)
)

if (!all(observed_exclusion_counts == scenario_definitions$n_excluded)) {
  stop("Internal scenario exclusion counts do not match their definitions.")
}

write_table_pair(scenario_definitions, "scenario_definitions", return_dir)

scenario_metadata <- list(primary_all_references = copy(meta))
metadata_change_list <- list()

for (sid in scenario_definitions[run_pipeline %in% TRUE, scenario_id]) {
  excluded_ids <- scenario_exclusions[[sid]]
  label <- scenario_definitions[scenario_id == sid, scenario_label]
  d <- copy(meta)

  before <- d[sample_id %in% excluded_ids, .(
    sample_id,
    country,
    site = if ("site" %in% names(d)) site else NA_character_,
    collection = if ("collection" %in% names(d)) collection else NA_character_,
    original_is_reference = is_reference,
    original_macroregion_3 = macroregion_3,
    original_subregion = subregion
  )]

  d[sample_id %in% excluded_ids, `:=`(
    macroregion_3 = "Others",
    subregion = "Other"
  )]

  after <- d[sample_id %in% excluded_ids, .(
    sample_id,
    scenario_is_reference = is_reference,
    scenario_macroregion_3 = macroregion_3,
    scenario_subregion = subregion
  )]

  changes <- merge(before, after, by = "sample_id", sort = FALSE)
  changes <- merge(
    changes,
    reference_evidence[, .(
      sample_id,
      evidence_category,
      evidence_source,
      evidence_detail,
      p1_markers,
      p1_genotypes_called,
      p1_call_fraction,
      p1_coverage_rank_low_to_high
    )],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  changes[, `:=`(
    scenario_id = sid,
    scenario_label = label
  )]

  setcolorder(
    changes,
    c(
      "scenario_id",
      "scenario_label",
      "sample_id",
      setdiff(names(changes), c("scenario_id", "scenario_label", "sample_id"))
    )
  )

  if (nrow(changes) != length(excluded_ids)) {
    stop("Scenario ", sid, " did not change every intended sample exactly once.")
  }

  if (d[sample_id %in% excluded_ids, any(!as_clean_logical(is_reference))]) {
    stop("Scenario ", sid, " accidentally changed reference status.")
  }

  if (d[sample_id %in% excluded_ids, any(
    macroregion_3 %in% c("Africa", "Asia") |
      subregion %in% c("C_Africa", "E_Africa", "S_Africa", "W_Africa", "E_Asia", "S_Asia", "SE_Asia")
  )]) {
    stop("Scenario ", sid, " failed to remove all excluded references from P1, P2, and P3.")
  }

  metadata_change_list[[sid]] <- changes
  scenario_metadata[[sid]] <- d
}

scenario_metadata_changes <- rbindlist(
  metadata_change_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(scenario_metadata_changes, scenario_id, sample_id)
write_table_pair(
  scenario_metadata_changes,
  "scenario_metadata_changes",
  return_dir
)

# =============================================================================
# Scenario paths and resumable execution
# =============================================================================

baseline_paths <- list(
  qc_dir = opt$qc_dir,
  step2_dir = opt$baseline_step2_dir,
  step4_dir = opt$baseline_step4_dir,
  step5_dir = opt$baseline_step5_dir
)

scenario_paths <- list(primary_all_references = baseline_paths)

for (sid in scenario_definitions[run_pipeline %in% TRUE, scenario_id]) {
  scenario_dir <- file.path(scenario_root, sid)
  qc_dir <- file.path(scenario_dir, "01_qc")
  step2_dir <- file.path(scenario_dir, "02_snp_panels")
  step4_dir <- file.path(scenario_dir, "04_origin_assignment")
  step5_dir <- file.path(scenario_dir, "05_loo_validation")

  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(step2_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(step4_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(step5_dir, recursive = TRUE, showWarnings = FALSE)

  scenario_meta_out <- copy(scenario_metadata[[sid]])
  scenario_meta_out[, is_reference_clean := NULL]

  fwrite(
    scenario_meta_out,
    file.path(qc_dir, "metadata_clean.tsv"),
    sep = "\t",
    na = "NA",
    quote = FALSE
  )

  saveRDS(scenario_meta_out, file.path(qc_dir, "metadata_clean.rds"))

  copied_ogs <- file.copy(
    ogs_pass_path,
    file.path(qc_dir, "ogs_pass_qc.txt"),
    overwrite = TRUE
  )

  if (!isTRUE(copied_ogs)) {
    stop("Could not copy ogs_pass_qc.txt for scenario ", sid, ".")
  }

  scenario_paths[[sid]] <- list(
    qc_dir = qc_dir,
    step2_dir = step2_dir,
    step4_dir = step4_dir,
    step5_dir = step5_dir
  )
}

rscript_bin <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

if (!file.exists(rscript_bin)) {
  stop("Cannot find Rscript at: ", rscript_bin)
}

run_log <- list()

run_external_step <- function(
    scenario_id,
    step_id,
    script_path,
    args,
    required_outputs,
    accept_nonzero_after_outputs = FALSE
) {
  started <- Sys.time()

  if (isTRUE(opt$resume) && nonempty_files_exist(required_outputs)) {
    message("\n", scenario_id, " ", step_id, ": complete outputs found; skipping.")

    run_log[[length(run_log) + 1L]] <<- data.table(
      scenario_id = scenario_id,
      step_id = step_id,
      status = "reused_complete_outputs",
      exit_status = 0L,
      started = as.character(started),
      finished = as.character(Sys.time()),
      elapsed_minutes = 0
    )

    return(invisible(0L))
  }

  message("\n============================================================")
  message("Scenario: ", scenario_id)
  message("Running:  ", step_id)
  message("============================================================")

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(repo_root)

  status <- system2(
    rscript_bin,
    args = c(shQuote(script_path), args),
    stdout = "",
    stderr = ""
  )

  status <- if (length(status) == 0L || is.null(status)) 0L else as.integer(status)
  outputs_ok <- nonempty_files_exist(required_outputs)
  finished <- Sys.time()

  if (status != 0L && !(isTRUE(accept_nonzero_after_outputs) && outputs_ok)) {
    stop(scenario_id, " ", step_id, " failed with status ", status, ".")
  }

  recorded_status <- if (status == 0L) {
    "completed"
  } else {
    warning(
      scenario_id,
      " ",
      step_id,
      " returned status ",
      status,
      " after all required tabular outputs were written. Continuing without ",
      "requiring its optional workbook."
    )
    "completed_required_outputs_nonzero_exit"
  }

  run_log[[length(run_log) + 1L]] <<- data.table(
    scenario_id = scenario_id,
    step_id = step_id,
    status = recorded_status,
    exit_status = status,
    started = as.character(started),
    finished = as.character(finished),
    elapsed_minutes = as.numeric(difftime(finished, started, units = "mins"))
  )

  invisible(status)
}

message("Repository root: ", repo_root)
message("Frozen references: ", n_references)
message("Frozen queries:    ", n_queries)
message("Published caution references: ", nrow(published_evidence))
message("Extreme low coverage references: ", length(low_coverage_ids))
message("Alternative scenarios to run: 3")
message("Step 02 cores option: ", opt$cores)

for (sid in scenario_definitions[run_pipeline %in% TRUE, scenario_id]) {
  paths <- scenario_paths[[sid]]

  step2_required <- c(
    file.path(paths$step2_dir, "panels_index.tsv"),
    file.path(paths$step2_dir, "panels", panel_ids, "snp_map.tsv"),
    file.path(paths$step2_dir, "panels", panel_ids, "snp_matrix.rds")
  )

  step2_args <- c(
    cli_value("--align_dir", opt$align_dir),
    cli_value("--qc_dir", paths$qc_dir),
    cli_value("--out_dir", paths$step2_dir),
    cli_value("--cores", opt$cores),
    cli_value("--chunk_size", opt$chunk_size),
    cli_value("--max_site_missing", 0.20),
    cli_value("--min_mac", 2L),
    cli_value("--min_ref_per_group", 3L),
    cli_value("--max_snps_per_og", 1L),
    cli_value("--panel_sizes", "20,50,100,200,500,1000,2000,5000")
  )

  run_external_step(
    scenario_id = sid,
    step_id = "Step_02",
    script_path = step2_script,
    args = step2_args,
    required_outputs = step2_required
  )

  step4_required <- c(
    file.path(paths$step4_dir, "final_assignment.tsv"),
    file.path(paths$step4_dir, "P1_macroregion_stability.tsv"),
    file.path(paths$step4_dir, "step4_K_grid_used.tsv")
  )

  step4_args <- c(
    cli_value("--qc_dir", paths$qc_dir),
    cli_value("--step3_dir", paths$step2_dir),
    cli_value("--out_dir", paths$step4_dir),
    cli_value("--out_xlsx", paste0("step4_", sid, ".xlsx"))
  )

  run_external_step(
    scenario_id = sid,
    step_id = "Step_04",
    script_path = step4_script,
    args = step4_args,
    required_outputs = step4_required,
    accept_nonzero_after_outputs = TRUE
  )

  step5_required <- c(
    file.path(paths$step5_dir, unname(loo_summary_files)),
    file.path(paths$step5_dir, unname(loo_perclass_files))
  )

  step5_args <- c(
    cli_value("--qc_dir", paths$qc_dir),
    cli_value("--step3_dir", paths$step2_dir),
    cli_value("--out_dir", paths$step5_dir)
  )

  run_external_step(
    scenario_id = sid,
    step_id = "Step_05",
    script_path = step5_script,
    args = step5_args,
    required_outputs = step5_required,
    accept_nonzero_after_outputs = TRUE
  )
}

run_status <- rbindlist(run_log, use.names = TRUE, fill = TRUE)

baseline_status <- data.table(
  scenario_id = "primary_all_references",
  step_id = c("Step_02", "Step_04", "Step_05"),
  status = "existing_corrected_primary_output",
  exit_status = 0L,
  started = NA_character_,
  finished = NA_character_,
  elapsed_minutes = NA_real_
)

run_status <- rbindlist(
  list(baseline_status, run_status),
  use.names = TRUE,
  fill = TRUE
)

write_table_pair(run_status, "scenario_run_status", return_dir)

# =============================================================================
# Reference counts under every scenario
# =============================================================================

make_panel_counts <- function(d, sid, label) {
  d <- copy(d)
  d[, is_reference_clean_internal := as_clean_logical(is_reference)]
  refs <- d[is_reference_clean_internal %in% TRUE]

  p1 <- refs[macroregion_3 %in% c("Africa", "Asia"), .N, by = .(
    analytical_group = macroregion_3
  )]
  p1[, panel_id := "P1_macroregion_africa_vs_asia"]

  p2 <- refs[macroregion_3 == "Africa" &
    subregion %in% c("C_Africa", "E_Africa", "S_Africa", "W_Africa"),
    .N,
    by = .(analytical_group = subregion)]
  p2[, panel_id := "P2_subregion_within_africa"]

  p3 <- refs[macroregion_3 == "Asia" &
    subregion %in% c("E_Asia", "S_Asia", "SE_Asia"),
    .N,
    by = .(analytical_group = subregion)]
  p3[, panel_id := "P3_subregion_within_asia"]

  out <- rbindlist(list(p1, p2, p3), use.names = TRUE, fill = TRUE)
  out[, `:=`(
    scenario_id = sid,
    scenario_label = label
  )]

  setcolorder(
    out,
    c("scenario_id", "scenario_label", "panel_id", "analytical_group", "N")
  )

  out
}

reference_count_list <- lapply(
  scenario_definitions$scenario_id,
  function(sid) {
    label <- scenario_definitions[scenario_id == sid, scenario_label]

    make_panel_counts(
      scenario_metadata[[sid]],
      sid,
      label
    )
  }
)

scenario_reference_counts <- rbindlist(
  reference_count_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  scenario_reference_counts,
  scenario_id,
  panel_id,
  analytical_group
)

write_table_pair(
  scenario_reference_counts,
  "scenario_reference_counts",
  return_dir
)

# =============================================================================
# Panel, LOO, and query assignment comparisons
# =============================================================================

panel_index_list <- list()
loo_summary_list <- list()
loo_perclass_list <- list()
assignment_list <- list()

for (sid in scenario_definitions$scenario_id) {
  paths <- scenario_paths[[sid]]
  label <- scenario_definitions[scenario_id == sid, scenario_label]

  panel_index <- fread(file.path(paths$step2_dir, "panels_index.tsv"))
  panel_index[, `:=`(
    scenario_id = sid,
    scenario_label = label
  )]
  setcolorder(
    panel_index,
    c("scenario_id", "scenario_label", setdiff(names(panel_index), c("scenario_id", "scenario_label")))
  )
  panel_index_list[[sid]] <- panel_index

  scenario_loo_summary <- rbindlist(
    lapply(panel_ids, function(pid) {
      path <- file.path(paths$step5_dir, loo_summary_files[[pid]])
      x <- fread(path)
      x[, requested_panel_id := pid]
      x
    }),
    use.names = TRUE,
    fill = TRUE
  )

  scenario_loo_summary[, `:=`(
    scenario_id = sid,
    scenario_label = label
  )]
  loo_summary_list[[sid]] <- scenario_loo_summary

  scenario_loo_perclass <- rbindlist(
    lapply(panel_ids, function(pid) {
      path <- file.path(paths$step5_dir, loo_perclass_files[[pid]])
      x <- fread(path)
      x[, requested_panel_id := pid]
      x
    }),
    use.names = TRUE,
    fill = TRUE
  )

  scenario_loo_perclass[, `:=`(
    scenario_id = sid,
    scenario_label = label
  )]
  loo_perclass_list[[sid]] <- scenario_loo_perclass

  assignment <- fread(file.path(paths$step4_dir, "final_assignment.tsv"))
  assignment[, `:=`(
    scenario_id = sid,
    scenario_label = label
  )]
  assignment_list[[sid]] <- assignment
}

panel_index_comparison <- rbindlist(
  panel_index_list,
  use.names = TRUE,
  fill = TRUE
)

loo_summary_comparison <- rbindlist(
  loo_summary_list,
  use.names = TRUE,
  fill = TRUE
)

loo_perclass_comparison <- rbindlist(
  loo_perclass_list,
  use.names = TRUE,
  fill = TRUE
)

query_assignment_long <- rbindlist(
  assignment_list,
  use.names = TRUE,
  fill = TRUE
)

write_table_pair(panel_index_comparison, "panel_index_comparison", return_dir)
write_table_pair(loo_summary_comparison, "loo_summary_comparison", return_dir)
write_table_pair(loo_perclass_comparison, "loo_perclass_comparison", return_dir)
write_table_pair(query_assignment_long, "query_assignment_long", return_dir)

assignment_fields <- c(
  "sample_id",
  "macroregion_call",
  "macroregion_confidence",
  "macroregion_reason",
  "macroregion_maxK_post",
  "macroregion_maxK_gap",
  "macroregion_agreement",
  "subregion_call",
  "subregion_confidence",
  "subregion_reason",
  "subregion_panel",
  "subregion_maxK_post",
  "subregion_maxK_gap",
  "subregion_agreement"
)

missing_assignment_fields <- setdiff(
  assignment_fields,
  names(assignment_list[["primary_all_references"]])
)

if (length(missing_assignment_fields) > 0L) {
  stop(
    "final_assignment.tsv is missing required comparison column(s): ",
    paste(missing_assignment_fields, collapse = ", ")
  )
}

assignment_change_list <- list()

baseline_assignment <- copy(assignment_list[["primary_all_references"]])
baseline_assignment <- baseline_assignment[, ..assignment_fields]

for (sid in setdiff(
  scenario_definitions$scenario_id,
  "primary_all_references"
)) {
  scenario_assignment <- copy(assignment_list[[sid]])
  scenario_assignment <- scenario_assignment[, ..assignment_fields]

  comparison <- merge(
    baseline_assignment,
    scenario_assignment,
    by = "sample_id",
    all = TRUE,
    suffixes = c("_baseline", "_scenario"),
    sort = TRUE
  )

  comparison[, `:=`(
    scenario_id = sid,
    scenario_label = scenario_definitions[scenario_id == sid, scenario_label],
    macroregion_call_changed = !same_with_na(
      macroregion_call_baseline,
      macroregion_call_scenario
    ),
    macroregion_confidence_changed = !same_with_na(
      macroregion_confidence_baseline,
      macroregion_confidence_scenario
    ),
    subregion_call_changed = !same_with_na(
      subregion_call_baseline,
      subregion_call_scenario
    ),
    subregion_confidence_changed = !same_with_na(
      subregion_confidence_baseline,
      subregion_confidence_scenario
    ),
    macroregion_posterior_abs_change = abs(
      safe_numeric(macroregion_maxK_post_scenario) -
        safe_numeric(macroregion_maxK_post_baseline)
    ),
    subregion_posterior_abs_change = abs(
      safe_numeric(subregion_maxK_post_scenario) -
        safe_numeric(subregion_maxK_post_baseline)
    )
  )]

  setcolorder(
    comparison,
    c(
      "scenario_id",
      "scenario_label",
      "sample_id",
      setdiff(names(comparison), c("scenario_id", "scenario_label", "sample_id"))
    )
  )

  assignment_change_list[[sid]] <- comparison
}

query_assignment_changes <- rbindlist(
  assignment_change_list,
  use.names = TRUE,
  fill = TRUE
)

query_assignment_change_summary <- query_assignment_changes[, .(
  n_queries = .N,
  n_macroregion_call_changed = sum(macroregion_call_changed, na.rm = TRUE),
  n_macroregion_confidence_changed = sum(
    macroregion_confidence_changed,
    na.rm = TRUE
  ),
  n_subregion_call_changed = sum(subregion_call_changed, na.rm = TRUE),
  n_subregion_confidence_changed = sum(
    subregion_confidence_changed,
    na.rm = TRUE
  ),
  max_macroregion_posterior_abs_change = if (all(is.na(
    macroregion_posterior_abs_change
  ))) NA_real_ else max(macroregion_posterior_abs_change, na.rm = TRUE),
  max_subregion_posterior_abs_change = if (all(is.na(
    subregion_posterior_abs_change
  ))) NA_real_ else max(subregion_posterior_abs_change, na.rm = TRUE)
), by = .(scenario_id, scenario_label)]

write_table_pair(query_assignment_changes, "query_assignment_changes", return_dir)
write_table_pair(
  query_assignment_change_summary,
  "query_assignment_change_summary",
  return_dir
)

# =============================================================================
# Marker ranking correlations and Top K overlap
# =============================================================================

top_k_values <- c(20L, 100L, 500L, 2000L, 5000L)
marker_overlap_list <- list()
marker_correlation_list <- list()

for (sid in setdiff(
  scenario_definitions$scenario_id,
  "primary_all_references"
)) {
  for (pid in panel_ids) {
    baseline_map <- fread(file.path(
      scenario_paths[["primary_all_references"]]$step2_dir,
      "panels",
      pid,
      "snp_map.tsv"
    ))

    scenario_map <- fread(file.path(
      scenario_paths[[sid]]$step2_dir,
      "panels",
      pid,
      "snp_map.tsv"
    ))

    required_map_cols <- c("og", "pos", "global_rank", "score")

    if (length(setdiff(required_map_cols, names(baseline_map))) > 0L ||
        length(setdiff(required_map_cols, names(scenario_map))) > 0L) {
      stop("A SNP map lacks og, pos, global_rank, or score for ", pid, ".")
    }

    baseline_map[, locus_key := paste(og, pos, sep = "::")]
    scenario_map[, locus_key := paste(og, pos, sep = "::")]

    if (anyDuplicated(baseline_map$locus_key) || anyDuplicated(scenario_map$locus_key)) {
      stop("A SNP map contains duplicated locus keys for ", pid, ".")
    }

    setorder(baseline_map, global_rank)
    setorder(scenario_map, global_rank)

    common <- merge(
      baseline_map[, .(
        locus_key,
        baseline_rank = global_rank,
        baseline_score = score
      )],
      scenario_map[, .(
        locus_key,
        scenario_rank = global_rank,
        scenario_score = score
      )],
      by = "locus_key",
      all = FALSE,
      sort = FALSE
    )

    marker_correlation_list[[length(marker_correlation_list) + 1L]] <- data.table(
      scenario_id = sid,
      scenario_label = scenario_definitions[scenario_id == sid, scenario_label],
      panel_id = pid,
      n_baseline_markers = nrow(baseline_map),
      n_scenario_markers = nrow(scenario_map),
      n_common_markers = nrow(common),
      common_fraction_of_baseline = nrow(common) / nrow(baseline_map),
      common_fraction_of_scenario = nrow(common) / nrow(scenario_map),
      spearman_rank_correlation_common = safe_cor(
        common$baseline_rank,
        common$scenario_rank,
        method = "spearman"
      ),
      spearman_score_correlation_common = safe_cor(
        common$baseline_score,
        common$scenario_score,
        method = "spearman"
      )
    )

    panel_k <- top_k_values[
      top_k_values <= max(nrow(baseline_map), nrow(scenario_map))
    ]

    for (k in panel_k) {
      baseline_n <- min(k, nrow(baseline_map))
      scenario_n <- min(k, nrow(scenario_map))
      baseline_keys <- baseline_map[seq_len(baseline_n), locus_key]
      scenario_keys <- scenario_map[seq_len(scenario_n), locus_key]
      overlap_n <- length(intersect(baseline_keys, scenario_keys))
      union_n <- length(union(baseline_keys, scenario_keys))

      marker_overlap_list[[length(marker_overlap_list) + 1L]] <- data.table(
        scenario_id = sid,
        scenario_label = scenario_definitions[scenario_id == sid, scenario_label],
        panel_id = pid,
        K_requested = k,
        K_baseline_used = baseline_n,
        K_scenario_used = scenario_n,
        overlap_n = overlap_n,
        overlap_fraction_of_baseline = overlap_n / baseline_n,
        overlap_fraction_of_scenario = overlap_n / scenario_n,
        jaccard = if (union_n > 0L) overlap_n / union_n else NA_real_
      )
    }
  }
}

marker_overlap_topK <- rbindlist(
  marker_overlap_list,
  use.names = TRUE,
  fill = TRUE
)

marker_rank_correlation <- rbindlist(
  marker_correlation_list,
  use.names = TRUE,
  fill = TRUE
)

write_table_pair(marker_overlap_topK, "marker_overlap_topK", return_dir)
write_table_pair(marker_rank_correlation, "marker_rank_correlation", return_dir)

# =============================================================================
# Compact source tables and human readable summary
# =============================================================================

for (sid in scenario_definitions$scenario_id) {
  paths <- scenario_paths[[sid]]

  file.copy(
    file.path(paths$step2_dir, "panels_index.tsv"),
    file.path(return_dir, paste0(sid, "_panels_index.tsv")),
    overwrite = TRUE
  )

  file.copy(
    file.path(paths$step4_dir, "final_assignment.tsv"),
    file.path(return_dir, paste0(sid, "_final_assignment.tsv")),
    overwrite = TRUE
  )

  for (pid in panel_ids) {
    short_panel <- sub("_.*$", "", pid)

    file.copy(
      file.path(paths$step5_dir, loo_summary_files[[pid]]),
      file.path(
        return_dir,
        paste0(sid, "_", short_panel, "_loo_summary.tsv")
      ),
      overwrite = TRUE
    )

    file.copy(
      file.path(paths$step5_dir, loo_perclass_files[[pid]]),
      file.path(
        return_dir,
        paste0(sid, "_", short_panel, "_loo_perclass.tsv")
      ),
      overwrite = TRUE
    )
  }
}

summary_lines <- c(
  "PESTFLY REFERENCE EXCLUSION SENSITIVITY",
  "",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Repository root: ", repo_root),
  "",
  "FROZEN SNAPSHOT",
  paste0("References: ", n_references),
  paste0("Queries: ", n_queries),
  paste0("Published caution references tested: ", nrow(published_evidence)),
  paste0("Extreme low coverage references tested: ", length(low_coverage_ids)),
  "",
  "SCENARIOS"
)

for (sid in scenario_definitions$scenario_id) {
  definition <- scenario_definitions[scenario_id == sid]
  summary_lines <- c(
    summary_lines,
    paste0(
      sid,
      ": ",
      definition$scenario_label,
      "; n excluded=",
      definition$n_excluded
    )
  )
}

summary_lines <- c(summary_lines, "", "LOO SUMMARY")

for (sid in scenario_definitions$scenario_id) {
  scenario_loo <- loo_summary_comparison[scenario_id == sid]

  for (pid in panel_ids) {
    row <- scenario_loo[requested_panel_id == pid]
    if (nrow(row) == 0L) next

    accuracy_value <- if ("accuracy_all" %in% names(row)) {
      signif(row$accuracy_all[1], 6)
    } else {
      NA_real_
    }

    uncertainty_value <- if ("uncertainty_rate" %in% names(row)) {
      signif(row$uncertainty_rate[1], 6)
    } else {
      NA_real_
    }

    n_value <- if ("n" %in% names(row)) row$n[1] else NA_integer_

    summary_lines <- c(
      summary_lines,
      paste0(
        sid,
        " | ",
        pid,
        " | n=", n_value,
        " | accuracy_all=", accuracy_value,
        " | uncertainty_rate=", uncertainty_value
      )
    )
  }
}

summary_lines <- c(summary_lines, "", "QUERY ASSIGNMENT CHANGES")

for (sid in setdiff(
  scenario_definitions$scenario_id,
  "primary_all_references"
)) {
  row <- query_assignment_change_summary[scenario_id == sid]

  summary_lines <- c(
    summary_lines,
    paste0(
      sid,
      " | macroregion calls changed=", row$n_macroregion_call_changed[1],
      "/", row$n_queries[1],
      " | subregion calls changed=", row$n_subregion_call_changed[1],
      "/", row$n_queries[1],
      " | macroregion confidence changed=", row$n_macroregion_confidence_changed[1],
      " | subregion confidence changed=", row$n_subregion_confidence_changed[1]
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  "INTERPRETATION",
  "This file reports numerical comparisons only. Biological interpretation should be made after reviewing the returned tables.",
  "",
  "SUMMARY EXPORTS",
  "Zip the complete return_bundle directory. Do not zip the much larger scenario directories."
)

writeLines(
  summary_lines,
  con = file.path(return_dir, "REFERENCE_EXCLUSION_SENSITIVITY_SUMMARY.txt"),
  useBytes = TRUE
)

run_info <- list(
  parameters = opt,
  repository_root = repo_root,
  scenario_definitions = scenario_definitions,
  scenario_exclusions = scenario_exclusions,
  scenario_paths = scenario_paths,
  run_status = run_status,
  n_references = n_references,
  n_queries = n_queries,
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(
  run_info,
  file.path(return_dir, "reference_exclusion_sensitivity_run_info.rds")
)

message("\nDone.")
message("Compact return folder: ", return_dir)
message("The summary export directory is available for inspection.")
