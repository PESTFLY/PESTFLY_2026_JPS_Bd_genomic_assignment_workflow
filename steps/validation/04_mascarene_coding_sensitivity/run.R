#!/usr/bin/env Rscript

# PESTFLY supplementary analysis: Mascarene coding sensitivity
#
# Purpose
# Compare the baseline Asian lineage coding of 40 Reunion and five Mauritius references with
# geographic Africa and East Africa coding, and with their exclusion from candidate classes as
# Other.
#
# Interpretation
# This tests the effect of analytical class coding. Baseline Asia and Southeast Asia labels
# describe the lineage interpretation used by the benchmark, rather than the islands physical
# geography.
#
# Technical notes
# Scenario runs repeat Steps 02, 04 and 05 in isolated directories and compare resulting
# assignments and fixed panel validation. FASTA alignments are required for scenario panel
# construction. Preexisting validated scenario outputs may be reused; consult the run record.
#
# Run from the repository root:
#   Rscript steps/validation/04_mascarene_coding_sensitivity/run.R
# See the adjacent README.md for inputs, outputs and complete CLI defaults.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# =============================================================================
# Repository and path helpers
# =============================================================================

find_repo_root <- function(max_up = 10L) {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  cand <- wd

  for (i in seq_len(max_up + 1L)) {
    if (dir.exists(file.path(cand, "data")) &&
        dir.exists(file.path(cand, "steps")) &&
        dir.exists(file.path(cand, "results"))) {
      return(cand)
    }

    parent <- normalizePath(
      file.path(cand, ".."),
      winslash = "/",
      mustWork = FALSE
    )

    if (identical(parent, cand)) break
    cand <- parent
  }

  stop(
    "Cannot identify the repository root. Run this script from within the ",
    "pipeline repository."
  )
}

repo_root <- find_repo_root()

is_abs_path <- function(p) {
  grepl("^(?:[A-Za-z]:[/\\\\]|/)", p)
}

resolve_path <- function(p, root = repo_root) {
  if (is.null(p) || !nzchar(p)) return(p)

  if (is_abs_path(p)) {
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }

  normalizePath(file.path(root, p), winslash = "/", mustWork = FALSE)
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

  y <- tolower(clean_text(x))
  out <- rep(NA, length(y))

  out[y %in% c("true", "t", "1", "yes", "y", "reference", "ref")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n", "query", "intercept")] <- FALSE

  out
}

country_key <- function(x) {
  y <- clean_text(x)
  y <- iconv(y, from = "", to = "ASCII//TRANSLIT")
  y <- tolower(y)
  gsub("[^a-z]", "", y)
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
    default = "results/validation/04_mascarene_coding_sensitivity",
    help = "Supplementary validation output directory [default %default]"
  ),
  make_option(
    "--cores",
    type = "integer",
    default = 0L,
    help = "Step 02 workers; 0 lets Step 02 use all detected cores minus one [default %default]"
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
    help = "Reuse complete scenario steps from an interrupted earlier run [default %default]"
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

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(return_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(scenario_root, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Core script and baseline preflight
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
    c(
      "P1_macroregion_africa_vs_asia",
      "P2_subregion_within_africa",
      "P3_subregion_within_asia"
    ),
    "snp_map.tsv"
  ),
  file.path(opt$baseline_step4_dir, "final_assignment.tsv"),
  file.path(
    opt$baseline_step5_dir,
    "LOO_P1_macroregion_loo_summary.tsv"
  ),
  file.path(
    opt$baseline_step5_dir,
    "LOO_P2_africa_subregion_loo_summary.tsv"
  ),
  file.path(
    opt$baseline_step5_dir,
    "LOO_P3_asia_subregion_loo_summary.tsv"
  ),
  file.path(
    opt$baseline_step5_dir,
    "LOO_P1_macroregion_loo_perclass.tsv"
  ),
  file.path(
    opt$baseline_step5_dir,
    "LOO_P2_africa_subregion_loo_perclass.tsv"
  ),
  file.path(
    opt$baseline_step5_dir,
    "LOO_P3_asia_subregion_loo_perclass.tsv"
  )
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
    "correction. Replace it with the validated implementation before running ",
    "this sensitivity analysis."
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

meta[, is_reference_clean := as_clean_logical(is_reference)]

if (anyNA(meta$is_reference_clean)) {
  stop("metadata_clean.tsv contains unrecognised is_reference values.")
}

n_references <- meta[is_reference_clean %in% TRUE, .N]
n_queries <- meta[is_reference_clean %in% FALSE, .N]

if (n_references != 330L || n_queries != 22L) {
  stop(
    "This script expects the frozen 330 reference and 22 query snapshot. ",
    "Observed references=", n_references, "; queries=", n_queries, "."
  )
}

meta[, country_key_internal := country_key(country)]

mascarene_idx <- which(
  meta$is_reference_clean %in% TRUE &
    meta$country_key_internal %in% c("reunion", "mauritius")
)

mascarene_counts <- meta[mascarene_idx, .N, by = country][order(country)]

mascarene_meta_check <- meta[mascarene_idx]
n_reunion <- mascarene_meta_check[country_key(country) == "reunion", .N]
n_mauritius <- mascarene_meta_check[country_key(country) == "mauritius", .N]

if (length(mascarene_idx) != 45L ||
    n_reunion != 40L ||
    n_mauritius != 5L) {
  stop(
    "Expected 45 Mascarene references: 40 Reunion and 5 Mauritius. ",
    "Observed: ",
    paste0(mascarene_counts$country, "=", mascarene_counts$N, collapse = "; "),
    "."
  )
}

if (meta[mascarene_idx, any(macroregion_3 != "Asia" | subregion != "SE_Asia")]) {
  stop(
    "The baseline snapshot no longer codes every Reunion and Mauritius ",
    "reference as Asia and SE_Asia. Review the snapshot before continuing."
  )
}

meta[, c("is_reference_clean", "country_key_internal") := NULL]

message("Repository root: ", repo_root)
message("Frozen references: ", n_references)
message("Frozen queries:    ", n_queries)
message("Mascarene references: 45 (Reunion=40; Mauritius=5)")
message("Step 02 cores option: ", opt$cores)

# =============================================================================
# Scenario definitions and metadata snapshots
# =============================================================================

scenario_definitions <- data.table(
  scenario_id = c("baseline_asia", "geo_africa", "excluded_other"),
  scenario_order = 1:3,
  scenario_label = c(
    "Current genomic lineage coding",
    "Geographic Africa coding",
    "Mascarene references outside Africa Asia panels"
  ),
  macroregion_3 = c("Asia", "Africa", "Others"),
  subregion = c("SE_Asia", "E_Africa", "Other"),
  run_pipeline = c(FALSE, TRUE, TRUE)
)

write_table_pair(
  scenario_definitions,
  "scenario_definitions",
  return_dir
)

metadata_change_list <- list()
scenario_metadata <- list(baseline_asia = copy(meta))

for (sid in scenario_definitions[run_pipeline %in% TRUE, scenario_id]) {
  cfg <- scenario_definitions[scenario_id == sid]
  d <- copy(meta)

  before <- d[mascarene_idx, .(
    sample_id,
    country,
    site = if ("site" %in% names(d)) site else NA_character_,
    collection = if ("collection" %in% names(d)) collection else NA_character_,
    original_macroregion_3 = macroregion_3,
    original_subregion = subregion
  )]

  d[mascarene_idx, `:=`(
    macroregion_3 = cfg$macroregion_3,
    subregion = cfg$subregion
  )]

  after <- d[mascarene_idx, .(
    sample_id,
    scenario_macroregion_3 = macroregion_3,
    scenario_subregion = subregion
  )]

  changes <- merge(before, after, by = "sample_id", sort = FALSE)
  changes[, `:=`(
    scenario_id = sid,
    scenario_label = cfg$scenario_label
  )]

  setcolorder(
    changes,
    c(
      "scenario_id",
      "scenario_label",
      "sample_id",
      "country",
      "site",
      "collection",
      "original_macroregion_3",
      "original_subregion",
      "scenario_macroregion_3",
      "scenario_subregion"
    )
  )

  metadata_change_list[[sid]] <- changes
  scenario_metadata[[sid]] <- d
}

metadata_changes <- rbindlist(
  metadata_change_list,
  use.names = TRUE,
  fill = TRUE
)

write_table_pair(
  metadata_changes,
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

scenario_paths <- list(baseline_asia = baseline_paths)

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

  fwrite(
    scenario_metadata[[sid]],
    file.path(qc_dir, "metadata_clean.tsv"),
    sep = "\t",
    na = "NA",
    quote = FALSE
  )

  saveRDS(
    scenario_metadata[[sid]],
    file.path(qc_dir, "metadata_clean.rds")
  )

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
    stop(
      scenario_id,
      " ",
      step_id,
      " failed with status ",
      status,
      "."
    )
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

  step5_required <- file.path(
    paths$step5_dir,
    unname(loo_summary_files)
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
  scenario_id = "baseline_asia",
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
# Reference counts under each coding scenario
# =============================================================================

make_panel_counts <- function(d, scenario_id, scenario_label) {
  d <- copy(d)
  d[, is_reference_clean := as_clean_logical(is_reference)]
  refs <- d[is_reference_clean %in% TRUE]

  p1 <- refs[macroregion_3 %in% c("Africa", "Asia"), .N, by = .(
    analytical_group = macroregion_3
  )]
  p1[, panel_id := "P1_macroregion_africa_vs_asia"]

  p2 <- refs[macroregion_3 == "Africa" & !is.na(subregion), .N, by = .(
    analytical_group = subregion
  )]
  p2[, panel_id := "P2_subregion_within_africa"]

  p3 <- refs[macroregion_3 == "Asia" & !is.na(subregion), .N, by = .(
    analytical_group = subregion
  )]
  p3[, panel_id := "P3_subregion_within_asia"]

  out <- rbindlist(list(p1, p2, p3), use.names = TRUE, fill = TRUE)
  out[, `:=`(
    scenario_id = scenario_id,
    scenario_label = scenario_label
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
    label <- scenario_definitions[
      scenario_id == sid,
      scenario_label
    ]

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
  label <- scenario_definitions[
    scenario_id == sid,
    scenario_label
  ]

  panel_index <- fread(file.path(paths$step2_dir, "panels_index.tsv"))
  panel_index[, `:=`(
    scenario_id = sid,
    scenario_label = label
  )]
  setcolorder(panel_index, c("scenario_id", "scenario_label", setdiff(names(panel_index), c("scenario_id", "scenario_label"))))
  panel_index_list[[sid]] <- panel_index

  scenario_loo_summary <- rbindlist(
    lapply(panel_ids, function(panel_id) {
      path <- file.path(paths$step5_dir, loo_summary_files[[panel_id]])
      x <- fread(path)
      x[, requested_panel_id := panel_id]
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
    lapply(panel_ids, function(panel_id) {
      path <- file.path(paths$step5_dir, loo_perclass_files[[panel_id]])
      x <- fread(path)
      x[, requested_panel_id := panel_id]
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

write_table_pair(
  panel_index_comparison,
  "panel_index_comparison",
  return_dir
)

write_table_pair(
  loo_summary_comparison,
  "loo_summary_comparison",
  return_dir
)

write_table_pair(
  loo_perclass_comparison,
  "loo_perclass_comparison",
  return_dir
)

write_table_pair(
  query_assignment_long,
  "query_assignment_long",
  return_dir
)

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

assignment_change_list <- list()

baseline_assignment <- copy(assignment_list[["baseline_asia"]])
baseline_assignment <- baseline_assignment[, ..assignment_fields]

for (sid in setdiff(scenario_definitions$scenario_id, "baseline_asia")) {
  scenario_assignment <- copy(assignment_list[[sid]])
  scenario_assignment <- scenario_assignment[, ..assignment_fields]

  cmp <- merge(
    baseline_assignment,
    scenario_assignment,
    by = "sample_id",
    all = TRUE,
    suffixes = c("_baseline", "_scenario"),
    sort = TRUE
  )

  cmp[, `:=`(
    scenario_id = sid,
    scenario_label = scenario_definitions[
      scenario_id == sid,
      scenario_label
    ],
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
    cmp,
    c(
      "scenario_id",
      "scenario_label",
      "sample_id",
      setdiff(names(cmp), c("scenario_id", "scenario_label", "sample_id"))
    )
  )

  assignment_change_list[[sid]] <- cmp
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

write_table_pair(
  query_assignment_changes,
  "query_assignment_changes",
  return_dir
)

write_table_pair(
  query_assignment_change_summary,
  "query_assignment_change_summary",
  return_dir
)

# =============================================================================
# Marker ranking overlap with the corrected baseline
# =============================================================================

top_k_values <- c(20L, 50L, 100L, 200L, 500L, 1000L, 2000L, 5000L, 10000L)
marker_overlap_list <- list()
marker_correlation_list <- list()

for (sid in setdiff(scenario_definitions$scenario_id, "baseline_asia")) {
  for (panel_id in panel_ids) {
    baseline_map <- fread(file.path(
      scenario_paths[["baseline_asia"]]$step2_dir,
      "panels",
      panel_id,
      "snp_map.tsv"
    ))

    scenario_map <- fread(file.path(
      scenario_paths[[sid]]$step2_dir,
      "panels",
      panel_id,
      "snp_map.tsv"
    ))

    baseline_map[, locus_key := paste(og, pos, sep = "::")]
    scenario_map[, locus_key := paste(og, pos, sep = "::")]

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
      scenario_label = scenario_definitions[
        scenario_id == sid,
        scenario_label
      ],
      panel_id = panel_id,
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

    panel_k <- sort(unique(c(
      top_k_values,
      min(nrow(baseline_map), nrow(scenario_map))
    )))

    panel_k <- panel_k[panel_k > 0L]

    for (k in panel_k) {
      b_n <- min(k, nrow(baseline_map))
      s_n <- min(k, nrow(scenario_map))
      b_keys <- baseline_map[seq_len(b_n), locus_key]
      s_keys <- scenario_map[seq_len(s_n), locus_key]
      overlap_n <- length(intersect(b_keys, s_keys))
      union_n <- length(union(b_keys, s_keys))

      marker_overlap_list[[length(marker_overlap_list) + 1L]] <- data.table(
        scenario_id = sid,
        scenario_label = scenario_definitions[
          scenario_id == sid,
          scenario_label
        ],
        panel_id = panel_id,
        K_requested = k,
        K_baseline_used = b_n,
        K_scenario_used = s_n,
        overlap_n = overlap_n,
        overlap_fraction_of_baseline = overlap_n / b_n,
        overlap_fraction_of_scenario = overlap_n / s_n,
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

write_table_pair(
  marker_overlap_topK,
  "marker_overlap_topK",
  return_dir
)

write_table_pair(
  marker_rank_correlation,
  "marker_rank_correlation",
  return_dir
)

# =============================================================================
# Compact source tables for return and audit
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

  for (panel_id in panel_ids) {
    short_panel <- sub("_.*$", "", panel_id)

    file.copy(
      file.path(paths$step5_dir, loo_summary_files[[panel_id]]),
      file.path(
        return_dir,
        paste0(sid, "_", short_panel, "_loo_summary.tsv")
      ),
      overwrite = TRUE
    )

    file.copy(
      file.path(paths$step5_dir, loo_perclass_files[[panel_id]]),
      file.path(
        return_dir,
        paste0(sid, "_", short_panel, "_loo_perclass.tsv")
      ),
      overwrite = TRUE
    )
  }
}

# =============================================================================
# Human readable summary
# =============================================================================

summary_lines <- c(
  "PESTFLY MASCARENE CODING SENSITIVITY",
  "",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Repository root: ", repo_root),
  "",
  "FROZEN SNAPSHOT",
  paste0("References: ", n_references),
  paste0("Queries: ", n_queries),
  "Mascarene references changed by the sensitivity test: 45",
  "Reunion references: 40",
  "Mauritius references: 5",
  "",
  "SCENARIOS",
  "baseline_asia: Asia and SE_Asia, using existing corrected primary outputs",
  "geo_africa: Africa and E_Africa, with Steps 02, 04, and 05 rerun",
  "excluded_other: Others and Other, with Steps 02, 04, and 05 rerun",
  "",
  "LOO SUMMARY"
)

for (sid in scenario_definitions$scenario_id) {
  x <- loo_summary_comparison[scenario_id == sid]

  for (panel_id in panel_ids) {
    y <- x[requested_panel_id == panel_id]
    if (nrow(y) == 0L) next

    summary_lines <- c(
      summary_lines,
      paste0(
        sid,
        " | ",
        panel_id,
        " | n=", y$n[1],
        " | accuracy_all=", signif(y$accuracy_all[1], 6),
        " | uncertainty_rate=", signif(y$uncertainty_rate[1], 6)
      )
    )
  }
}

summary_lines <- c(summary_lines, "", "QUERY ASSIGNMENT CHANGES")

for (sid in setdiff(scenario_definitions$scenario_id, "baseline_asia")) {
  x <- query_assignment_change_summary[scenario_id == sid]

  summary_lines <- c(
    summary_lines,
    paste0(
      sid,
      " | macroregion calls changed=", x$n_macroregion_call_changed[1],
      "/", x$n_queries[1],
      " | subregion calls changed=", x$n_subregion_call_changed[1],
      "/", x$n_queries[1],
      " | macroregion confidence changed=", x$n_macroregion_confidence_changed[1],
      " | subregion confidence changed=", x$n_subregion_confidence_changed[1]
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
  con = file.path(return_dir, "MASCARENE_CODING_SENSITIVITY_SUMMARY.txt"),
  useBytes = TRUE
)

run_info <- list(
  parameters = opt,
  repository_root = repo_root,
  scenario_definitions = scenario_definitions,
  scenario_paths = scenario_paths,
  run_status = run_status,
  n_references = n_references,
  n_queries = n_queries,
  n_mascarene_references = length(mascarene_idx),
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(
  run_info,
  file.path(return_dir, "mascarene_coding_sensitivity_run_info.rds")
)

message("\nDone.")
message("Compact return folder: ", return_dir)
message("The summary export directory is available for inspection.")
