#!/usr/bin/env Rscript

# PESTFLY: Geographic out of scope challenge
#
# Treat five Hawaiian and five Papua New Guinean Other references as unknown queries
# against Africa and Asia classes.
# Support compares represented classes; it cannot detect an absent source.
#
# Run from the repository root:
#   Rscript steps/validation/01_out_of_scope_challenge/run.R
# Inputs, outputs and options: adjacent README.md. Full methods: associated manuscript.

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
        dir.exists(file.path(cand, "steps"))) {
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

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  median(x)
}

# =============================================================================
# Command line options
# =============================================================================

option_list <- list(
  make_option(
    "--metadata",
    type = "character",
    default = "results/01_qc/metadata_clean.tsv",
    help = "Original cleaned metadata [default %default]"
  ),
  make_option(
    "--step3_dir",
    type = "character",
    default = "results/02_snp_panels",
    help = "Corrected Step 02 panel directory [default %default]"
  ),
  make_option(
    "--step4_script",
    type = "character",
    default = "steps/04_hierarchical_origin_assignment/run.R",
    help = "Existing Step 04 assignment script [default %default]"
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/validation/01_out_of_scope_challenge",
    help = "Separate challenge output directory [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$metadata <- resolve_path(opt$metadata)
opt$step3_dir <- resolve_path(opt$step3_dir)
opt$step4_script <- resolve_path(opt$step4_script)
opt$out_dir <- resolve_path(opt$out_dir)

challenge_input_dir <- file.path(opt$out_dir, "challenge_input")
step4_out_dir <- file.path(opt$out_dir, "step04_closed_set")

dir.create(challenge_input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(step4_out_dir, recursive = TRUE, showWarnings = FALSE)

message("Repository root:  ", repo_root)
message("Original metadata: ", opt$metadata)
message("Step 02 panels:    ", opt$step3_dir)
message("Step 04 script:    ", opt$step4_script)
message("Challenge output:  ", opt$out_dir)

if (!file.exists(opt$metadata)) {
  stop("Missing metadata file: ", opt$metadata)
}

if (!dir.exists(opt$step3_dir)) {
  stop("Missing Step 02 panel directory: ", opt$step3_dir)
}

if (!file.exists(opt$step4_script)) {
  stop("Missing Step 04 script: ", opt$step4_script)
}

# =============================================================================
# Identify the challenge samples from the unmodified metadata
# =============================================================================

meta <- fread(opt$metadata)
setnames(meta, tolower(names(meta)))

required_meta <- c(
  "sample_id",
  "is_reference",
  "macroregion_3",
  "subregion",
  "country",
  "site",
  "population"
)

missing_meta <- setdiff(required_meta, names(meta))

if (length(missing_meta) > 0L) {
  stop(
    "Metadata is missing required column(s): ",
    paste(missing_meta, collapse = ", ")
  )
}

for (nm in names(meta)) {
  if (!inherits(meta[[nm]], c("numeric", "integer", "logical", "Date", "POSIXct"))) {
    meta[[nm]] <- clean_text(meta[[nm]])
  }
}

meta[, sample_id := clean_text(sample_id)]
meta[, is_reference := as_clean_logical(is_reference)]
meta[, macroregion_3 := clean_text(macroregion_3)]
meta[, subregion := clean_text(subregion)]
meta[, country := clean_text(country)]
meta[, site := clean_text(site)]
meta[, population := clean_text(population)]

challenge_truth <- meta[
  is_reference %in% TRUE &
    macroregion_3 == "Others" &
    (
      site == "Hawaii" |
        country == "Papua_New_Guinea"
    ),
  .(
    sample_id,
    withheld_collection_class = macroregion_3,
    collection_subregion_label = subregion,
    collection_country = country,
    collection_site = site,
    collection_population = population,
    original_is_reference = is_reference
  )
]

setorder(challenge_truth, collection_site, sample_id)

if (nrow(challenge_truth) != 10L) {
  stop(
    "Expected 10 challenge samples but found ", nrow(challenge_truth),
    ". Review the Hawaii and Papua New Guinea metadata."
  )
}

if (uniqueN(challenge_truth$sample_id) != nrow(challenge_truth)) {
  stop("Challenge sample identifiers are not unique.")
}

expected_locations <- c("Hawaii", "Port_Moresby")

if (!setequal(unique(challenge_truth$collection_site), expected_locations)) {
  stop(
    "Unexpected challenge sites: ",
    paste(sort(unique(challenge_truth$collection_site)), collapse = ", ")
  )
}

message("Challenge samples identified: ", nrow(challenge_truth))
message(
  "  Hawaii: ",
  challenge_truth[collection_site == "Hawaii", .N]
)
message(
  "  Papua New Guinea: ",
  challenge_truth[collection_country == "Papua_New_Guinea", .N]
)

fwrite(
  challenge_truth,
  file.path(opt$out_dir, "challenge_truth.tsv"),
  sep = "\t"
)

saveRDS(
  challenge_truth,
  file.path(opt$out_dir, "challenge_truth.rds")
)

# =============================================================================
# Create an isolated metadata snapshot for the challenge
# =============================================================================
#
# Only original reference samples are retained.  Africa and Asia remain
# references.  The ten known Others samples become pseudo unknown queries.
# The 22 interception queries are deliberately omitted from this snapshot.
#

challenge_ids <- challenge_truth$sample_id

meta_challenge <- copy(meta[is_reference %in% TRUE])
meta_challenge[, is_reference := !(sample_id %in% challenge_ids)]

if (meta_challenge[is_reference %in% FALSE, .N] != 10L) {
  stop("The challenge metadata does not contain exactly 10 pseudo queries.")
}

if (meta_challenge[
  is_reference %in% TRUE & macroregion_3 == "Africa",
  .N
] != 105L) {
  stop("Expected 105 African training references.")
}

if (meta_challenge[
  is_reference %in% TRUE & macroregion_3 == "Asia",
  .N
] != 215L) {
  stop("Expected 215 Asian training references.")
}

challenge_meta_path <- file.path(challenge_input_dir, "metadata_clean.tsv")
fwrite(meta_challenge, challenge_meta_path, sep = "\t")

# =============================================================================
# Run the unchanged Step 04 assignment workflow with all parameters stated
# explicitly.  The original Step 04 Excel writer does not support a branch
# with zero rows.  This can legitimately occur in a challenge analysis when
# every query is assigned to the same macroregion.  The supplementary validation
# therefore accepts a nonzero Step 04 exit only when the complete tabular
# assignment output exists and contains exactly the expected challenge IDs.
# It can also reuse such validated output after an interrupted reporting stage.
# =============================================================================

rscript_bin <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

if (!file.exists(rscript_bin)) {
  stop("Cannot find Rscript at: ", rscript_bin)
}

step4_args <- c(
  shQuote(opt$step4_script),
  "--qc_dir", shQuote(challenge_input_dir),
  "--step3_dir", shQuote(opt$step3_dir),
  "--out_dir", shQuote(step4_out_dir),
  "--out_xlsx", shQuote("other_regions_closed_set_assignment.xlsx"),
  "--query_is_reference_value", "FALSE",
  "--min_snps_query_macroregion", "200",
  "--min_snps_query_subregion", "300",
  "--pseudocount", "0.5",
  "--epsilon", "0.02",
  "--high_posterior", "0.95",
  "--moderate_posterior", "0.85",
  "--min_gap_macroregion", "0.20",
  "--min_gap_subregion", "0.25",
  "--base_K", "1,2,3,4,5,10,20,50,100,200,500",
  "--extra_K", "1000,2000,5000,10000",
  "--auto_extend_K", "TRUE",
  "--include_all_snps_K", "TRUE",
  "--max_K_points", "30",
  "--min_agreement", "0.90",
  "--min_K_available", "6",
  "--tail_fraction", "0.50",
  "--min_tail_agreement", "1.00"
)

old_wd <- getwd()
step4_final_path <- file.path(step4_out_dir, "final_assignment.tsv")
step4_xlsx_path <- file.path(
  step4_out_dir,
  "other_regions_closed_set_assignment.xlsx"
)

valid_step4_tabular_output <- function(path, expected_ids) {
  if (!file.exists(path)) return(FALSE)

  x <- tryCatch(
    fread(path),
    error = function(e) NULL
  )

  if (is.null(x) || !("sample_id" %in% names(x))) return(FALSE)

  ids <- clean_text(x$sample_id)

  nrow(x) == length(expected_ids) &&
    uniqueN(ids) == length(expected_ids) &&
    setequal(ids, expected_ids)
}

step4_status <- NA_integer_
step4_execution_mode <- NA_character_

if (valid_step4_tabular_output(step4_final_path, challenge_ids)) {
  step4_execution_mode <- "reused_valid_existing_tabular_output"
  message(
    "Reusing validated Step 04 tabular output from the preceding run: ",
    step4_final_path
  )
} else {
  setwd(repo_root)

  step4_status <- tryCatch(
    system2(rscript_bin, args = step4_args),
    finally = setwd(old_wd)
  )

  if (identical(as.integer(step4_status), 0L)) {
    step4_execution_mode <- "new_step04_run_completed"
  } else if (valid_step4_tabular_output(step4_final_path, challenge_ids)) {
    step4_execution_mode <- paste0(
      "new_step04_tabular_output_recovered_after_reporting_exit_",
      as.integer(step4_status)
    )

    warning(
      "Step 04 completed the assignment tables but its optional Excel writer ",
      "stopped because one conditional branch contained zero rows. ",
      "The validated assignment table will be used for the challenge summary."
    )
  } else {
    stop(
      "The Step 04 challenge run failed with status ", step4_status,
      " and did not produce a complete assignment table."
    )
  }
}

if (!file.exists(step4_final_path)) {
  stop("Step 04 did not produce final_assignment.tsv.")
}

# =============================================================================
# Create reader facing challenge summaries
# =============================================================================

assignment <- fread(step4_final_path)

if (nrow(assignment) != 10L || uniqueN(assignment$sample_id) != 10L) {
  stop(
    "Expected 10 unique Step 04 challenge assignments but found ",
    nrow(assignment), " rows and ", uniqueN(assignment$sample_id),
    " unique identifiers."
  )
}

challenge <- merge(
  challenge_truth,
  assignment,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("_truth", "_step4")
)

if (challenge[is.na(macroregion_confidence), .N] > 0L) {
  stop("At least one challenge sample has no macroregion result.")
}

challenge[, closed_set_macroregion_affinity := macroregion_call]
challenge[, macroregion_affinity_accepted :=
  macroregion_confidence %in% c("High", "Moderate")]
challenge[, flagged_uncertain := !macroregion_affinity_accepted]
challenge[, out_of_scope_not_flagged := macroregion_affinity_accepted]

challenge[
  macroregion_affinity_accepted %in% TRUE,
  accepted_macroregion_affinity := macroregion_call
]

challenge[
  macroregion_affinity_accepted %in% FALSE,
  accepted_macroregion_affinity := NA_character_
]

challenge[, subregion_affinity_accepted :=
  macroregion_affinity_accepted %in% TRUE &
    subregion_confidence %in% c("High", "Moderate")]
challenge[, closed_set_subregion_affinity := subregion_call]

challenge[, reporting_level := fifelse(
  !macroregion_affinity_accepted,
  "Uncertain",
  fifelse(subregion_affinity_accepted, "Subregion", "Macroregion")
)]

challenge[, reported_affinity := fifelse(
  reporting_level == "Subregion",
  closed_set_subregion_affinity,
  fifelse(
    reporting_level == "Macroregion",
    macroregion_call,
    "Uncertain"
  )
)]

challenge[, challenge_interpretation := fifelse(
  flagged_uncertain,
  "Flagged uncertain by existing confidence rules",
  paste0(
    "Accepted ", macroregion_call, " genomic affinity; collection ",
    "population absent from candidate classes"
  )
)]

setorder(challenge, collection_site, sample_id)

challenge_columns <- c(
  "sample_id",
  "withheld_collection_class",
  "collection_subregion_label",
  "collection_country",
  "collection_site",
  "collection_population",
  "closed_set_macroregion_affinity",
  "macroregion_affinity_accepted",
  "accepted_macroregion_affinity",
  "flagged_uncertain",
  "out_of_scope_not_flagged",
  "macroregion_confidence",
  "macroregion_reason",
  "macroregion_maxK_post",
  "macroregion_maxK_gap",
  "macroregion_agreement",
  "macroregion_tail_agreement",
  "macroregion_stable_from_K",
  "macroregion_nK",
  "macroregion_maxK_nsnps",
  "closed_set_subregion_affinity",
  "subregion_confidence",
  "subregion_reason",
  "subregion_affinity_accepted",
  "reporting_level",
  "reported_affinity",
  "challenge_interpretation"
)

missing_output_columns <- setdiff(challenge_columns, names(challenge))

if (length(missing_output_columns) > 0L) {
  stop(
    "Challenge output is missing column(s): ",
    paste(missing_output_columns, collapse = ", ")
  )
}

challenge_summary <- challenge[, ..challenge_columns]

fwrite(
  challenge_summary,
  file.path(opt$out_dir, "other_regions_challenge_results.tsv"),
  sep = "\t"
)

saveRDS(
  challenge_summary,
  file.path(opt$out_dir, "other_regions_challenge_results.rds")
)

location_summary <- challenge_summary[
  ,
  .(
    n_samples = .N,
    n_flagged_uncertain = sum(flagged_uncertain, na.rm = TRUE),
    flagged_uncertain_rate = mean(flagged_uncertain, na.rm = TRUE),
    n_out_of_scope_not_flagged = sum(out_of_scope_not_flagged, na.rm = TRUE),
    out_of_scope_not_flagged_rate = mean(
      out_of_scope_not_flagged,
      na.rm = TRUE
    ),
    n_closed_set_africa_calls = sum(
      closed_set_macroregion_affinity == "Africa",
      na.rm = TRUE
    ),
    n_closed_set_asia_calls = sum(
      closed_set_macroregion_affinity == "Asia",
      na.rm = TRUE
    ),
    n_accepted_africa_affinities = sum(
      accepted_macroregion_affinity == "Africa",
      na.rm = TRUE
    ),
    n_accepted_asia_affinities = sum(
      accepted_macroregion_affinity == "Asia",
      na.rm = TRUE
    ),
    median_maxK_posterior = safe_median(macroregion_maxK_post),
    median_maxK_gap = safe_median(macroregion_maxK_gap),
    median_K_agreement = safe_median(macroregion_agreement)
  ),
  by = .(
    collection_country,
    collection_site
  )
]

overall_summary <- challenge_summary[
  ,
  .(
    collection_country = "All",
    collection_site = "All",
    n_samples = .N,
    n_flagged_uncertain = sum(flagged_uncertain, na.rm = TRUE),
    flagged_uncertain_rate = mean(flagged_uncertain, na.rm = TRUE),
    n_out_of_scope_not_flagged = sum(out_of_scope_not_flagged, na.rm = TRUE),
    out_of_scope_not_flagged_rate = mean(
      out_of_scope_not_flagged,
      na.rm = TRUE
    ),
    n_closed_set_africa_calls = sum(
      closed_set_macroregion_affinity == "Africa",
      na.rm = TRUE
    ),
    n_closed_set_asia_calls = sum(
      closed_set_macroregion_affinity == "Asia",
      na.rm = TRUE
    ),
    n_accepted_africa_affinities = sum(
      accepted_macroregion_affinity == "Africa",
      na.rm = TRUE
    ),
    n_accepted_asia_affinities = sum(
      accepted_macroregion_affinity == "Asia",
      na.rm = TRUE
    ),
    median_maxK_posterior = safe_median(macroregion_maxK_post),
    median_maxK_gap = safe_median(macroregion_maxK_gap),
    median_K_agreement = safe_median(macroregion_agreement)
  )
]

location_summary <- rbindlist(
  list(location_summary, overall_summary),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  location_summary,
  file.path(opt$out_dir, "other_regions_challenge_by_location.tsv"),
  sep = "\t"
)

saveRDS(
  location_summary,
  file.path(opt$out_dir, "other_regions_challenge_by_location.rds")
)

run_info <- list(
  analysis = "out_of_scope_geographic_challenge",
  purpose = paste(
    "Test the behaviour of a closed Africa versus Asia classifier when",
    "collection populations are absent from the candidate set"
  ),
  withheld_collection_class = "Others",
  challenge_locations = c("Hawaii", "Papua_New_Guinea"),
  candidate_training_classes = c("Africa", "Asia"),
  n_challenge_samples = nrow(challenge_truth),
  challenge_sample_ids = challenge_ids,
  n_africa_training_references = 105L,
  n_asia_training_references = 215L,
  feature_panels = "existing corrected Step 02 fixed panels",
  assignment_workflow = "unchanged Step 04 hierarchy and thresholds",
  terminology_version = "collection_location_vs_source_affinity_v2",
  primary_endpoint = paste(
    "rate at which geographically withheld samples are flagged uncertain",
    "by the closed reference confidence rules"
  ),
  interpretation = paste(
    "A reported class is genomic affinity within the available candidate set",
    "and is not evidence of direct collection origin when the collection",
    "population is absent from the model"
  ),
  original_metadata = opt$metadata,
  challenge_metadata = challenge_meta_path,
  step3_dir = opt$step3_dir,
  step4_script = opt$step4_script,
  step4_output = step4_out_dir,
  step4_execution_mode = step4_execution_mode,
  step4_exit_status = step4_status,
  step4_workbook_created = file.exists(step4_xlsx_path),
  step4_reporting_note = paste(
    "The original Step 04 Excel writer cannot write a conditional branch",
    "with zero columns. Validated TSV and RDS assignments remain complete",
    "and are the authoritative inputs to this challenge summary."
  ),
  challenge_summary = file.path(
    opt$out_dir,
    "other_regions_challenge_results.tsv"
  ),
  location_summary = file.path(
    opt$out_dir,
    "other_regions_challenge_by_location.tsv"
  ),
  step4_parameters = list(
    min_snps_query_macroregion = 200L,
    min_snps_query_subregion = 300L,
    pseudocount = 0.5,
    epsilon = 0.02,
    high_posterior = 0.95,
    moderate_posterior = 0.85,
    min_gap_macroregion = 0.20,
    min_gap_subregion = 0.25,
    base_K = c(1L, 2L, 3L, 4L, 5L, 10L, 20L, 50L, 100L, 200L, 500L),
    extra_K = c(1000L, 2000L, 5000L, 10000L),
    auto_extend_K = TRUE,
    include_all_snps_K = TRUE,
    max_K_points = 30L,
    min_agreement = 0.90,
    min_K_available = 6L,
    tail_fraction = 0.50,
    min_tail_agreement = 1.00
  ),
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(
  run_info,
  file.path(opt$out_dir, "challenge_run_info.rds")
)

message("\nDone.")
message(
  "Sample results:   ",
  file.path(opt$out_dir, "other_regions_challenge_results.tsv")
)
message(
  "Location summary: ",
  file.path(opt$out_dir, "other_regions_challenge_by_location.tsv")
)
message("Full Step 04 output: ", step4_out_dir)
message("\nOverall outcome:")
print(overall_summary)
