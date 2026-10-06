#!/usr/bin/env Rscript

# PESTFLY: Reports for plant health authorities
# Format saved Step 04 assignments with six criterion outcomes.
# Read thresholds from the matching parameter record and check saved decisions.
# RF and commodity comparisons remain separate; no model is refitted.
# Run from the repository root:
#   Rscript steps/07_authority_facing_report/run.R
# Inputs, outputs and options: adjacent README.md. Full interpretation: manuscript.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(openxlsx)
})

# =============================================================================
# Repo-root detection
# =============================================================================

find_repo_root <- function(max_up = 10L) {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  cand <- wd
  
  for (i in seq_len(max_up + 1L)) {
    if (dir.exists(file.path(cand, "data")) && dir.exists(file.path(cand, "steps"))) {
      return(cand)
    }
    
    parent <- normalizePath(file.path(cand, ".."), winslash = "/", mustWork = FALSE)
    if (identical(parent, cand)) break
    cand <- parent
  }
  
  wd
}

repo_root <- find_repo_root()

is_abs_path <- function(p) {
  grepl("^(?:[A-Za-z]:[/\\\\]|/)", p)
}

resolve_path <- function(p, root = repo_root) {
  if (is.null(p) || !nzchar(p)) return(p)
  if (is_abs_path(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(root, p), winslash = "/", mustWork = FALSE)
}

# =============================================================================
# CLI
# =============================================================================

option_list <- list(
  make_option(
    "--step4_final",
    type = "character",
    default = "results/04_origin_assignment/final_assignment.tsv",
    help = "Step 04 final assignment table [default %default]"
  ),
  make_option(
    "--step4_params",
    type = "character",
    default = "",
    help = "Step 04 parameter TSV; empty uses step4_params.tsv beside --step4_final"
  ),
  make_option(
    "--step5_dir",
    type = "character",
    default = "results/06_rf_corroboration",
    help = "Step 06 output directory [default %default]"
  ),
  make_option(
    "--step4b_dir",
    type = "character",
    default = "results/05_loo_validation",
    help = "Step 05 output directory [default %default]"
  ),
  make_option(
    "--step6b_dir",
    type = "character",
    default = "",
    help = "Optional internal diagnostic-development Step 09 directory; empty skips reduced-panel sheets [default empty]"
  ),
  make_option(
    "--intercept_origin_file",
    type = "character",
    default = "",
    help = "Optional TSV with sample_id and sample_origin_country [default empty: use built-in FAVV mapping]"
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/07_authority_report",
    help = "Step 07 output directory [default %default]"
  ),
  make_option(
    "--out_xlsx",
    type = "character",
    default = "FINAL_origin_tracing_authority_REPORT.xlsx",
    help = "Authority-facing Excel workbook filename [default %default]"
  ),
  make_option(
    "--include_validation_sheets",
    type = "character",
    default = "TRUE",
    help = "Include validation sheets if available: TRUE/FALSE [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$step4_final <- resolve_path(opt$step4_final)
opt$step4_params <- if (nzchar(opt$step4_params)) resolve_path(opt$step4_params) else
  file.path(dirname(opt$step4_final), "step4_params.tsv")
source(file.path(repo_root, "steps", "07_authority_facing_report", "reporting_criteria.R"))
opt$step5_dir <- resolve_path(opt$step5_dir)
opt$step4b_dir <- resolve_path(opt$step4b_dir)
if (nzchar(opt$step6b_dir)) {
  opt$step6b_dir <- resolve_path(opt$step6b_dir)
}
opt$out_dir <- resolve_path(opt$out_dir)

if (nzchar(opt$intercept_origin_file)) {
  opt$intercept_origin_file <- resolve_path(opt$intercept_origin_file)
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Generic helpers
# =============================================================================

clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("[\u200B-\u200D\uFEFF]", "", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x
}

parse_bool1 <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L) return(as.logical(default)[1])
  
  y <- tolower(clean_text(x[1]))
  
  if (is.na(y) || !nzchar(y)) return(as.logical(default)[1])
  if (y %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (y %in% c("false", "f", "0", "no", "n")) return(FALSE)
  
  stop("Cannot parse boolean value: ", x)
}

safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  setnames(dt, tolower(names(dt)))
  if ("sample_id" %in% names(dt)) {
    dt[, sample_id := clean_text(sample_id)]
  }
  dt
}

safe_sheet <- function(wb, name) {
  nm <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", name)
  if (nchar(nm) > 31) nm <- substr(nm, 1, 31)
  
  existing <- names(wb)
  if (!(nm %in% existing)) return(nm)
  
  i <- 2
  repeat {
    suffix <- paste0("_", i)
    base <- substr(nm, 1, max(1, 31 - nchar(suffix)))
    nm2 <- paste0(base, suffix)
    
    if (!(nm2 %in% existing)) return(nm2)
    
    i <- i + 1
  }
}

first_existing_col <- function(dt, candidates) {
  candidates <- tolower(candidates)
  hit <- intersect(candidates, names(dt))
  if (length(hit) > 0L) hit[1] else NA_character_
}

get_col_or_na <- function(dt, col, type = "character") {
  if (is.na(col) || !(col %in% names(dt))) {
    if (type == "numeric") return(rep(NA_real_, nrow(dt)))
    if (type == "integer") return(rep(NA_integer_, nrow(dt)))
    if (type == "logical") return(rep(NA, nrow(dt)))
    return(rep(NA_character_, nrow(dt)))
  }
  
  x <- dt[[col]]
  
  if (type == "numeric") return(suppressWarnings(as.numeric(x)))
  if (type == "integer") return(suppressWarnings(as.integer(x)))
  if (type == "logical") return(as.logical(x))
  
  clean_text(x)
}

is_accepted_conf <- function(x) {
  tolower(clean_text(x)) %in% c("high", "moderate")
}

is_high_conf <- function(x) {
  tolower(clean_text(x)) == "high"
}

is_moderate_conf <- function(x) {
  tolower(clean_text(x)) == "moderate"
}

same_nonmissing <- function(a, b) {
  aa <- tolower(clean_text(a))
  bb <- tolower(clean_text(b))
  
  ok <- !is.na(aa) & !is.na(bb)
  out <- rep(NA, length(aa))
  out[ok] <- aa[ok] == bb[ok]
  
  out
}

norm_country <- function(x) {
  y <- clean_text(x)
  y <- gsub("Côte d'Ivoire", "Ivory Coast", y, fixed = TRUE)
  y <- gsub("Cote d'Ivoire", "Ivory Coast", y, fixed = TRUE)
  y <- gsub("Côte d’Ivoire", "Ivory Coast", y, fixed = TRUE)
  y <- gsub("Cote d’Ivoire", "Ivory Coast", y, fixed = TRUE)
  y
}

safe_ifelse_text <- function(condition, yes, no) {
  fifelse(condition, yes, no, na = no)
}

# =============================================================================
# Origin-country mapping for intercepted samples
# =============================================================================

built_in_intercept_origin <- function() {
  data.table(
    sample_id = c(
      "FAVV_1", "FAVV_2", "FAVV_3", "FAVV_4", "FAVV_5", "FAVV_6", "FAVV_7",
      "FAVV_11", "FAVV_12", "FAVV_13", "FAVV_14"
    ),
    sample_origin_country = c(
      "Cameroon", "Cameroon", "Ivory Coast", "Ivory Coast", "Senegal", "Ivory Coast", "Senegal",
      "Burkina Faso", "Senegal", "Senegal", "Bangladesh"
    ),
    commodity = "Mango",
    sample_origin_country_source = "Built-in FAVV/FASFC origin metadata from Vanbergen et al. 2025 Table 1"
  )
}

read_or_build_intercept_origin <- function(path) {
  if (!is.null(path) && nzchar(path) && !file.exists(path)) {
    stop("Missing explicitly supplied intercept origin table: ", path)
  }
  if (!is.null(path) && nzchar(path) && file.exists(path)) {
    dt <- fread(path)
    setnames(dt, tolower(names(dt)))
    
    if (!all(c("sample_id", "sample_origin_country") %in% names(dt))) {
      stop("--intercept_origin_file must contain sample_id and sample_origin_country")
    }
    
    dt[, sample_id := clean_text(sample_id)]
    dt[, sample_origin_country := norm_country(sample_origin_country)]
    if (anyNA(dt$sample_id) || anyDuplicated(dt$sample_id)) {
      stop("Intercept origin sample_id values must be populated and unique")
    }
    
    if (!("sample_origin_country_source" %in% names(dt))) {
      dt[, sample_origin_country_source := "User-provided intercept origin table"]
    }
    
    if (!("commodity" %in% names(dt))) {
      dt[, commodity := NA_character_]
    }
    
    return(dt[, .(sample_id, sample_origin_country, commodity, sample_origin_country_source)])
  }
  
  built_in_intercept_origin()
}

country_region_map <- data.table(
  sample_origin_country = c(
    "Cameroon",
    "Ivory Coast",
    "Senegal",
    "Burkina Faso",
    "Bangladesh"
  ),
  expected_macroregion_from_origin_country = c(
    "Africa",
    "Africa",
    "Africa",
    "Africa",
    "Asia"
  ),
  expected_subregion_from_origin_country = c(
    "C_Africa",
    "W_Africa",
    "W_Africa",
    "W_Africa",
    "S_Asia"
  )
)

# =============================================================================
# Vanbergen et al. 2025 macroregion reference interpretation
# =============================================================================

vanbergen_macroregion_table <- data.table(
  sample_id = c(
    "Bd_2023_1", "Bd_2023_2", "Bd_2023_3", "Bd_2023_4", "Bd_2023_5", "Bd_2023_6", "Bd_2023_7",
    "Bd_2024_1", "Bd_2024_2", "Bd_2024_3", "Bd_2024_4",
    "FAVV_1", "FAVV_2", "FAVV_3", "FAVV_4", "FAVV_5", "FAVV_6", "FAVV_7",
    "FAVV_11", "FAVV_12", "FAVV_13", "FAVV_14"
  ),
  vanbergen_macroregion = c(
    "Asia", "Africa", "Asia", "Asia", "Asia", "Asia", "Africa",
    "Asia", "Africa", "Asia", "Asia",
    "Africa", "Africa", "Africa", "Africa", "Africa", "Africa", "Africa",
    "Africa", "Africa", "Africa", "Asia"
  ),
  vanbergen_note = "Macroregion-level interpretation from Vanbergen et al. 2025 Table 2 / discussion"
)

# =============================================================================
# RF reading and standardisation
# =============================================================================

read_rf_all_or_panels <- function(step5_dir) {
  all_file <- file.path(step5_dir, "RF_all_panels_predictions_queries.tsv")
  
  if (file.exists(all_file)) {
    rf <- fread(all_file)
    setnames(rf, tolower(names(rf)))
    if (!("sample_id" %in% names(rf))) stop("RF_all_panels_predictions_queries.tsv missing sample_id")
    if (!("panel_id" %in% names(rf))) stop("RF_all_panels_predictions_queries.tsv missing panel_id")
    rf[, sample_id := clean_text(sample_id)]
    rf[, panel_id := clean_text(panel_id)]
    return(rf)
  }
  
  panel_files <- c(
    file.path(step5_dir, "RF_P1_macroregion_africa_vs_asia_predictions_queries.tsv"),
    file.path(step5_dir, "RF_P2_subregion_within_africa_predictions_queries.tsv"),
    file.path(step5_dir, "RF_P3_subregion_within_asia_predictions_queries.tsv")
  )
  
  panel_ids <- c(
    "P1_macroregion_africa_vs_asia",
    "P2_subregion_within_africa",
    "P3_subregion_within_asia"
  )
  
  out <- list()
  
  for (i in seq_along(panel_files)) {
    f <- panel_files[i]
    if (!file.exists(f)) next
    
    dt <- fread(f)
    setnames(dt, tolower(names(dt)))
    
    if (!("sample_id" %in% names(dt))) {
      stop("RF panel file missing sample_id: ", f)
    }
    
    dt[, sample_id := clean_text(sample_id)]
    
    if (!("panel_id" %in% names(dt))) {
      dt[, panel_id := panel_ids[i]]
    }
    
    out[[length(out) + 1L]] <- dt
  }
  
  if (length(out) == 0L) return(NULL)
  
  rbindlist(out, fill = TRUE)
}

standardise_rf_by_panel <- function(rf) {
  if (is.null(rf) || nrow(rf) == 0L) return(NULL)
  
  call_col <- first_existing_col(rf, c("top_group", "rf_call", "call", "predicted_group"))
  prob_col <- first_existing_col(rf, c("top_prob", "rf_top_prob", "prob", "top_posterior"))
  gap_col <- first_existing_col(rf, c("prob_gap", "rf_gap", "posterior_gap", "gap"))
  conf_col <- first_existing_col(rf, c("confidence", "rf_confidence", "final_confidence"))
  status_col <- first_existing_col(rf, c("status", "rf_note", "note"))
  vs_col <- first_existing_col(rf, c("rf_vs_step4"))
  
  rf_std <- data.table(
    sample_id = rf$sample_id,
    panel_id = rf$panel_id,
    rf_call = get_col_or_na(rf, call_col, "character"),
    rf_top_prob = get_col_or_na(rf, prob_col, "numeric"),
    rf_gap = get_col_or_na(rf, gap_col, "numeric"),
    rf_confidence = get_col_or_na(rf, conf_col, "character"),
    rf_status = get_col_or_na(rf, status_col, "character"),
    rf_vs_step4 = get_col_or_na(rf, vs_col, "character")
  )
  
  rf_std[, panel_class := fifelse(
    grepl("^P1_", panel_id), "macroregion",
    fifelse(grepl("^P2_", panel_id), "africa_subregion",
            fifelse(grepl("^P3_", panel_id), "asia_subregion", "other")
    )
  )]
  
  rf_std
}

wide_rf <- function(rf_std) {
  if (is.null(rf_std) || nrow(rf_std) == 0L) {
    return(data.table(sample_id = character()))
  }
  
  make_one <- function(class_name, prefix) {
    x <- rf_std[panel_class == class_name]
    
    if (nrow(x) == 0L) {
      return(data.table(sample_id = character()))
    }
    
    x <- x[, .SD[1], by = sample_id]
    
    out <- x[
      ,
      .(
        sample_id,
        call = rf_call,
        top_prob = rf_top_prob,
        gap = rf_gap,
        confidence = rf_confidence,
        status = rf_status,
        rf_vs_step4 = rf_vs_step4
      )
    ]
    
    setnames(
      out,
      old = c("call", "top_prob", "gap", "confidence", "status", "rf_vs_step4"),
      new = paste0(prefix, c("call", "top_prob", "gap", "confidence", "status", "rf_vs_step4"))
    )
    
    out
  }
  
  m <- make_one("macroregion", "rf_macro_")
  a <- make_one("africa_subregion", "rf_africa_sub_")
  s <- make_one("asia_subregion", "rf_asia_sub_")
  
  Reduce(
    function(x, y) merge(x, y, by = "sample_id", all = TRUE),
    list(m, a, s)
  )
}

# =============================================================================
# Load Step 04 primary results
# =============================================================================

include_validation_sheets <- parse_bool1(opt$include_validation_sheets, default = TRUE)

if (!file.exists(opt$step4_final)) {
  stop("Missing Step 04 final assignment table: ", opt$step4_final)
}

s4 <- fread(opt$step4_final)
setnames(s4, tolower(names(s4)))

if (!("sample_id" %in% names(s4))) {
  stop("Step 04 final assignment table missing sample_id")
}

s4[, sample_id := clean_text(sample_id)]

macro_call_col <- first_existing_col(s4, c("macroregion_call", "final_macroregion", "macroregion"))
macro_conf_col <- first_existing_col(s4, c("macroregion_confidence", "macroregion_final_confidence", "confidence"))
macro_reason_col <- first_existing_col(s4, c("macroregion_reason", "macroregion_final_reason"))
macro_post_col <- first_existing_col(s4, c("macroregion_maxk_post", "macroregion_posterior", "maxk_top_posterior"))
macro_gap_col <- first_existing_col(s4, c("macroregion_maxk_gap", "macroregion_gap", "maxk_gap"))
macro_agree_col <- first_existing_col(s4, c("macroregion_agreement", "agreement_frac"))
macro_tail_col <- first_existing_col(s4, c("macroregion_tail_agreement", "tail_agreement_frac"))
macro_nk_col <- first_existing_col(s4, c("macroregion_nk_usable_for_stability", "macroregion_nk", "n_k_usable_for_stability"))

sub_call_col <- first_existing_col(s4, c("subregion_call", "final_subregion", "subregion_final"))
sub_conf_col <- first_existing_col(s4, c("subregion_confidence", "subregion_final_confidence"))
sub_reason_col <- first_existing_col(s4, c("subregion_reason", "subregion_final_reason"))
sub_panel_col <- first_existing_col(s4, c("subregion_panel"))
sub_post_col <- first_existing_col(s4, c("subregion_maxk_post", "subregion_posterior"))
sub_gap_col <- first_existing_col(s4, c("subregion_maxk_gap", "subregion_gap"))
sub_agree_col <- first_existing_col(s4, c("subregion_agreement"))
sub_tail_col <- first_existing_col(s4, c("subregion_tail_agreement"))
sub_nk_col <- first_existing_col(s4, c("subregion_nk_usable_for_stability", "subregion_nk"))

if (is.na(macro_call_col) || is.na(macro_conf_col)) {
  stop(
    "Step 04 final table does not contain recognizable macroregion call/confidence columns.\n",
    "Found columns: ", paste(names(s4), collapse = ", ")
  )
}

reporting_parameters <- read_reporting_parameters(opt$step4_params)
criterion_report <- assess_reporting_criteria(s4, reporting_parameters)
criteria_detail <- as.data.table(criterion_report$detail)
criteria_summary <- as.data.table(criterion_report$summary)
reporting_thresholds <- as.data.table(criterion_report$thresholds)

extended <- copy(s4)

extended[, step4_macroregion_call := get_col_or_na(s4, macro_call_col, "character")]
extended[, step4_macroregion_confidence := get_col_or_na(s4, macro_conf_col, "character")]
extended[, step4_macroregion_reason := get_col_or_na(s4, macro_reason_col, "character")]
extended[, step4_macroregion_posterior := get_col_or_na(s4, macro_post_col, "numeric")]
extended[, step4_macroregion_gap := get_col_or_na(s4, macro_gap_col, "numeric")]
extended[, step4_macroregion_agreement := get_col_or_na(s4, macro_agree_col, "numeric")]
extended[, step4_macroregion_tail_agreement := get_col_or_na(s4, macro_tail_col, "numeric")]
extended[, step4_macroregion_nK_usable := get_col_or_na(s4, macro_nk_col, "integer")]

extended[, step4_subregion_call := get_col_or_na(s4, sub_call_col, "character")]
extended[, step4_subregion_confidence := get_col_or_na(s4, sub_conf_col, "character")]
extended[, step4_subregion_reason := get_col_or_na(s4, sub_reason_col, "character")]
extended[, step4_subregion_panel := get_col_or_na(s4, sub_panel_col, "character")]
extended[, step4_subregion_posterior := get_col_or_na(s4, sub_post_col, "numeric")]
extended[, step4_subregion_gap := get_col_or_na(s4, sub_gap_col, "numeric")]
extended[, step4_subregion_agreement := get_col_or_na(s4, sub_agree_col, "numeric")]
extended[, step4_subregion_tail_agreement := get_col_or_na(s4, sub_tail_col, "numeric")]
extended[, step4_subregion_nK_usable := get_col_or_na(s4, sub_nk_col, "integer")]

# =============================================================================
# Add sample-origin country metadata
# =============================================================================

intercept_origin <- read_or_build_intercept_origin(opt$intercept_origin_file)

extended <- merge(
  extended,
  intercept_origin,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

extended[, sample_origin_country := norm_country(sample_origin_country)]
extended[is.na(sample_origin_country), sample_origin_country := "Unknown"]
extended[is.na(sample_origin_country_source), sample_origin_country_source := "Unknown for trapped adults or unavailable metadata"]
extended[is.na(commodity), commodity := NA_character_]

extended <- merge(
  extended,
  country_region_map,
  by = "sample_origin_country",
  all.x = TRUE,
  sort = FALSE
)

extended[
  sample_origin_country == "Unknown",
  `:=`(
    expected_macroregion_from_origin_country = NA_character_,
    expected_subregion_from_origin_country = NA_character_
  )
]

# =============================================================================
# Add Vanbergen et al. 2025 macroregion consistency
# =============================================================================

extended <- merge(
  extended,
  vanbergen_macroregion_table,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

extended[, vanbergen_macroregion_consistency := fifelse(
  is.na(vanbergen_macroregion),
  "No Vanbergen comparison available",
  fifelse(
    tolower(step4_macroregion_call) == tolower(vanbergen_macroregion),
    "Consistent with Vanbergen et al. 2025 macroregion",
    "Different from Vanbergen et al. 2025 macroregion"
  )
)]

# =============================================================================
# Merge Step 06 RF corroboration
# =============================================================================

rf_raw <- read_rf_all_or_panels(opt$step5_dir)
rf_std <- standardise_rf_by_panel(rf_raw)
rf_wide <- wide_rf(rf_std)

if (nrow(rf_wide) > 0L) {
  extended <- merge(extended, rf_wide, by = "sample_id", all.x = TRUE, sort = FALSE)
} else {
  extended[, `:=`(
    rf_macro_call = NA_character_,
    rf_macro_top_prob = NA_real_,
    rf_macro_gap = NA_real_,
    rf_macro_confidence = NA_character_,
    rf_macro_status = NA_character_,
    rf_macro_rf_vs_step4 = NA_character_,
    rf_africa_sub_call = NA_character_,
    rf_africa_sub_top_prob = NA_real_,
    rf_africa_sub_gap = NA_real_,
    rf_africa_sub_confidence = NA_character_,
    rf_africa_sub_status = NA_character_,
    rf_africa_sub_rf_vs_step4 = NA_character_,
    rf_asia_sub_call = NA_character_,
    rf_asia_sub_top_prob = NA_real_,
    rf_asia_sub_gap = NA_real_,
    rf_asia_sub_confidence = NA_character_,
    rf_asia_sub_status = NA_character_,
    rf_asia_sub_rf_vs_step4 = NA_character_
  )]
}

# A missing optional panel must not prevent primary criterion reporting.
rf_expected <- c("call", "top_prob", "gap", "confidence", "status", "rf_vs_step4")
for (prefix in c("rf_macro_", "rf_africa_sub_", "rf_asia_sub_")) {
  for (suffix in rf_expected) {
    key <- paste0(prefix, suffix)
    if (!(key %in% names(extended))) {
      extended[, (key) := if (suffix %in% c("top_prob", "gap")) NA_real_ else NA_character_]
    }
  }
}

extended[, rf_relevant_subregion_call := fifelse(
  tolower(step4_macroregion_call) == "africa",
  rf_africa_sub_call,
  fifelse(tolower(step4_macroregion_call) == "asia", rf_asia_sub_call, NA_character_)
)]

extended[, rf_relevant_subregion_top_prob := fifelse(
  tolower(step4_macroregion_call) == "africa",
  rf_africa_sub_top_prob,
  fifelse(tolower(step4_macroregion_call) == "asia", rf_asia_sub_top_prob, NA_real_)
)]

extended[, rf_relevant_subregion_gap := fifelse(
  tolower(step4_macroregion_call) == "africa",
  rf_africa_sub_gap,
  fifelse(tolower(step4_macroregion_call) == "asia", rf_asia_sub_gap, NA_real_)
)]

extended[, rf_relevant_subregion_confidence := fifelse(
  tolower(step4_macroregion_call) == "africa",
  rf_africa_sub_confidence,
  fifelse(tolower(step4_macroregion_call) == "asia", rf_asia_sub_confidence, NA_character_)
)]

extended[, rf_macro_agrees_step4 := same_nonmissing(step4_macroregion_call, rf_macro_call)]
extended[, rf_subregion_agrees_step4 := same_nonmissing(step4_subregion_call, rf_relevant_subregion_call)]

extended[, rf_macro_strong_discordance := (
  rf_macro_agrees_step4 %in% FALSE &
    is_accepted_conf(rf_macro_confidence)
)]

extended[, rf_subregion_strong_discordance := (
  rf_subregion_agrees_step4 %in% FALSE &
    is_accepted_conf(rf_relevant_subregion_confidence)
)]

extended[, rf_corroboration := fifelse(
  is.na(rf_macro_call) & is.na(rf_relevant_subregion_call),
  "No RF result",
  fifelse(
    rf_macro_strong_discordance %in% TRUE,
    "RF caution: macroregion discordance",
    fifelse(
      rf_subregion_strong_discordance %in% TRUE,
      "RF caution: subregion discordance",
      fifelse(
        rf_subregion_agrees_step4 %in% TRUE & is_accepted_conf(rf_relevant_subregion_confidence),
        "RF supports subregion",
        fifelse(
          rf_macro_agrees_step4 %in% TRUE & is_accepted_conf(rf_macro_confidence),
          "RF supports macroregion",
          fifelse(
            rf_macro_agrees_step4 %in% TRUE,
            "RF weakly supports macroregion",
            "RF inconclusive"
          )
        )
      )
    )
  )
)]

extended[, rf_corroboration_detail := fifelse(
  rf_corroboration == "RF supports subregion",
  paste0(
    "RF relevant subregion call agrees with Step 04 subregion call (RF=", rf_relevant_subregion_call,
    ", Step4=", step4_subregion_call,
    ") with accepted RF confidence (", rf_relevant_subregion_confidence, ")."
  ),
  fifelse(
    rf_corroboration == "RF supports macroregion",
    paste0(
      "RF macroregion call agrees with Step 04 macroregion call (RF=", rf_macro_call,
      ", Step4=", step4_macroregion_call,
      ") with accepted RF confidence (", rf_macro_confidence, ")."
    ),
    fifelse(
      rf_corroboration == "RF weakly supports macroregion",
      paste0(
        "RF macroregion top class agrees with Step 04 macroregion call (RF=", rf_macro_call,
        ", Step4=", step4_macroregion_call,
        "), but RF confidence is ", rf_macro_confidence,
        ". This is directional agreement only, not strong independent corroboration. ",
        "If subregion is not reported, this is due to Step 04 subregion filtering, e.g. ",
        ifelse(!is.na(step4_subregion_reason), step4_subregion_reason, "no accepted Step 04 subregion call"),
        ", not because RF caused the subregion to be withheld."
      ),
      fifelse(
        rf_corroboration == "RF caution: macroregion discordance",
        paste0(
          "RF macroregion call conflicts with Step 04 macroregion call with accepted RF confidence: RF=",
          rf_macro_call, ", Step4=", step4_macroregion_call, "."
        ),
        fifelse(
          rf_corroboration == "RF caution: subregion discordance",
          paste0(
            "RF relevant subregion call conflicts with Step 04 subregion call with accepted RF confidence: RF=",
            rf_relevant_subregion_call, ", Step4=", step4_subregion_call, "."
          ),
          fifelse(
            rf_corroboration == "No RF result",
            "No RF result was available for this sample/panel.",
            "RF did not provide a clear accepted-confidence corroboration or contradiction."
          )
        )
      )
    )
  )
)]

# =============================================================================
# Step 07 conservative reporting decision
# =============================================================================

for (res in c("macroregion", "subregion")) {
  one <- criteria_summary[resolution == res,
    .(sample_id, status, criteria_met, criteria_unmet, unmet_criteria, step4_confidence_consistent)]
  setnames(one, setdiff(names(one), "sample_id"),
           paste0(res, "_", setdiff(names(one), "sample_id")))
  extended <- merge(extended, one, by = "sample_id", all.x = TRUE, sort = FALSE)
}

extended[, macro_accepted := macroregion_status == "PASS"]
extended[, subregion_accepted := subregion_status == "PASS"]
extended[, reported_level := fifelse(
  subregion_accepted, "subregion", fifelse(macro_accepted, "macroregion", "uncertain"))]
extended[, reported_origin := fifelse(
  subregion_accepted, step4_subregion_call, fifelse(macro_accepted, step4_macroregion_call, "Uncertain"))]
extended[, reported_confidence := fifelse(
  subregion_accepted, step4_subregion_confidence,
  fifelse(macro_accepted, step4_macroregion_confidence, "Uncertain"))]
extended[, reporting_basis := fifelse(
  reported_level == "subregion", "Macroregion PASS 6/6 and subregion PASS 6/6.",
  fifelse(reported_level == "macroregion",
    paste0("Macroregion PASS 6/6. Subregion FAIL ", subregion_criteria_met,
           "/6. Unmet criteria: ", subregion_unmet_criteria, "."),
    paste0("Macroregion FAIL ", macroregion_criteria_met, "/6. Unmet criteria: ",
           macroregion_unmet_criteria, ". Subregion NOT_EVALUATED.")))]
extended[, subregion_withheld_explanation := fifelse(
  subregion_status == "PASS", "None. All six subregion criteria are met.",
  fifelse(subregion_status == "NOT_EVALUATED", "Subregion was not evaluated because macroregion failed.",
          paste0(subregion_criteria_unmet, " of six subregion criteria unmet: ", subregion_unmet_criteria, ".")))]
extended[, action_note := fifelse(
  reported_level == "uncertain", "Report unresolved genomic assignment.",
  fifelse(reported_level == "macroregion", "Report macroregion affinity only.",
          "Report conditional subregion affinity."))]
extended[, report_statement := paste0(sample_id, ": ",
  fifelse(reported_level == "uncertain", "genomic assignment unresolved. ",
          paste0("genomic affinity to ", reported_origin, " (", reported_level, "). ")),
  reporting_basis)]

# =============================================================================
# Commodity-origin mismatch logic
# =============================================================================

extended[, commodity_origin_vs_assignment := fifelse(
  sample_origin_country == "Unknown" | is.na(sample_origin_country),
  "Not assessable: sample-origin country unknown",
  fifelse(
    reported_level == "uncertain",
    "Not assessable: genetic assignment uncertain",
    fifelse(
      reported_level == "macroregion" &
        !is.na(expected_macroregion_from_origin_country) &
        tolower(reported_origin) == tolower(expected_macroregion_from_origin_country),
      "Consistent at macroregion level; subregion not reported",
      fifelse(
        reported_level == "macroregion" &
          !is.na(expected_macroregion_from_origin_country) &
          tolower(reported_origin) != tolower(expected_macroregion_from_origin_country),
        "Mismatch at macroregion level",
        fifelse(
          reported_level == "subregion" &
            !is.na(expected_subregion_from_origin_country) &
            tolower(reported_origin) == tolower(expected_subregion_from_origin_country),
          "Consistent at subregion level",
          fifelse(
            reported_level == "subregion" &
              !is.na(expected_macroregion_from_origin_country) &
              tolower(step4_macroregion_call) != tolower(expected_macroregion_from_origin_country),
            "Mismatch at macroregion level",
            fifelse(
              reported_level == "subregion" &
                !is.na(expected_subregion_from_origin_country) &
                tolower(reported_origin) != tolower(expected_subregion_from_origin_country),
              "Mismatch at subregion level",
              "Not assessable"
            )
          )
        )
      )
    )
  )
)]

extended[, commodity_origin_mismatch_detail := fifelse(
  sample_origin_country == "Unknown" | is.na(sample_origin_country),
  "No known commodity-origin country is available for this sample. This is expected for trapped adults.",
  fifelse(
    reported_level == "uncertain",
    paste0(
      "The commodity-origin country is ", sample_origin_country,
      ", but the genetic assignment is uncertain; mismatch cannot be assessed."
    ),
    fifelse(
      commodity_origin_vs_assignment == "Consistent at subregion level",
      paste0(
        "Known commodity-origin country is ", sample_origin_country,
        ", expected subregion is ", expected_subregion_from_origin_country,
        ", and the reported genetic origin is ", reported_origin, "."
      ),
      fifelse(
        commodity_origin_vs_assignment == "Consistent at macroregion level; subregion not reported",
        paste0(
          "Known commodity-origin country is ", sample_origin_country,
          ", expected macroregion is ", expected_macroregion_from_origin_country,
          ", and the reported genetic origin is ", reported_origin,
          ". Subregion is not reported because Step 04 subregion confidence/stability was insufficient",
          ifelse(!is.na(step4_subregion_reason), paste0(" (", step4_subregion_reason, ")."), ".")
        ),
        fifelse(
          commodity_origin_vs_assignment == "Mismatch at subregion level",
          paste0(
            "Potential subregional mismatch: sample came from ", sample_origin_country,
            " (expected ", expected_subregion_from_origin_country,
            "), but the reported genetic origin is ", reported_origin,
            ". This may reflect limited subregional resolution, incomplete reference coverage, gene flow/admixture, or a true commodity/origin discrepancy."
          ),
          fifelse(
            commodity_origin_vs_assignment == "Mismatch at macroregion level",
            paste0(
              "Potential macroregional mismatch: sample came from ", sample_origin_country,
              " (expected ", expected_macroregion_from_origin_country,
              "), but the reported genetic macroregion is ", step4_macroregion_call,
              ". This should be checked carefully."
            ),
            "Mismatch assessment not available."
          )
        )
      )
    )
  )
)]

# =============================================================================
# Simplified authority-facing table
# =============================================================================

meta_cols <- intersect(
  c("sample_id", "country", "population", "macroregion_3", "subregion"),
  names(extended)
)

simplified <- extended[
  ,
  .(
    sample_id,
    sample_origin_country,
    sample_origin_country_source,
    commodity,
    expected_macroregion_from_origin_country,
    expected_subregion_from_origin_country,
    reported_origin,
    reported_level,
    macroregion_status,
    macroregion_criteria_met,
    macroregion_criteria_unmet,
    macroregion_unmet_criteria,
    subregion_status,
    subregion_criteria_met,
    subregion_criteria_unmet,
    subregion_unmet_criteria,
    reported_confidence,
    commodity_origin_vs_assignment,
    commodity_origin_mismatch_detail,
    rf_corroboration,
    rf_corroboration_detail,
    subregion_withheld_explanation,
    action_note,
    report_statement,
    vanbergen_macroregion,
    vanbergen_macroregion_consistency,
    step4_macroregion_call,
    step4_macroregion_confidence,
    step4_macroregion_reason,
    step4_subregion_call,
    step4_subregion_confidence,
    step4_subregion_reason,
    rf_macro_call,
    rf_macro_confidence,
    rf_relevant_subregion_call,
    rf_relevant_subregion_confidence
  )
]

# Add non-reference metadata columns except is_reference.
for (mc in setdiff(meta_cols, "sample_id")) {
  if (!(mc %in% names(simplified))) {
    simplified[, (mc) := extended[[mc]]]
  }
}

front_cols <- c(
  "sample_id",
  "reported_origin",
  "reported_level",
  "macroregion_status",
  "macroregion_criteria_met",
  "subregion_status",
  "subregion_criteria_met",
  "reported_confidence",
  "macroregion_unmet_criteria",
  "subregion_unmet_criteria",
  "rf_corroboration",
  "commodity_origin_vs_assignment",
  setdiff(meta_cols, "sample_id"),
  "sample_origin_country",
  "sample_origin_country_source",
  "commodity",
  "expected_macroregion_from_origin_country",
  "expected_subregion_from_origin_country",
  "reported_origin",
  "reported_level",
  "macroregion_status",
  "macroregion_criteria_met",
  "macroregion_criteria_unmet",
  "macroregion_unmet_criteria",
  "subregion_status",
  "subregion_criteria_met",
  "subregion_criteria_unmet",
  "subregion_unmet_criteria",
  "reported_confidence",
  "commodity_origin_vs_assignment",
  "commodity_origin_mismatch_detail",
  "rf_corroboration",
  "rf_corroboration_detail",
  "subregion_withheld_explanation",
  "action_note",
  "report_statement",
  "vanbergen_macroregion",
  "vanbergen_macroregion_consistency"
)

front_cols <- unique(front_cols)
setcolorder(simplified, c(front_cols, setdiff(names(simplified), front_cols)))

# Explicitly remove is_reference if it somehow exists.
if ("is_reference" %in% names(simplified)) {
  simplified[, is_reference := NULL]
}

setorder(simplified, sample_id)
setorder(extended, sample_id)

# =============================================================================
# Optional validation and reduced-panel summaries
# =============================================================================

step4b_summary <- NULL
step4b_perclass <- NULL
step5_cv_summary <- NULL
step5_cv_perclass <- NULL
step6b_minK <- NULL
step6b_recommended <- NULL

if (isTRUE(include_validation_sheets)) {
  step4b_summary <- safe_read(file.path(opt$step4b_dir, "step4b_loo_summary_all_panels.tsv"))
  step4b_perclass <- safe_read(file.path(opt$step4b_dir, "step4b_loo_perclass_all_panels.tsv"))
  
  step5_cv_summary <- safe_read(file.path(opt$step5_dir, "RF_all_panels_cv_summary.tsv"))
  step5_cv_perclass <- safe_read(file.path(opt$step5_dir, "RF_all_panels_cv_perclass.tsv"))
  
  if (nzchar(opt$step6b_dir)) {
    step6b_minK <- safe_read(file.path(opt$step6b_dir, "Tables", "step6b_empirical_minK_by_target.tsv"))
    step6b_recommended <- safe_read(file.path(opt$step6b_dir, "Tables", "step6b_recommended_panels_index.tsv"))
  }
}

# =============================================================================
# Legend and column dictionary
# =============================================================================

legend_dt <- data.table(
  section = c("Primary decision",
    "Subregion reporting",
    "NOT_EVALUATED",
    "Threshold source",
    "Posterior category",
    "RF corroboration",
    "Commodity consistency",
    "Genomic affinity",
    "Missing observations",
    "Historical comparison"),
  explanation = c("PASS requires all six criteria at the evaluated resolution. FAIL requires at least one unmet criterion.",
    "Report subregion only when macroregion and its conditional subregion both PASS. If only macroregion passes, report macroregion.",
    "A failed macroregion prevents subregion evaluation. Criterion counts and subregion observations are unavailable, rather than zero.",
    "Thresholds are read from the Step 04 parameter table associated with the assignment file.",
    "High is at least the recorded high_posterior threshold. Moderate is at least moderate_posterior and below high_posterior. Both require six criteria.",
    "RF is a separate classifier comparison. It does not add a seventh criterion or alter the six criterion outcome.",
    "Compare genomic affinity with declared commodity country where available. Country of commodity dispatch does not establish biological source.",
    "Assignments are conditional on represented reference classes and do not demonstrate actual collection origin or a transport pathway.",
    "A missing metric in an evaluated branch fails its criterion. Missing required fields or inconsistent saved decisions stop report generation.",
    "Vanbergen et al. 2025 provides an additional macroregion comparison."))

column_key <- data.table(
  column = c("sample_id",
    "reported_origin",
    "reported_level",
    "reported_confidence",
    "macroregion_status",
    "subregion_status",
    "macroregion_criteria_met",
    "subregion_criteria_met",
    "macroregion_criteria_unmet",
    "subregion_criteria_unmet",
    "macroregion_unmet_criteria",
    "subregion_unmet_criteria",
    "criterion",
    "observed_field",
    "observed_value",
    "threshold",
    "result",
    "reason",
    "sample_origin_country",
    "sample_origin_country_source",
    "expected_macroregion_from_origin_country",
    "expected_subregion_from_origin_country",
    "commodity_origin_vs_assignment",
    "commodity_origin_mismatch_detail",
    "rf_corroboration",
    "rf_corroboration_detail",
    "subregion_withheld_explanation",
    "vanbergen_macroregion",
    "vanbergen_macroregion_consistency"),
  meaning = c("Unique query or intercept specimen identifier.",
    "Finest class that meets all six reporting criteria within an accepted hierarchy.",
    "Reported resolution: subregion, macroregion, or uncertain.",
    "Step 04 posterior support category for the accepted class. Uncertain when macroregion fails.",
    "PASS or FAIL for six macroregion criteria.",
    "PASS, FAIL, or NOT_EVALUATED for the conditional subregion.",
    "Number of the six macroregion criteria met.",
    "Number of subregion criteria met. Missing when the branch was not evaluated.",
    "Six minus macroregion_criteria_met.",
    "Six minus subregion_criteria_met. Missing when the branch was not evaluated.",
    "All unmet macroregion criterion identifiers. None when PASS.",
    "All unmet subregion criterion identifiers. None when PASS, missing when NOT_EVALUATED.",
    "Criterion identifier in Criteria. Every specimen and resolution has six rows.",
    "Lowercase Step 04 column name supplying the observation.",
    "Saved metric value. Unavailable for a branch that was not evaluated.",
    "Threshold recorded in the saved Step 04 parameter table.",
    "Criterion PASS, FAIL, or NOT_EVALUATED. Compare observed_value with threshold using operator.",
    "criterion_met, below_threshold, missing_value, or macroregion_failed.",
    "Declared commodity country where available. Unknown for trapped adults.",
    "Source of the declared commodity country.",
    "Macroregion assigned to the commodity country for the consistency comparison.",
    "Subregion assigned to the commodity country for the consistency comparison.",
    "Consistency comparison with commodity country, separate from criterion acceptance.",
    "Explanation of the consistency comparison.",
    "Separate RF support or caution relative to Step 04 candidate classes.",
    "Explanation of RF agreement and confidence. Does not alter criterion counts.",
    "All unmet subregion criteria, or the reason no branch was evaluated.",
    "Published macroregion interpretation from Vanbergen et al. 2025.",
    "Comparison of the Step 04 candidate with the published interpretation."))

# =============================================================================
# Write TSV/RDS outputs
# =============================================================================

simplified_file <- file.path(opt$out_dir, "FINAL_origin_tracing_authority_simplified.tsv")
extended_file <- file.path(opt$out_dir, "FINAL_origin_tracing_authority_extended.tsv")

fwrite(simplified, simplified_file, sep = "\t")
fwrite(extended, extended_file, sep = "\t")
fwrite(criteria_detail, file.path(opt$out_dir, "FINAL_origin_tracing_authority_criteria.tsv"), sep = "\t")
fwrite(reporting_thresholds, file.path(opt$out_dir, "reporting_thresholds.tsv"), sep = "\t")

saveRDS(simplified, file.path(opt$out_dir, "FINAL_origin_tracing_authority_simplified.rds"))
saveRDS(extended, file.path(opt$out_dir, "FINAL_origin_tracing_authority_extended.rds"))

# =============================================================================
# Write Excel workbook
# =============================================================================

out_xlsx <- file.path(opt$out_dir, opt$out_xlsx)

wb <- createWorkbook()

header_style <- createStyle(
  textDecoration = "bold",
  fgFill = "#D9EAF7",
  border = "Bottom",
  halign = "center",
  valign = "center"
)

amber_style <- createStyle(fgFill = "#FFEB9C")
red_style <- createStyle(fgFill = "#FFC7CE")
wrap_style <- createStyle(wrapText = TRUE, valign = "top")

add_table_sheet <- function(wb, sheet_name, dt, freeze_first_row = TRUE) {
  sn <- safe_sheet(wb, sheet_name)
  addWorksheet(wb, sn)
  
  writeDataTable(
    wb,
    sn,
    dt,
    tableStyle = "TableStyleMedium2",
    withFilter = TRUE
  )
  
  addStyle(
    wb,
    sn,
    header_style,
    rows = 1,
    cols = seq_len(ncol(dt)),
    gridExpand = TRUE
  )
  
  addStyle(
    wb,
    sn,
    wrap_style,
    rows = seq_len(nrow(dt) + 1L),
    cols = seq_len(ncol(dt)),
    gridExpand = TRUE,
    stack = TRUE
  )
  
  setColWidths(wb, sn, cols = seq_len(ncol(dt)), widths = "auto")
  # Bound text widths and fit wrapped rows for long criterion explanations.
  long_text <- grepl("explanation|detail|statement|basis|note|unmet_criteria|source|meaning", names(dt))
  if (any(long_text)) setColWidths(wb, sn, cols = which(long_text), widths = 48)
  setRowHeights(wb, sn, rows = 1, heights = 36)
  if (nrow(dt)) setRowHeights(wb, sn, rows = seq_len(nrow(dt)) + 1L, heights = 54)
  
  if (freeze_first_row) {
    freezePane(wb, sn, firstRow = TRUE)
  }
  
  sn
}

simp_sheet <- add_table_sheet(wb, "Simplified", simplified)

for (status_col in c("macroregion_status", "subregion_status")) {
  col <- which(names(simplified) == status_col)
  rows <- seq_len(nrow(simplified)) + 1L
  conditionalFormatting(wb, simp_sheet, cols = col, rows = rows,
                        rule = paste0(int2col(col), '2="FAIL"'), style = red_style)
}

if ("commodity_origin_vs_assignment" %in% names(simplified)) {
  mm_col <- which(names(simplified) == "commodity_origin_vs_assignment")
  rows <- seq_len(nrow(simplified)) + 1L
  
  conditionalFormatting(
    wb,
    simp_sheet,
    cols = mm_col,
    rows = rows,
    rule = 'NOT(ISERROR(SEARCH("Mismatch",INDIRECT(ADDRESS(ROW(),COLUMN())))))',
    style = amber_style
  )
}

add_table_sheet(wb, "Extended", extended)
add_table_sheet(wb, "Criteria", criteria_detail)
add_table_sheet(wb, "Reporting_Thresholds", reporting_thresholds)
add_table_sheet(wb, "Legend", legend_dt)
add_table_sheet(wb, "Column_Key", column_key)
add_table_sheet(wb, "Intercept_Origin_Map", intercept_origin)
add_table_sheet(wb, "Country_Region_Map", country_region_map)
add_table_sheet(wb, "Vanbergen_Reference", vanbergen_macroregion_table)

if (isTRUE(include_validation_sheets)) {
  if (!is.null(step4b_summary)) add_table_sheet(wb, "Step4b_Validation", step4b_summary)
  if (!is.null(step4b_perclass)) add_table_sheet(wb, "Step4b_PerClass", step4b_perclass)
  if (!is.null(step5_cv_summary)) add_table_sheet(wb, "Step5_RF_Validation", step5_cv_summary)
  if (!is.null(step5_cv_perclass)) add_table_sheet(wb, "Step5_RF_PerClass", step5_cv_perclass)
  if (!is.null(step6b_minK)) add_table_sheet(wb, "Step6b_MinK", step6b_minK)
  if (!is.null(step6b_recommended)) add_table_sheet(wb, "Step6b_Reduced_Panels", step6b_recommended)
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

# =============================================================================
# Run info
# =============================================================================

run_info <- list(
  step = "step7_reporting_authorities",
  repo_root = repo_root,
  step4_final = opt$step4_final,
  step4_params = opt$step4_params,
  reporting_parameters = reporting_parameters,
  step5_dir = opt$step5_dir,
  step4b_dir = opt$step4b_dir,
  step6b_dir = opt$step6b_dir,
  intercept_origin_file = opt$intercept_origin_file,
  out_dir = opt$out_dir,
  out_xlsx = out_xlsx,
  n_samples = nrow(simplified),
  criterion_counts = criteria_summary[, .N, by = .(resolution, status)],
  decision_consistency_with_step4 = TRUE,
  commodity_origin_counts = simplified[, .N, by = sample_origin_country],
  commodity_origin_vs_assignment_counts = simplified[, .N, by = commodity_origin_vs_assignment],
  vanbergen_consistency_counts = simplified[, .N, by = vanbergen_macroregion_consistency],
  subregion_withheld_explanation_blank_n = simplified[
    is.na(subregion_withheld_explanation) | subregion_withheld_explanation == "",
    .N
  ],
  include_validation_sheets = include_validation_sheets,
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(run_info, file.path(opt$out_dir, "step7_run_info.rds"))

# Keep the checkpoint reconstruction record as history after a successful native run.
reconstruction_record <- file.path(opt$out_dir, "reporting_reconstruction_info.json")
if (file.exists(reconstruction_record)) {
  history_record <- file.path(opt$out_dir, "reporting_reconstruction_history.json")
  if (file.exists(history_record)) {
    history_record <- file.path(opt$out_dir, paste0("reporting_reconstruction_history_",
      format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
  }
  if (!file.rename(reconstruction_record, history_record)) {
    warning("Could not move the checkpoint reconstruction record to history: ", reconstruction_record)
  }
}

message("\nDone.")
message("Simplified TSV: ", simplified_file)
message("Extended TSV:   ", extended_file)
message("Excel report:   ", out_xlsx)
message("Run info:        ", file.path(opt$out_dir, "step7_run_info.rds"))
