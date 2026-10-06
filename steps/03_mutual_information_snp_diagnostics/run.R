#!/usr/bin/env Rscript

# PESTFLY: Mutual information diagnostics
#
# Purpose
# Calculate reference only mutual information between the consensus allele and target class, in
# bits, as a diagnostic comparison with Hudson ranking.
#
# Interpretation
# Mutual information does not select, reorder or replace the Step 02 marker panels and does not
# enter the Step 04 likelihood. It is a descriptive comparison between two informativeness
# statistics.
#
# Technical notes
# The script removes unused class levels and includes MI self checks. MI output tables are
# written next to their Step 02 panels; the global index is written under
# results/03_mi_diagnostics.
#
# Run from the repository root:
#   Rscript steps/03_mutual_information_snp_diagnostics/run.R
# Detailed inputs, outputs and parameters are in the adjacent README.md.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# =============================================================================
# Repository-root detection
# =============================================================================

find_repo_root <- function(max_up = 10L) {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  cand <- wd
  
  for (i in seq_len(max_up + 1L)) {
    if (
      dir.exists(file.path(cand, "data")) &&
      dir.exists(file.path(cand, "steps"))
    ) {
      return(cand)
    }
    
    parent <- normalizePath(
      file.path(cand, ".."),
      winslash = "/",
      mustWork = FALSE
    )
    
    if (identical(parent, cand)) {
      break
    }
    
    cand <- parent
  }
  
  wd
}

repo_root <- find_repo_root()

is_abs_path <- function(p) {
  grepl("^(?:[A-Za-z]:[/\\\\]|/)", p)
}

resolve_path <- function(p, root = repo_root) {
  if (is.null(p) || !nzchar(p)) {
    return(p)
  }
  
  if (is_abs_path(p)) {
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  
  normalizePath(
    file.path(root, p),
    winslash = "/",
    mustWork = FALSE
  )
}

# =============================================================================
# Command-line options
# =============================================================================

option_list <- list(
  make_option(
    "--meta_path",
    type = "character",
    default = "results/01_qc/metadata_clean.tsv",
    help = "Path to metadata_clean.tsv [default %default]"
  ),
  make_option(
    "--step3_dir",
    type = "character",
    default = "results/02_snp_panels",
    help = paste(
      "Step 02 output directory containing panels_index.tsv",
      "[default %default]"
    )
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/03_mi_diagnostics",
    help = paste(
      "Step 03 output directory for MI index and run information",
      "[default %default]"
    )
  ),
  make_option(
    "--min_class_n",
    type = "integer",
    default = 3,
    help = paste(
      "Drop target classes with fewer than this many reference samples",
      "[default %default]"
    )
  ),
  make_option(
    "--min_snps_keep",
    type = "integer",
    default = 50,
    help = "Skip panels with fewer than this many SNPs [default %default]"
  ),
  make_option(
    "--include_p4_if_present",
    action = "store_true",
    default = FALSE,
    help = "Also compute MI for optional P4 if present [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$meta_path <- resolve_path(opt$meta_path)
opt$step3_dir <- resolve_path(opt$step3_dir)
opt$out_dir <- resolve_path(opt$out_dir)

dir.create(
  opt$out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

panels_index_path <- file.path(
  opt$step3_dir,
  "panels_index.tsv"
)

message("Repo root:    ", repo_root)
message("Metadata:     ", opt$meta_path)
message("Step 02 dir:  ", opt$step3_dir)
message("Output dir:   ", opt$out_dir)
message("Panels index: ", panels_index_path)
message("Min class n:  ", opt$min_class_n)

if (!file.exists(opt$meta_path)) {
  stop("Missing metadata file: ", opt$meta_path)
}

if (!file.exists(panels_index_path)) {
  stop("Missing panels_index.tsv: ", panels_index_path)
}

# =============================================================================
# Helper functions
# =============================================================================

clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("[\u200B-\u200D\uFEFF]", "", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x
}

as_clean_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  
  y <- tolower(clean_text(x))
  
  out <- rep(NA, length(y))
  
  out[
    y %in% c(
      "true",
      "t",
      "1",
      "yes",
      "y",
      "reference",
      "ref"
    )
  ] <- TRUE
  
  out[
    y %in% c(
      "false",
      "f",
      "0",
      "no",
      "n",
      "query",
      "intercept"
    )
  ] <- FALSE
  
  out
}

entropy_bits <- function(vec) {
  # Converting factors to character removes unused factor levels.
  # Retaining unused levels would introduce zero probabilities and make
  # 0 * log2(0) evaluate to NaN.
  
  vec <- as.character(vec)
  vec <- vec[!is.na(vec)]
  
  if (length(vec) == 0L) {
    return(NA_real_)
  }
  
  counts <- table(vec, useNA = "no")
  counts <- counts[counts > 0L]
  
  p <- as.numeric(counts) / sum(counts)
  
  -sum(p * log2(p))
}

mutual_information_bits <- function(x, y) {
  ok <- !(is.na(x) | is.na(y))
  
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  hy <- entropy_bits(y)
  hyx <- 0
  
  # Iterate only over observed genotype states.
  # This avoids zero-frequency conditional classes.
  
  for (lv in unique(x)) {
    idx <- x == lv
    p_lv <- sum(idx) / length(x)
    
    hyx <- hyx + p_lv * entropy_bits(y[idx])
  }
  
  mi <- hy - hyx
  
  # Mutual information cannot be negative.
  # Clamp only negligible floating-point error.
  
  if (
    is.finite(mi) &&
    mi < 0 &&
    mi > -1e-12
  ) {
    mi <- 0
  }
  
  mi
}

run_mi_self_test <- function() {
  # Unused target levels must not produce NaN entropy.
  
  y_unused <- factor(
    c("A", "A"),
    levels = c("A", "B")
  )
  
  stopifnot(
    isTRUE(
      all.equal(
        entropy_bits(y_unused),
        0,
        tolerance = 1e-12
      )
    )
  )
  
  # A perfectly associated balanced binary SNP must contain one bit
  # of information, even when the factor contains an unused level.
  
  y_perfect <- factor(
    c("A", "A", "B", "B"),
    levels = c("A", "B", "C")
  )
  
  mi_perfect <- mutual_information_bits(
    c(0L, 0L, 1L, 1L),
    y_perfect
  )
  
  stopifnot(
    isTRUE(
      all.equal(
        mi_perfect,
        1,
        tolerance = 1e-12
      )
    )
  )
  
  # A balanced independent SNP and target must contain zero information.
  
  y_independent <- factor(
    c("A", "A", "B", "B"),
    levels = c("A", "B", "C")
  )
  
  mi_independent <- mutual_information_bits(
    c(0L, 1L, 0L, 1L),
    y_independent
  )
  
  stopifnot(
    isTRUE(
      all.equal(
        mi_independent,
        0,
        tolerance = 1e-12
      )
    )
  )
  
  message("MI helper self-test passed.")
  
  invisible(TRUE)
}

run_mi_self_test()

panel_to_groupcol <- function(panel_id) {
  if (grepl("^P1_", panel_id)) {
    return("macroregion_3")
  }
  
  if (grepl("^P2_", panel_id)) {
    return("subregion")
  }
  
  if (grepl("^P3_", panel_id)) {
    return("subregion")
  }
  
  if (grepl("^P4_", panel_id)) {
    return("population")
  }
  
  if (grepl("macroregion", panel_id, ignore.case = TRUE)) {
    return("macroregion_3")
  }
  
  if (grepl("subregion", panel_id, ignore.case = TRUE)) {
    return("subregion")
  }
  
  if (grepl("population", panel_id, ignore.case = TRUE)) {
    return("population")
  }
  
  "macroregion_3"
}

fix_panel_dir <- function(panel_dir, step3_dir) {
  if (is.na(panel_dir) || !nzchar(panel_dir)) {
    return(NA_character_)
  }
  
  if (dir.exists(panel_dir)) {
    return(
      normalizePath(
        panel_dir,
        winslash = "/",
        mustWork = FALSE
      )
    )
  }
  
  p2 <- file.path(
    step3_dir,
    panel_dir
  )
  
  if (dir.exists(p2)) {
    return(
      normalizePath(
        p2,
        winslash = "/",
        mustWork = FALSE
      )
    )
  }
  
  p3 <- file.path(
    step3_dir,
    "panels",
    basename(panel_dir)
  )
  
  if (dir.exists(p3)) {
    return(
      normalizePath(
        p3,
        winslash = "/",
        mustWork = FALSE
      )
    )
  }
  
  normalizePath(
    panel_dir,
    winslash = "/",
    mustWork = FALSE
  )
}

panel_reference_filter <- function(meta_ref, panel_id) {
  if (grepl("^P1_", panel_id)) {
    return(
      meta_ref[
        macroregion_3 %in% c("Africa", "Asia")
      ]
    )
  }
  
  if (grepl("^P2_", panel_id)) {
    return(
      meta_ref[
        macroregion_3 == "Africa"
      ]
    )
  }
  
  if (grepl("^P3_", panel_id)) {
    return(
      meta_ref[
        macroregion_3 == "Asia"
      ]
    )
  }
  
  meta_ref
}

# =============================================================================
# Load metadata
# =============================================================================

meta <- fread(opt$meta_path)

setnames(
  meta,
  tolower(names(meta))
)

required_meta <- c(
  "sample_id",
  "is_reference",
  "macroregion_3",
  "subregion",
  "population"
)

missing_meta <- setdiff(
  required_meta,
  names(meta)
)

if (length(missing_meta) > 0L) {
  stop(
    "metadata_clean.tsv is missing required column(s): ",
    paste(missing_meta, collapse = ", ")
  )
}

for (nm in names(meta)) {
  if (
    !inherits(
      meta[[nm]],
      c(
        "numeric",
        "integer",
        "logical",
        "Date",
        "POSIXct"
      )
    )
  ) {
    meta[[nm]] <- clean_text(meta[[nm]])
  }
}

meta[, sample_id := clean_text(sample_id)]
meta[, is_reference := as_clean_logical(is_reference)]
meta[, macroregion_3 := clean_text(macroregion_3)]
meta[, subregion := clean_text(subregion)]
meta[, population := clean_text(population)]

# =============================================================================
# Load Step 02 panel index
# =============================================================================

panels_index <- fread(panels_index_path)

setnames(
  panels_index,
  tolower(names(panels_index))
)

required_pi <- c(
  "panel_id",
  "panel_dir"
)

missing_pi <- setdiff(
  required_pi,
  names(panels_index)
)

if (length(missing_pi) > 0L) {
  stop(
    "panels_index.tsv is missing required column(s): ",
    paste(missing_pi, collapse = ", ")
  )
}

if (!isTRUE(opt$include_p4_if_present)) {
  panels_index <- panels_index[
    !grepl("^P4_", panel_id)
  ]
}

if (nrow(panels_index) == 0L) {
  stop("No Step 02 panels are available for MI calculation.")
}

# =============================================================================
# Process panels
# =============================================================================

mi_index <- list()

for (i in seq_len(nrow(panels_index))) {
  panel_id <- as.character(
    panels_index$panel_id[i]
  )
  
  panel_dir <- fix_panel_dir(
    as.character(panels_index$panel_dir[i]),
    opt$step3_dir
  )
  
  message("")
  message("============================================================")
  message("Panel: ", panel_id)
  message("Dir:   ", panel_dir)
  message("============================================================")
  
  snp_matrix_path <- file.path(
    panel_dir,
    "snp_matrix.rds"
  )
  
  snp_map_path <- file.path(
    panel_dir,
    "snp_map.tsv"
  )
  
  if (!file.exists(snp_matrix_path)) {
    message("Skipping panel: missing snp_matrix.rds")
    next
  }
  
  if (!file.exists(snp_map_path)) {
    message("Skipping panel: missing snp_map.tsv")
    next
  }
  
  X <- readRDS(snp_matrix_path)
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  storage.mode(X) <- "integer"
  
  if (
    is.null(rownames(X)) ||
    is.null(colnames(X))
  ) {
    message(
      "Skipping panel: the SNP matrix is missing row names or column names"
    )
    
    next
  }
  
  if (nrow(X) < opt$min_snps_keep) {
    message(
      "Skipping panel: too few SNPs: ",
      nrow(X)
    )
    
    next
  }
  
  snp_map <- fread(snp_map_path)
  
  if (!("snp_id" %in% names(snp_map))) {
    message(
      "Skipping panel: snp_map.tsv lacks the snp_id column"
    )
    
    next
  }
  
  group_col <- panel_to_groupcol(panel_id)
  
  if (!(group_col %in% names(meta))) {
    message(
      "Skipping panel: metadata lacks group column ",
      group_col
    )
    
    next
  }
  
  ids_in_matrix <- colnames(X)
  
  meta_ref <- meta[
    is_reference %in% TRUE &
      sample_id %in% ids_in_matrix
  ]
  
  meta_ref <- meta_ref[
    !is.na(get(group_col)) &
      clean_text(get(group_col)) != ""
  ]
  
  meta_ref <- panel_reference_filter(
    meta_ref,
    panel_id
  )
  
  if (nrow(meta_ref) < 2L) {
    message(
      "Skipping panel: too few reference samples after panel filtering"
    )
    
    next
  }
  
  class_counts <- meta_ref[
    ,
    .N,
    by = group_col
  ]
  
  setnames(
    class_counts,
    group_col,
    "class"
  )
  
  keep_classes <- class_counts[
    N >= opt$min_class_n,
    class
  ]
  
  meta_ref <- meta_ref[
    get(group_col) %in% keep_classes
  ]
  
  if (
    nrow(meta_ref) < 2L ||
    uniqueN(meta_ref[[group_col]]) < 2L
  ) {
    message(
      "Skipping panel: fewer than two classes remain after ",
      "the min_class_n filter"
    )
    
    next
  }
  
  ref_ids <- meta_ref$sample_id
  
  y <- as.factor(
    as.character(meta_ref[[group_col]])
  )
  
  X_ref <- X[
    ,
    ref_ids,
    drop = FALSE
  ]
  
  res <- data.table(
    snp_id = rownames(X_ref),
    mi_bits = NA_real_,
    n_nonmissing = as.integer(
      rowSums(!is.na(X_ref))
    ),
    n_ref_samples = length(ref_ids),
    n_classes = nlevels(y),
    target = group_col
  )
  
  for (s in seq_len(nrow(X_ref))) {
    x <- X_ref[s, ]
    
    if (length(unique(x[!is.na(x)])) < 2L) {
      res$mi_bits[s] <- 0
    } else {
      res$mi_bits[s] <- mutual_information_bits(
        x,
        y
      )
    }
  }
  
  # Stop explicitly if the calculation produces invalid MI values.
  
  bad_nonfinite <- which(
    !is.finite(res$mi_bits)
  )
  
  if (length(bad_nonfinite) > 0L) {
    stop(
      "Panel ",
      panel_id,
      " produced ",
      length(bad_nonfinite),
      " non-finite MI values."
    )
  }
  
  bad_negative <- which(
    res$mi_bits < -1e-12
  )
  
  if (length(bad_negative) > 0L) {
    stop(
      "Panel ",
      panel_id,
      " produced ",
      length(bad_negative),
      " negative MI values."
    )
  }
  
  panel_entropy <- entropy_bits(y)
  
  bad_upper <- which(
    res$mi_bits > panel_entropy + 1e-12
  )
  
  if (length(bad_upper) > 0L) {
    stop(
      "Panel ",
      panel_id,
      " produced ",
      length(bad_upper),
      " MI values greater than the target entropy."
    )
  }
  
  setorder(
    res,
    -mi_bits,
    -n_nonmissing
  )
  
  res <- merge(
    res,
    snp_map,
    by = "snp_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  setorder(
    res,
    -mi_bits,
    -n_nonmissing
  )
  
  out_file <- file.path(
    panel_dir,
    "snp_mi.tsv"
  )
  
  fwrite(
    res,
    out_file,
    sep = "\t"
  )
  
  saveRDS(
    res,
    file.path(
      panel_dir,
      "snp_mi.rds"
    )
  )
  
  message("Wrote: ", out_file)
  
  mi_index[[length(mi_index) + 1L]] <- data.table(
    panel_id = panel_id,
    target = group_col,
    panel_dir = panel_dir,
    snp_matrix_rds = snp_matrix_path,
    mi_file = out_file,
    n_snps = nrow(res),
    n_ref_samples = length(ref_ids),
    n_classes = nlevels(y),
    classes = paste(
      levels(y),
      collapse = ","
    )
  )
}

# =============================================================================
# Write global Step 03 outputs
# =============================================================================

mi_index_dt <- rbindlist(
  mi_index,
  fill = TRUE
)

out_index <- file.path(
  opt$out_dir,
  "panels_mi_index.tsv"
)

fwrite(
  mi_index_dt,
  out_index,
  sep = "\t"
)

saveRDS(
  mi_index_dt,
  file.path(
    opt$out_dir,
    "panels_mi_index.rds"
  )
)

run_info <- list(
  step = "step03_mutual_information_snp_diagnostics",
  mi_implementation = paste(
    "discrete entropy with unused factor levels removed"
  ),
  mi_self_test_passed = TRUE,
  mi_validation = paste(
    "finite, non-negative, and no greater than target entropy",
    "for every reported locus"
  ),
  repo_root = repo_root,
  metadata = opt$meta_path,
  step3_dir = opt$step3_dir,
  out_dir = opt$out_dir,
  min_class_n = opt$min_class_n,
  min_snps_keep = opt$min_snps_keep,
  include_p4_if_present = opt$include_p4_if_present,
  mi_index = mi_index_dt,
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(
  run_info,
  file.path(
    opt$out_dir,
    "step03_run_info.rds"
  )
)

message("")
message("Done.")
message("MI index: ", out_index)
message("Rows in MI index: ", nrow(mi_index_dt))
