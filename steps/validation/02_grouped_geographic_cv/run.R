#!/usr/bin/env Rscript

# PESTFLY supplementary analysis: Grouped geographic validation
#
# Purpose
# Withhold all reference specimens from each country or country nested collection site, then
# evaluate predictions for the held out geographic unit.
#
# Interpretation
# Each evaluable fold filters the retained SNP set, recalculates corrected Hudson scores,
# reranks markers and estimates allele frequencies using training references only. The retained
# one SNP per ortholog candidates were discovered beforehand. Alternative sites are not
# rediscovered from alignments, so this is not fully nested de novo feature discovery.
#
# Technical notes
# A fold is not evaluable if the true class is absent after withholding, or training
# representation fails required support. The benchmark Central Africa country and site folds and
# East Asia country fold illustrate missing class coverage. Report evaluable sample counts
# alongside accuracy.
#
# Run from the repository root:
#   Rscript steps/validation/02_grouped_geographic_cv/run.R
# See the adjacent README.md for inputs, outputs and complete CLI defaults.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(openxlsx)
  library(parallel)
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

# =============================================================================
# General helpers
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
  if (is.logical(x)) return(x)

  y <- tolower(clean_text(x))
  out <- rep(NA, length(y))

  out[y %in% c("true", "t", "1", "yes", "y", "reference", "ref")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n", "query", "intercept")] <- FALSE

  out
}

parse_bool1 <- function(x, default = FALSE, argname = "") {
  if (is.null(x) || length(x) == 0L) return(as.logical(default)[1])

  y <- tolower(clean_text(x[1]))

  if (is.na(y) || !nzchar(y)) return(as.logical(default)[1])
  if (y %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (y %in% c("false", "f", "0", "no", "n")) return(FALSE)

  stop("Cannot parse boolean for ", argname, ": '", x, "'. Use TRUE or FALSE.")
}

parse_integer_vector <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) {
    return(integer())
  }

  out <- suppressWarnings(as.integer(clean_text(strsplit(x, ",")[[1]])))
  out <- out[!is.na(out) & out > 0L]
  sort(unique(out))
}

parse_character_vector <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) {
    return(character())
  }

  out <- clean_text(strsplit(x, ",")[[1]])
  sort(unique(out[!is.na(out) & nzchar(out)]))
}

collapse_values <- function(x) {
  x <- sort(unique(clean_text(x)))
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0L) return(NA_character_)
  paste(x, collapse = ";")
}

safe_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  median(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

conf_label <- function(pmax, high = 0.95, mod = 0.85) {
  if (is.na(pmax)) return("NA")
  if (pmax >= high) return("High")
  if (pmax >= mod) return("Moderate")
  "Low"
}

# =============================================================================
# Command line options
# =============================================================================

option_list <- list(
  make_option(
    "--qc_dir",
    type = "character",
    default = "results/01_qc",
    help = "Directory containing metadata_clean.tsv [default %default]"
  ),
  make_option(
    "--step2_dir",
    type = "character",
    default = "results/02_snp_panels",
    help = "Corrected Step 02 panel directory [default %default]"
  ),
  make_option(
    "--loo_dir",
    type = "character",
    default = "results/05_loo_validation",
    help = "Existing individual LOO directory used for comparison [default %default]"
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/validation/02_grouped_geographic_cv",
    help = "Separate grouped validation output directory [default %default]"
  ),
  make_option(
    "--panels",
    type = "character",
    default = paste(
      c(
        "P1_macroregion_africa_vs_asia",
        "P2_subregion_within_africa",
        "P3_subregion_within_asia"
      ),
      collapse = ","
    ),
    help = "Comma separated panel identifiers [default %default]"
  ),
  make_option(
    "--fold_units",
    type = "character",
    default = "country,site",
    help = "Grouped validation units: country, site, or both [default %default]"
  ),
  make_option(
    "--n_cores",
    type = "integer",
    default = 0L,
    help = "Parallel workers. Zero uses up to eight detected cores minus one [default %default]"
  ),
  make_option(
    "--min_ref_per_group",
    type = "integer",
    default = 3L,
    help = "Minimum training references required for a candidate class [default %default]"
  ),
  make_option(
    "--max_site_missing",
    type = "double",
    default = 0.20,
    help = "Maximum training missing fraction for a retained SNP [default %default]"
  ),
  make_option(
    "--min_mac",
    type = "integer",
    default = 2L,
    help = "Minimum training minor allele count for a retained SNP [default %default]"
  ),
  make_option(
    "--foldwise_reranking",
    type = "character",
    default = "TRUE",
    help = "Recalculate Hudson FST and rerank retained SNPs per fold [default %default]"
  ),
  make_option(
    "--pseudocount",
    type = "double",
    default = 0.5,
    help = "Allele frequency pseudocount [default %default]"
  ),
  make_option(
    "--epsilon",
    type = "double",
    default = 0.02,
    help = "Per site allele flip probability [default %default]"
  ),
  make_option(
    "--min_snps_query_macroregion",
    type = "integer",
    default = 200L,
    help = "Minimum usable SNPs for P1 [default %default]"
  ),
  make_option(
    "--min_snps_query_subregion",
    type = "integer",
    default = 300L,
    help = "Minimum usable SNPs for P2 and P3 [default %default]"
  ),
  make_option(
    "--high_posterior",
    type = "double",
    default = 0.95,
    help = "High posterior threshold [default %default]"
  ),
  make_option(
    "--moderate_posterior",
    type = "double",
    default = 0.85,
    help = "Moderate posterior threshold [default %default]"
  ),
  make_option(
    "--min_gap_macroregion",
    type = "double",
    default = 0.20,
    help = "Minimum P1 posterior gap [default %default]"
  ),
  make_option(
    "--min_gap_subregion",
    type = "double",
    default = 0.25,
    help = "Minimum P2 and P3 posterior gap [default %default]"
  ),
  make_option(
    "--base_K",
    type = "character",
    default = "1,2,3,4,5,10,20,50,100,200,500",
    help = "Comma separated base K values [default %default]"
  ),
  make_option(
    "--extra_K",
    type = "character",
    default = "1000,2000,5000,10000",
    help = "Additional K values when available [default %default]"
  ),
  make_option(
    "--auto_extend_K",
    type = "character",
    default = "TRUE",
    help = "Use extra K values when available [default %default]"
  ),
  make_option(
    "--include_all_snps_K",
    type = "character",
    default = "TRUE",
    help = "Include all eligible SNPs as the final K [default %default]"
  ),
  make_option(
    "--max_K_points",
    type = "integer",
    default = 30L,
    help = "Maximum number of K values [default %default]"
  ),
  make_option(
    "--min_agreement",
    type = "double",
    default = 0.90,
    help = "Minimum agreement across usable K values [default %default]"
  ),
  make_option(
    "--min_K_available",
    type = "integer",
    default = 6L,
    help = "Minimum number of usable K values [default %default]"
  ),
  make_option(
    "--tail_fraction",
    type = "double",
    default = 0.50,
    help = "Largest K fraction used for tail stability [default %default]"
  ),
  make_option(
    "--min_tail_agreement",
    type = "double",
    default = 1.00,
    help = "Required agreement in the largest K tail [default %default]"
  ),
  make_option(
    "--top_markers_per_fold",
    type = "integer",
    default = 20L,
    help = "Number of fold ranked markers retained for audit [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$qc_dir <- resolve_path(opt$qc_dir)
opt$step2_dir <- resolve_path(opt$step2_dir)
opt$loo_dir <- resolve_path(opt$loo_dir)
opt$out_dir <- resolve_path(opt$out_dir)

panels_requested <- parse_character_vector(opt$panels)
fold_units <- parse_character_vector(opt$fold_units)

allowed_fold_units <- c("country", "site")
bad_fold_units <- setdiff(fold_units, allowed_fold_units)

if (length(fold_units) == 0L || length(bad_fold_units) > 0L) {
  stop(
    "--fold_units must contain country, site, or both. Invalid value(s): ",
    paste(bad_fold_units, collapse = ", ")
  )
}

foldwise_reranking <- parse_bool1(
  opt$foldwise_reranking,
  default = TRUE,
  argname = "--foldwise_reranking"
)

auto_extend_K <- parse_bool1(
  opt$auto_extend_K,
  default = TRUE,
  argname = "--auto_extend_K"
)

include_all_snps_K <- parse_bool1(
  opt$include_all_snps_K,
  default = TRUE,
  argname = "--include_all_snps_K"
)

baseK <- parse_integer_vector(opt$base_K)
extraK <- parse_integer_vector(opt$extra_K)

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))

if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}

if (is.na(opt$n_cores) || opt$n_cores <= 0L) {
  opt$n_cores <- max(1L, min(8L, detected_cores - 1L))
}

opt$n_cores <- max(1L, as.integer(opt$n_cores))

settings <- list(
  min_ref_per_group = as.integer(opt$min_ref_per_group),
  max_site_missing = as.numeric(opt$max_site_missing),
  min_mac = as.integer(opt$min_mac),
  foldwise_reranking = foldwise_reranking,
  pseudocount = as.numeric(opt$pseudocount),
  epsilon = as.numeric(opt$epsilon),
  high_posterior = as.numeric(opt$high_posterior),
  moderate_posterior = as.numeric(opt$moderate_posterior),
  min_agreement = as.numeric(opt$min_agreement),
  min_K_available = as.integer(opt$min_K_available),
  tail_fraction = as.numeric(opt$tail_fraction),
  min_tail_agreement = as.numeric(opt$min_tail_agreement),
  auto_extend_K = auto_extend_K,
  include_all_snps_K = include_all_snps_K,
  max_K_points = as.integer(opt$max_K_points),
  top_markers_per_fold = max(0L, as.integer(opt$top_markers_per_fold))
)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Repository root:    ", repo_root)
message("QC directory:       ", opt$qc_dir)
message("Corrected panels:   ", opt$step2_dir)
message("Individual LOO:     ", opt$loo_dir)
message("Grouped CV output:  ", opt$out_dir)
message("Panels:             ", paste(panels_requested, collapse = ", "))
message("Fold units:         ", paste(fold_units, collapse = ", "))
message("Parallel workers:   ", opt$n_cores)
message("Foldwise reranking: ", foldwise_reranking)

# =============================================================================
# Panel configuration
# =============================================================================

panel_config <- function(panel_id) {
  if (panel_id == "P1_macroregion_africa_vs_asia") {
    return(list(
      panel_id = panel_id,
      panel_label = "P1 macroregion",
      panel_type = "binary",
      group_col = "macroregion_3",
      macroregion_filter = c("Africa", "Asia"),
      min_snps_query = as.integer(opt$min_snps_query_macroregion),
      min_gap = as.numeric(opt$min_gap_macroregion)
    ))
  }

  if (panel_id == "P2_subregion_within_africa") {
    return(list(
      panel_id = panel_id,
      panel_label = "P2 Africa subregion",
      panel_type = "multiclass",
      group_col = "subregion",
      macroregion_filter = "Africa",
      min_snps_query = as.integer(opt$min_snps_query_subregion),
      min_gap = as.numeric(opt$min_gap_subregion)
    ))
  }

  if (panel_id == "P3_subregion_within_asia") {
    return(list(
      panel_id = panel_id,
      panel_label = "P3 Asia subregion",
      panel_type = "multiclass",
      group_col = "subregion",
      macroregion_filter = "Asia",
      min_snps_query = as.integer(opt$min_snps_query_subregion),
      min_gap = as.numeric(opt$min_gap_subregion)
    ))
  }

  stop("Unknown panel identifier: ", panel_id)
}

# =============================================================================
# Standard Hudson FST scoring
# =============================================================================

hudson_fst_vector <- function(p1, p2, n1, n2) {
  out <- rep(NA_real_, length(p1))

  ok <- is.finite(p1) & is.finite(p2) &
    is.finite(n1) & is.finite(n2) &
    n1 >= 2 & n2 >= 2

  if (!any(ok)) return(out)

  numerator <- (p1 - p2)^2 -
    p1 * (1 - p1) / (n1 - 1) -
    p2 * (1 - p2) / (n2 - 1)

  denominator <- p1 * (1 - p2) + p2 * (1 - p1)

  ok <- ok & is.finite(numerator) & is.finite(denominator) & denominator > 0
  out[ok] <- numerator[ok] / denominator[ok]

  # Negative finite estimates are intentionally retained.
  out
}

score_fold_snps <- function(
    X,
    train_meta,
    group_col,
    panel_type,
    snp_map,
    max_site_missing,
    min_mac,
    foldwise_reranking
) {
  train_ids <- intersect(train_meta$sample_id, colnames(X))

  if (length(train_ids) < 2L) return(data.table())

  X_train <- X[, train_ids, drop = FALSE]
  n_nonmissing <- rowSums(!is.na(X_train))
  missing_frac <- 1 - n_nonmissing / length(train_ids)

  alt_count <- rowSums(X_train == 1L, na.rm = TRUE)
  ref_count <- rowSums(X_train == 0L, na.rm = TRUE)
  mac <- pmin(alt_count, ref_count)

  eligible_basic <- is.finite(missing_frac) &
    missing_frac <= max_site_missing &
    mac >= min_mac

  groups <- sort(unique(clean_text(train_meta[[group_col]])))
  groups <- groups[!is.na(groups) & nzchar(groups)]

  if (length(groups) < 2L) return(data.table())

  pmat <- matrix(
    NA_real_,
    nrow = nrow(X),
    ncol = length(groups),
    dimnames = list(rownames(X), groups)
  )

  nmat <- matrix(
    0L,
    nrow = nrow(X),
    ncol = length(groups),
    dimnames = list(rownames(X), groups)
  )

  for (g in groups) {
    ids_g <- intersect(
      train_meta[get(group_col) == g, sample_id],
      colnames(X)
    )

    if (length(ids_g) == 0L) next

    Xg <- X[, ids_g, drop = FALSE]
    ng <- rowSums(!is.na(Xg))
    ag <- rowSums(Xg == 1L, na.rm = TRUE)

    pg <- rep(NA_real_, length(ng))
    has_data <- ng > 0L
    pg[has_data] <- ag[has_data] / ng[has_data]

    pmat[, g] <- pg
    nmat[, g] <- ng
  }

  fold_score <- rep(NA_real_, nrow(X))
  best_pair <- rep(NA_character_, nrow(X))

  if (isTRUE(foldwise_reranking)) {
    if (panel_type == "binary") {
      g1 <- groups[1]
      g2 <- groups[2]

      fold_score <- hudson_fst_vector(
        pmat[, g1],
        pmat[, g2],
        nmat[, g1],
        nmat[, g2]
      )

      best_pair[is.finite(fold_score)] <- paste(g1, g2, sep = "__vs__")

    } else {
      for (a in seq_len(length(groups) - 1L)) {
        for (b in (a + 1L):length(groups)) {
          g1 <- groups[a]
          g2 <- groups[b]

          score_ab <- hudson_fst_vector(
            pmat[, g1],
            pmat[, g2],
            nmat[, g1],
            nmat[, g2]
          )

          take <- is.finite(score_ab) &
            (!is.finite(fold_score) | score_ab > fold_score)

          fold_score[take] <- score_ab[take]
          best_pair[take] <- paste(g1, g2, sep = "__vs__")
        }
      }
    }

  } else {
    global_scores <- snp_map$score[match(rownames(X), snp_map$snp_id)]
    fold_score <- as.numeric(global_scores)
    best_pair <- snp_map$best_pair[match(rownames(X), snp_map$snp_id)]
  }

  keep <- eligible_basic & is.finite(fold_score)

  if (!any(keep)) return(data.table())

  out <- data.table(
    snp_id = rownames(X)[keep],
    fold_score = fold_score[keep],
    fold_best_pair = best_pair[keep],
    fold_missing_frac = missing_frac[keep],
    fold_mac = as.integer(mac[keep]),
    fold_n_nonmissing = as.integer(n_nonmissing[keep])
  )

  idx_map <- match(out$snp_id, snp_map$snp_id)

  out[, `:=`(
    og = if ("og" %in% names(snp_map)) snp_map$og[idx_map] else NA_character_,
    original_score = if ("score" %in% names(snp_map)) as.numeric(snp_map$score[idx_map]) else NA_real_,
    original_global_rank = if ("global_rank" %in% names(snp_map)) as.integer(snp_map$global_rank[idx_map]) else NA_integer_
  )]

  setorder(
    out,
    -fold_score,
    fold_missing_frac,
    -fold_mac,
    snp_id
  )

  out[, fold_rank := seq_len(.N)]
  out
}

# =============================================================================
# Likelihood and multi K helpers
# =============================================================================

loglik_one_group <- function(qvec, pvec, epsilon = 0.02) {
  ok <- !is.na(qvec) & !is.na(pvec)

  if (!any(ok)) return(NA_real_)

  q <- qvec[ok]
  p <- pvec[ok]

  p1 <- p * (1 - epsilon) + (1 - p) * epsilon
  p0 <- (1 - p) * (1 - epsilon) + p * epsilon

  p1 <- pmin(pmax(p1, 1e-12), 1 - 1e-12)
  p0 <- pmin(pmax(p0, 1e-12), 1 - 1e-12)

  sum(ifelse(q == 1L, log(p1), log(p0)))
}

post_from_ll <- function(lls) {
  if (all(is.na(lls))) return(rep(NA_real_, length(lls)))

  m <- max(lls, na.rm = TRUE)
  w <- exp(lls - m)
  w[is.na(w)] <- 0

  if (sum(w) <= 0) return(rep(NA_real_, length(lls)))

  out <- w / sum(w)
  names(out) <- names(lls)
  out
}

make_Ks <- function(
    n_snps,
    baseK,
    extraK,
    auto_extend,
    include_all,
    max_points
) {
  if (is.na(n_snps) || n_snps <= 0L) return(integer())

  baseK <- sort(unique(as.integer(baseK)))
  baseK <- baseK[!is.na(baseK) & baseK > 0L]

  extraK <- sort(unique(as.integer(extraK)))
  extraK <- extraK[!is.na(extraK) & extraK > 0L]

  Ks <- baseK[baseK <= n_snps]

  if (isTRUE(auto_extend)) {
    Ks <- sort(unique(c(Ks, extraK[extraK <= n_snps])))
  }

  if (isTRUE(include_all)) {
    Ks <- sort(unique(c(Ks, as.integer(n_snps))))
  }

  Ks <- sort(unique(Ks[Ks <= n_snps]))

  if (length(Ks) > max_points) {
    small_K <- Ks[Ks <= 500L]
    large_K <- Ks[Ks > 500L]

    if (length(large_K) > 0L) {
      keep_large <- unique(round(seq(
        from = min(large_K),
        to = max(large_K),
        length.out = max(1L, max_points - length(small_K))
      )))

      Ks <- sort(unique(c(small_K, keep_large, max(Ks))))
    }

    if (length(Ks) > max_points) {
      Ks <- sort(unique(c(head(Ks, max_points - 1L), max(Ks))))
    }
  }

  Ks
}

build_ref_freqs <- function(snp_mat, meta_ref, group_col, pseudocount) {
  ref_ids <- intersect(meta_ref$sample_id, colnames(snp_mat))

  if (length(ref_ids) < 2L) {
    return(list(freq = NULL, groups = character(), ref_ids = ref_ids))
  }

  meta_ref2 <- meta_ref[sample_id %in% ref_ids]
  meta_ref2 <- meta_ref2[
    !is.na(get(group_col)) & clean_text(get(group_col)) != ""
  ]

  groups <- sort(unique(meta_ref2[[group_col]]))
  groups <- groups[!is.na(groups) & nzchar(groups)]

  if (length(groups) < 2L) {
    return(list(freq = NULL, groups = groups, ref_ids = ref_ids))
  }

  freq <- matrix(
    NA_real_,
    nrow = nrow(snp_mat),
    ncol = length(groups),
    dimnames = list(rownames(snp_mat), groups)
  )

  for (g in groups) {
    ids <- intersect(
      meta_ref2[get(group_col) == g, sample_id],
      colnames(snp_mat)
    )

    if (length(ids) == 0L) next

    Xg <- snp_mat[, ids, drop = FALSE]
    alt_count <- rowSums(Xg == 1L, na.rm = TRUE)
    ref_count <- rowSums(Xg == 0L, na.rm = TRUE)
    total_count <- alt_count + ref_count

    p <- rep(NA_real_, length(total_count))
    ok <- total_count > 0L

    p[ok] <- (alt_count[ok] + pseudocount) /
      (total_count[ok] + 2 * pseudocount)

    freq[, g] <- p
  }

  list(freq = freq, groups = groups, ref_ids = ref_ids)
}

assign_one_K <- function(
    snp_mat,
    freq,
    groups,
    sample_id,
    panel_id,
    group_col,
    K,
    epsilon,
    min_snps_query,
    high_thr,
    mod_thr,
    min_gap
) {
  qvec <- snp_mat[, sample_id]
  any_frequency <- rowSums(!is.na(freq)) > 0L
  n_used <- sum(!is.na(qvec) & any_frequency)

  if (is.na(n_used) || n_used < min_snps_query) {
    return(data.table(
      panel_id = panel_id,
      group_col = group_col,
      K = as.integer(K),
      sample_id = sample_id,
      top_group = NA_character_,
      top_posterior = NA_real_,
      second_group = NA_character_,
      second_posterior = NA_real_,
      posterior_gap = NA_real_,
      n_snps_used = as.integer(n_used),
      confidence = "Low",
      status = "too_few_snps",
      note = paste0("Need at least ", min_snps_query, "; have ", n_used)
    ))
  }

  lls <- vapply(
    groups,
    function(g) loglik_one_group(qvec, freq[, g], epsilon = epsilon),
    numeric(1)
  )

  post <- post_from_ll(lls)

  if (all(is.na(post))) {
    return(data.table(
      panel_id = panel_id,
      group_col = group_col,
      K = as.integer(K),
      sample_id = sample_id,
      top_group = NA_character_,
      top_posterior = NA_real_,
      second_group = NA_character_,
      second_posterior = NA_real_,
      posterior_gap = NA_real_,
      n_snps_used = as.integer(n_used),
      confidence = "Low",
      status = "no_likelihood",
      note = "All likelihoods are missing"
    ))
  }

  ord <- order(post, decreasing = TRUE)
  top1 <- ord[1]
  top2 <- if (length(ord) >= 2L) ord[2] else NA_integer_

  top_group <- names(post)[top1]
  top_posterior <- unname(post[top1])
  second_group <- if (!is.na(top2)) names(post)[top2] else NA_character_
  second_posterior <- if (!is.na(top2)) unname(post[top2]) else NA_real_
  gap <- if (!is.na(second_posterior)) {
    top_posterior - second_posterior
  } else {
    NA_real_
  }

  confidence <- conf_label(
    top_posterior,
    high = high_thr,
    mod = mod_thr
  )

  status <- "ok"
  note <- ""

  if (!is.na(gap) && gap < min_gap) {
    status <- "low_gap"
    note <- paste0("posterior_gap<", min_gap)
  }

  data.table(
    panel_id = panel_id,
    group_col = group_col,
    K = as.integer(K),
    sample_id = sample_id,
    top_group = top_group,
    top_posterior = top_posterior,
    second_group = second_group,
    second_posterior = second_posterior,
    posterior_gap = gap,
    n_snps_used = as.integer(n_used),
    confidence = confidence,
    status = status,
    note = note
  )
}

stable_from_threshold <- function(dt, final_group) {
  dt2 <- dt[status %in% c("ok", "low_gap") & !is.na(top_group)]

  if (nrow(dt2) == 0L || is.na(final_group)) return(NA_real_)

  setorder(dt2, K)
  Ks_sorted <- sort(unique(dt2$K))

  for (k in Ks_sorted) {
    later <- dt2[K >= k]

    if (nrow(later) > 0L && all(later$top_group == final_group)) {
      return(as.numeric(k))
    }
  }

  NA_real_
}

tail_agreement_with_final <- function(dt, final_group, tail_fraction) {
  dt2 <- dt[status %in% c("ok", "low_gap") & !is.na(top_group)]

  if (nrow(dt2) == 0L || is.na(final_group)) return(NA_real_)

  setorder(dt2, K)
  n_tail <- max(1L, ceiling(nrow(dt2) * tail_fraction))
  tail_dt <- tail(dt2, n_tail)

  mean(tail_dt$top_group == final_group)
}

summarise_multiK_one_sample <- function(raw_dt, tail_fraction) {
  dt_all <- raw_dt
  dt_ok <- raw_dt[status %in% c("ok", "low_gap") & !is.na(top_group)]

  nK_total <- length(unique(dt_all$K[!is.na(dt_all$K)]))
  nK <- length(unique(dt_ok$K))
  maxK_used <- if (nK > 0L) max(dt_ok$K) else NA_real_

  final_group <- if (!is.na(maxK_used)) {
    dt_ok[K == maxK_used][1, top_group]
  } else {
    NA_character_
  }

  agreement <- if (!is.na(final_group) && nK > 0L) {
    mean(dt_ok$top_group == final_group)
  } else {
    NA_real_
  }

  stable_from_K <- stable_from_threshold(dt_ok, final_group)
  tail_agreement <- tail_agreement_with_final(
    dt_ok,
    final_group,
    tail_fraction
  )

  maxrow <- if (!is.na(maxK_used)) {
    dt_ok[K == maxK_used][1]
  } else {
    NULL
  }

  data.table(
    final_group = final_group,
    stable_from_K = as.numeric(stable_from_K),
    agreement_frac = as.numeric(agreement),
    tail_agreement_frac = as.numeric(tail_agreement),
    tail_fraction = as.numeric(tail_fraction),
    n_K_total_tested = as.integer(nK_total),
    n_K_usable_for_stability = as.integer(nK),
    n_K_available = as.integer(nK),
    maxK_used = as.numeric(maxK_used),
    maxK_top_posterior = if (!is.null(maxrow)) as.numeric(maxrow$top_posterior) else NA_real_,
    maxK_second_posterior = if (!is.null(maxrow)) as.numeric(maxrow$second_posterior) else NA_real_,
    maxK_gap = if (!is.null(maxrow)) as.numeric(maxrow$posterior_gap) else NA_real_,
    maxK_confidence = if (!is.null(maxrow)) as.character(maxrow$confidence) else NA_character_,
    maxK_status = if (!is.null(maxrow)) as.character(maxrow$status) else NA_character_,
    maxK_n_snps_used = if (!is.null(maxrow)) as.integer(maxrow$n_snps_used) else NA_integer_,
    status = if (nK > 0L) "ok" else "no_valid_K"
  )
}

final_call <- function(
    stability_row,
    high_thr,
    mod_thr,
    min_gap,
    min_agreement,
    min_K,
    min_tail_agreement
) {
  if (is.null(stability_row) || nrow(stability_row) == 0L) {
    return(list(call = NA_character_, label = "Uncertain", reason = "no_result"))
  }

  s <- stability_row[1]

  if (is.na(s$final_group)) {
    return(list(call = NA_character_, label = "Uncertain", reason = "no_final_group"))
  }

  if (is.na(s$maxK_n_snps_used) || s$maxK_n_snps_used <= 0L) {
    return(list(call = NA_character_, label = "Uncertain", reason = "no_snps_used"))
  }

  p <- s$maxK_top_posterior
  gap <- s$maxK_gap

  if (is.na(p)) {
    return(list(call = NA_character_, label = "Uncertain", reason = "posterior_missing"))
  }

  confidence <- conf_label(p, high = high_thr, mod = mod_thr)

  if (!(confidence %in% c("High", "Moderate"))) {
    return(list(call = s$final_group, label = "Uncertain", reason = "low_posterior"))
  }

  if (!is.na(gap) && gap < min_gap) {
    return(list(
      call = s$final_group,
      label = "Uncertain",
      reason = paste0("gap<", min_gap)
    ))
  }

  if (is.na(s$n_K_available) || s$n_K_available < min_K) {
    return(list(
      call = s$final_group,
      label = "Uncertain",
      reason = "too_few_usable_K"
    ))
  }

  if (is.na(s$agreement_frac) || s$agreement_frac < min_agreement) {
    return(list(
      call = s$final_group,
      label = "Uncertain",
      reason = "low_global_K_agreement"
    ))
  }

  if (is.na(s$tail_agreement_frac) ||
      s$tail_agreement_frac < min_tail_agreement) {
    return(list(
      call = s$final_group,
      label = "Uncertain",
      reason = "unstable_large_K_tail"
    ))
  }

  list(call = s$final_group, label = confidence, reason = "ok")
}

# =============================================================================
# Fold construction and evaluation
# =============================================================================

make_fold_tasks <- function(ref_meta, fold_units) {
  tasks <- list()

  for (unit in fold_units) {
    if (unit == "country") {
      fold_keys <- ref_meta$fold_country
    } else if (unit == "site") {
      fold_keys <- paste(
        ref_meta$fold_country,
        ref_meta$fold_site,
        sep = "::"
      )
    } else {
      stop("Unsupported fold unit: ", unit)
    }

    keys <- sort(unique(fold_keys))

    for (key in keys) {
      idx <- which(fold_keys == key)

      tasks[[length(tasks) + 1L]] <- list(
        validation_unit = unit,
        fold_id = paste(unit, key, sep = "::"),
        heldout_ids = ref_meta$sample_id[idx]
      )
    }
  }

  tasks
}

make_non_evaluable_result <- function(
    task,
    cfg,
    test_meta,
    fold_row,
    training_counts,
    reason,
    n_panel_snps
) {
  predictions <- test_meta[, .(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label,
    group_col = cfg$group_col,
    validation_unit = task$validation_unit,
    fold_id = task$fold_id,
    heldout_country = collapse_values(country),
    heldout_site = collapse_values(site),
    sample_id,
    sample_country = country,
    sample_site = site,
    true_group = get(cfg$group_col),
    fold_evaluable = FALSE,
    fold_reason = reason,
    predicted_group = NA_character_,
    final_confidence = "Not_evaluable",
    final_reason = reason,
    accepted = NA,
    correct_call = NA,
    correct_accepted = NA,
    incorrect_accepted = NA,
    unresolved = NA,
    final_group = NA_character_,
    maxK_used = NA_real_,
    maxK_top_posterior = NA_real_,
    maxK_second_posterior = NA_real_,
    maxK_gap = NA_real_,
    maxK_confidence = NA_character_,
    maxK_status = NA_character_,
    maxK_n_snps_used = NA_integer_,
    agreement_frac = NA_real_,
    tail_agreement_frac = NA_real_,
    tail_fraction = NA_real_,
    stable_from_K = NA_real_,
    n_K_total_tested = NA_integer_,
    n_K_usable_for_stability = NA_integer_,
    n_K_available = NA_integer_,
    status = "fold_not_evaluable"
  )]

  raw <- predictions[, .(
    panel_id,
    panel_label,
    group_col,
    validation_unit,
    fold_id,
    sample_id,
    true_group,
    fold_evaluable,
    fold_reason,
    K = NA_integer_,
    top_group = NA_character_,
    top_posterior = NA_real_,
    second_group = NA_character_,
    second_posterior = NA_real_,
    posterior_gap = NA_real_,
    n_snps_used = NA_integer_,
    confidence = "Not_evaluable",
    status = "fold_not_evaluable",
    note = reason
  )]

  marker_stats <- data.table(
    panel_id = cfg$panel_id,
    validation_unit = task$validation_unit,
    fold_id = task$fold_id,
    fold_evaluable = FALSE,
    fold_reason = reason,
    n_panel_snps = as.integer(n_panel_snps),
    n_fold_eligible_snps = NA_integer_,
    n_positive_scores = NA_integer_,
    n_zero_scores = NA_integer_,
    n_negative_scores = NA_integer_,
    maximum_fold_score = NA_real_,
    median_fold_score = NA_real_,
    minimum_fold_score = NA_real_,
    ranking_scope = if (isTRUE(settings$foldwise_reranking)) {
      "foldwise_reranking_of_step02_retained_snps"
    } else {
      "fixed_step02_ranking"
    }
  )

  list(
    predictions = predictions,
    raw = raw,
    fold = fold_row,
    training_counts = training_counts,
    marker_stats = marker_stats,
    top_markers = data.table(),
    K_grid = data.table()
  )
}

process_one_fold <- function(
    task,
    X,
    snp_map,
    ref_meta,
    cfg,
    settings,
    baseK,
    extraK
) {
  test_meta <- ref_meta[sample_id %in% task$heldout_ids]
  train_meta <- ref_meta[!(sample_id %in% task$heldout_ids)]

  full_groups <- sort(unique(ref_meta[[cfg$group_col]]))
  true_groups <- sort(unique(test_meta[[cfg$group_col]]))

  train_counts_observed <- train_meta[, .N, by = c(cfg$group_col)]
  setnames(train_counts_observed, cfg$group_col, "group")

  training_counts <- merge(
    data.table(group = full_groups),
    train_counts_observed,
    by = "group",
    all.x = TRUE,
    sort = TRUE
  )

  training_counts[is.na(N), N := 0L]
  training_counts[, supported := N >= settings$min_ref_per_group]
  training_counts[, `:=`(
    panel_id = cfg$panel_id,
    validation_unit = task$validation_unit,
    fold_id = task$fold_id
  )]

  setcolorder(
    training_counts,
    c(
      "panel_id",
      "validation_unit",
      "fold_id",
      "group",
      "N",
      "supported"
    )
  )

  supported_groups <- training_counts[supported %in% TRUE, group]

  fold_reason <- "ok"

  if (length(supported_groups) < 2L) {
    fold_reason <- "fewer_than_two_supported_training_classes"
  } else if (!all(true_groups %in% supported_groups)) {
    fold_reason <- "true_class_absent_or_undersupported_after_holdout"
  }

  fold_evaluable <- identical(fold_reason, "ok")

  fold_row <- data.table(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label,
    group_col = cfg$group_col,
    validation_unit = task$validation_unit,
    fold_id = task$fold_id,
    heldout_country = collapse_values(test_meta$country),
    heldout_site = collapse_values(test_meta$site),
    heldout_true_groups = collapse_values(true_groups),
    n_test = nrow(test_meta),
    n_train = nrow(train_meta),
    n_training_classes_observed = nrow(train_counts_observed),
    n_supported_training_classes = length(supported_groups),
    supported_training_classes = collapse_values(supported_groups),
    fold_evaluable = fold_evaluable,
    fold_reason = fold_reason
  )

  if (!fold_evaluable) {
    return(make_non_evaluable_result(
      task = task,
      cfg = cfg,
      test_meta = test_meta,
      fold_row = fold_row,
      training_counts = training_counts,
      reason = fold_reason,
      n_panel_snps = nrow(X)
    ))
  }

  train_supported <- train_meta[get(cfg$group_col) %in% supported_groups]

  ranked_snps <- score_fold_snps(
    X = X,
    train_meta = train_supported,
    group_col = cfg$group_col,
    panel_type = cfg$panel_type,
    snp_map = snp_map,
    max_site_missing = settings$max_site_missing,
    min_mac = settings$min_mac,
    foldwise_reranking = settings$foldwise_reranking
  )

  if (nrow(ranked_snps) == 0L) {
    fold_reason <- "no_fold_eligible_snps"
    fold_row[, `:=`(
      fold_evaluable = FALSE,
      fold_reason = fold_reason
    )]

    return(make_non_evaluable_result(
      task = task,
      cfg = cfg,
      test_meta = test_meta,
      fold_row = fold_row,
      training_counts = training_counts,
      reason = fold_reason,
      n_panel_snps = nrow(X)
    ))
  }

  snp_order <- ranked_snps$snp_id

  Ks <- make_Ks(
    n_snps = length(snp_order),
    baseK = baseK,
    extraK = extraK,
    auto_extend = settings$auto_extend_K,
    include_all = settings$include_all_snps_K,
    max_points = settings$max_K_points
  )

  if (length(Ks) == 0L) {
    fold_reason <- "empty_K_grid"
    fold_row[, `:=`(
      fold_evaluable = FALSE,
      fold_reason = fold_reason
    )]

    return(make_non_evaluable_result(
      task = task,
      cfg = cfg,
      test_meta = test_meta,
      fold_row = fold_row,
      training_counts = training_counts,
      reason = fold_reason,
      n_panel_snps = nrow(X)
    ))
  }

  X_ranked <- X[snp_order, , drop = FALSE]

  ref_freqs <- build_ref_freqs(
    snp_mat = X_ranked,
    meta_ref = train_supported,
    group_col = cfg$group_col,
    pseudocount = settings$pseudocount
  )

  if (is.null(ref_freqs$freq) || length(ref_freqs$groups) < 2L) {
    fold_reason <- "frequency_estimation_failed"
    fold_row[, `:=`(
      fold_evaluable = FALSE,
      fold_reason = fold_reason
    )]

    return(make_non_evaluable_result(
      task = task,
      cfg = cfg,
      test_meta = test_meta,
      fold_row = fold_row,
      training_counts = training_counts,
      reason = fold_reason,
      n_panel_snps = nrow(X)
    ))
  }

  raw_list <- vector("list", nrow(test_meta))
  prediction_list <- vector("list", nrow(test_meta))

  for (i in seq_len(nrow(test_meta))) {
    sample_id <- test_meta$sample_id[i]
    true_group <- test_meta[[cfg$group_col]][i]
    raw_K_list <- vector("list", length(Ks))

    for (j in seq_along(Ks)) {
      K <- Ks[j]
      snp_ids_K <- snp_order[seq_len(min(K, length(snp_order)))]

      raw_K_list[[j]] <- assign_one_K(
        snp_mat = X_ranked[snp_ids_K, , drop = FALSE],
        freq = ref_freqs$freq[snp_ids_K, , drop = FALSE],
        groups = ref_freqs$groups,
        sample_id = sample_id,
        panel_id = cfg$panel_id,
        group_col = cfg$group_col,
        K = K,
        epsilon = settings$epsilon,
        min_snps_query = cfg$min_snps_query,
        high_thr = settings$high_posterior,
        mod_thr = settings$moderate_posterior,
        min_gap = cfg$min_gap
      )
    }

    raw_i <- rbindlist(raw_K_list, use.names = TRUE, fill = TRUE)
    stability_i <- summarise_multiK_one_sample(
      raw_i,
      tail_fraction = settings$tail_fraction
    )

    final_i <- final_call(
      stability_row = stability_i,
      high_thr = settings$high_posterior,
      mod_thr = settings$moderate_posterior,
      min_gap = cfg$min_gap,
      min_agreement = settings$min_agreement,
      min_K = settings$min_K_available,
      min_tail_agreement = settings$min_tail_agreement
    )

    accepted_i <- final_i$label %in% c("High", "Moderate") &&
      identical(final_i$reason, "ok")

    correct_i <- !is.na(final_i$call) && identical(final_i$call, true_group)

    stability_i[, `:=`(
      panel_id = cfg$panel_id,
      panel_label = cfg$panel_label,
      group_col = cfg$group_col,
      validation_unit = task$validation_unit,
      fold_id = task$fold_id,
      heldout_country = collapse_values(test_meta$country),
      heldout_site = collapse_values(test_meta$site),
      sample_id = sample_id,
      sample_country = test_meta$country[i],
      sample_site = test_meta$site[i],
      true_group = true_group,
      fold_evaluable = TRUE,
      fold_reason = "ok",
      predicted_group = final_i$call,
      final_confidence = final_i$label,
      final_reason = final_i$reason,
      accepted = accepted_i,
      correct_call = correct_i,
      correct_accepted = accepted_i && correct_i,
      incorrect_accepted = accepted_i && !correct_i,
      unresolved = !accepted_i
    )]

    raw_i[, `:=`(
      panel_label = cfg$panel_label,
      validation_unit = task$validation_unit,
      fold_id = task$fold_id,
      true_group = true_group,
      fold_evaluable = TRUE,
      fold_reason = "ok"
    )]

    raw_list[[i]] <- raw_i
    prediction_list[[i]] <- stability_i
  }

  predictions <- rbindlist(prediction_list, use.names = TRUE, fill = TRUE)
  raw <- rbindlist(raw_list, use.names = TRUE, fill = TRUE)

  marker_stats <- ranked_snps[, .(
    panel_id = cfg$panel_id,
    validation_unit = task$validation_unit,
    fold_id = task$fold_id,
    fold_evaluable = TRUE,
    fold_reason = "ok",
    n_panel_snps = as.integer(nrow(X)),
    n_fold_eligible_snps = as.integer(.N),
    n_positive_scores = as.integer(sum(fold_score > 0, na.rm = TRUE)),
    n_zero_scores = as.integer(sum(fold_score == 0, na.rm = TRUE)),
    n_negative_scores = as.integer(sum(fold_score < 0, na.rm = TRUE)),
    maximum_fold_score = safe_max(fold_score),
    median_fold_score = safe_median(fold_score),
    minimum_fold_score = safe_min(fold_score),
    ranking_scope = if (isTRUE(settings$foldwise_reranking)) {
      "foldwise_reranking_of_step02_retained_snps"
    } else {
      "fixed_step02_ranking"
    }
  )]

  top_markers <- if (settings$top_markers_per_fold > 0L) {
    head(ranked_snps, settings$top_markers_per_fold)
  } else {
    data.table()
  }

  if (nrow(top_markers) > 0L) {
    top_markers[, `:=`(
      panel_id = cfg$panel_id,
      validation_unit = task$validation_unit,
      fold_id = task$fold_id
    )]

    setcolorder(
      top_markers,
      c(
        "panel_id",
        "validation_unit",
        "fold_id",
        setdiff(
          names(top_markers),
          c("panel_id", "validation_unit", "fold_id")
        )
      )
    )
  }

  K_grid <- data.table(
    panel_id = cfg$panel_id,
    validation_unit = task$validation_unit,
    fold_id = task$fold_id,
    K = as.integer(Ks),
    n_fold_eligible_snps = as.integer(length(snp_order))
  )

  list(
    predictions = predictions,
    raw = raw,
    fold = fold_row,
    training_counts = training_counts,
    marker_stats = marker_stats,
    top_markers = top_markers,
    K_grid = K_grid
  )
}

# =============================================================================
# Load metadata
# =============================================================================

metadata_path <- file.path(opt$qc_dir, "metadata_clean.tsv")

if (!file.exists(metadata_path)) {
  stop("Missing cleaned metadata: ", metadata_path)
}

meta <- fread(metadata_path)
setnames(meta, tolower(names(meta)))

required_meta <- c(
  "sample_id",
  "is_reference",
  "macroregion_3",
  "subregion",
  "country",
  "site"
)

missing_meta <- setdiff(required_meta, names(meta))

if (length(missing_meta) > 0L) {
  stop(
    "metadata_clean.tsv is missing required column(s): ",
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

if (any(is.na(meta$sample_id))) {
  stop("metadata_clean.tsv contains missing sample identifiers.")
}

if (anyDuplicated(meta$sample_id)) {
  stop("metadata_clean.tsv contains duplicated sample identifiers.")
}

meta[, fold_country := ifelse(
  is.na(country),
  "[country_missing]",
  country
)]

meta[, fold_site := ifelse(
  is.na(site),
  "[site_missing]",
  site
)]

# =============================================================================
# Run one panel
# =============================================================================

run_one_panel <- function(cfg) {
  message("\n============================================================")
  message("Panel: ", cfg$panel_id)
  message("============================================================")

  panel_dir <- file.path(opt$step2_dir, "panels", cfg$panel_id)
  matrix_path <- file.path(panel_dir, "snp_matrix.rds")
  map_path <- file.path(panel_dir, "snp_map.tsv")

  if (!file.exists(matrix_path)) {
    stop("Missing SNP matrix for ", cfg$panel_id, ": ", matrix_path)
  }

  if (!file.exists(map_path)) {
    stop("Missing SNP map for ", cfg$panel_id, ": ", map_path)
  }

  X <- readRDS(matrix_path)

  if (!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "integer"

  if (is.null(rownames(X)) || is.null(colnames(X))) {
    stop(
      "snp_matrix.rds must have SNP row names and sample column names for ",
      cfg$panel_id
    )
  }

  snp_map <- fread(map_path)

  required_map <- c("snp_id", "score")
  missing_map <- setdiff(required_map, names(snp_map))

  if (length(missing_map) > 0L) {
    stop(
      "SNP map for ", cfg$panel_id,
      " is missing column(s): ", paste(missing_map, collapse = ", ")
    )
  }

  if (!("best_pair" %in% names(snp_map))) {
    snp_map[, best_pair := NA_character_]
  }

  if ("global_rank" %in% names(snp_map)) {
    setorder(snp_map, global_rank)
  } else {
    setorder(snp_map, -score)
  }

  common_snps <- snp_map[snp_id %in% rownames(X), snp_id]

  if (length(common_snps) == 0L) {
    stop("No SNP identifiers match between SNP map and matrix for ", cfg$panel_id)
  }

  X <- X[common_snps, , drop = FALSE]
  snp_map <- snp_map[match(common_snps, snp_id)]

  ref_meta <- meta[
    is_reference %in% TRUE &
      macroregion_3 %in% cfg$macroregion_filter &
      !is.na(get(cfg$group_col))
  ]

  missing_samples <- setdiff(ref_meta$sample_id, colnames(X))

  if (length(missing_samples) > 0L) {
    stop(
      "SNP matrix for ", cfg$panel_id,
      " is missing reference sample(s): ",
      paste(head(missing_samples, 20L), collapse = ", ")
    )
  }

  full_group_counts <- ref_meta[, .N, by = c(cfg$group_col)]
  setnames(full_group_counts, cfg$group_col, "group")

  keep_groups <- full_group_counts[
    N >= settings$min_ref_per_group,
    group
  ]

  ref_meta <- ref_meta[get(cfg$group_col) %in% keep_groups]
  full_group_counts <- full_group_counts[group %in% keep_groups]

  if (nrow(full_group_counts) < 2L) {
    stop("Fewer than two eligible classes for ", cfg$panel_id)
  }

  tasks <- make_fold_tasks(ref_meta, fold_units)

  message("Reference samples: ", nrow(ref_meta))
  message(
    "Classes: ",
    paste(
      full_group_counts$group,
      full_group_counts$N,
      sep = "=",
      collapse = ", "
    )
  )
  message("Retained Step 02 SNPs: ", nrow(X))
  message("Grouped folds: ", length(tasks))

  n_workers <- min(opt$n_cores, length(tasks))

  if (n_workers <= 1L) {
    message("Execution mode: serial")

    fold_results <- lapply(
      tasks,
      process_one_fold,
      X = X,
      snp_map = snp_map,
      ref_meta = ref_meta,
      cfg = cfg,
      settings = settings,
      baseK = baseK,
      extraK = extraK
    )

  } else if (.Platform$OS.type == "windows") {
    message("Execution mode: Windows PSOCK with ", n_workers, " workers")

    cl <- parallel::makeCluster(n_workers)

    fold_results <- tryCatch(
      {
        parallel::clusterEvalQ(cl, {
          library(data.table)
          data.table::setDTthreads(1L)
          NULL
        })

        parallel::clusterExport(
          cl,
          varlist = c(
            "X",
            "snp_map",
            "ref_meta",
            "cfg",
            "settings",
            "baseK",
            "extraK",
            "clean_text",
            "collapse_values",
            "safe_mean",
            "safe_median",
            "safe_min",
            "safe_max",
            "conf_label",
            "hudson_fst_vector",
            "score_fold_snps",
            "loglik_one_group",
            "post_from_ll",
            "make_Ks",
            "build_ref_freqs",
            "assign_one_K",
            "stable_from_threshold",
            "tail_agreement_with_final",
            "summarise_multiK_one_sample",
            "final_call",
            "make_non_evaluable_result",
            "process_one_fold"
          ),
          envir = environment()
        )

        parallel::parLapplyLB(
          cl,
          tasks,
          function(task) {
            process_one_fold(
              task = task,
              X = X,
              snp_map = snp_map,
              ref_meta = ref_meta,
              cfg = cfg,
              settings = settings,
              baseK = baseK,
              extraK = extraK
            )
          }
        )
      },
      finally = parallel::stopCluster(cl)
    )

  } else {
    message("Execution mode: forked parallel with ", n_workers, " workers")

    fold_results <- parallel::mclapply(
      tasks,
      function(task) {
        process_one_fold(
          task = task,
          X = X,
          snp_map = snp_map,
          ref_meta = ref_meta,
          cfg = cfg,
          settings = settings,
          baseK = baseK,
          extraK = extraK
        )
      },
      mc.cores = n_workers,
      mc.preschedule = FALSE
    )
  }

  message("Completed folds: ", length(fold_results))

  list(
    predictions = rbindlist(
      lapply(fold_results, `[[`, "predictions"),
      use.names = TRUE,
      fill = TRUE
    ),
    raw = rbindlist(
      lapply(fold_results, `[[`, "raw"),
      use.names = TRUE,
      fill = TRUE
    ),
    folds = rbindlist(
      lapply(fold_results, `[[`, "fold"),
      use.names = TRUE,
      fill = TRUE
    ),
    training_counts = rbindlist(
      lapply(fold_results, `[[`, "training_counts"),
      use.names = TRUE,
      fill = TRUE
    ),
    marker_stats = rbindlist(
      lapply(fold_results, `[[`, "marker_stats"),
      use.names = TRUE,
      fill = TRUE
    ),
    top_markers = rbindlist(
      lapply(fold_results, `[[`, "top_markers"),
      use.names = TRUE,
      fill = TRUE
    ),
    K_grid = rbindlist(
      lapply(fold_results, `[[`, "K_grid"),
      use.names = TRUE,
      fill = TRUE
    ),
    full_group_counts = full_group_counts
  )
}

# =============================================================================
# Run requested panels
# =============================================================================

panel_results <- vector("list", length(panels_requested))
names(panel_results) <- panels_requested

for (i in seq_along(panels_requested)) {
  cfg <- panel_config(panels_requested[i])
  panel_results[[i]] <- run_one_panel(cfg)
}

predictions <- rbindlist(
  lapply(panel_results, `[[`, "predictions"),
  use.names = TRUE,
  fill = TRUE
)

raw <- rbindlist(
  lapply(panel_results, `[[`, "raw"),
  use.names = TRUE,
  fill = TRUE
)

folds <- rbindlist(
  lapply(panel_results, `[[`, "folds"),
  use.names = TRUE,
  fill = TRUE
)

training_counts <- rbindlist(
  lapply(panel_results, `[[`, "training_counts"),
  use.names = TRUE,
  fill = TRUE
)

marker_stats <- rbindlist(
  lapply(panel_results, `[[`, "marker_stats"),
  use.names = TRUE,
  fill = TRUE
)

top_markers <- rbindlist(
  lapply(panel_results, `[[`, "top_markers"),
  use.names = TRUE,
  fill = TRUE
)

K_grid <- rbindlist(
  lapply(panel_results, `[[`, "K_grid"),
  use.names = TRUE,
  fill = TRUE
)

setorder(predictions, panel_id, validation_unit, fold_id, sample_id)
setorder(raw, panel_id, validation_unit, fold_id, sample_id, K)
setorder(folds, panel_id, validation_unit, fold_id)
setorder(training_counts, panel_id, validation_unit, fold_id, group)
setorder(marker_stats, panel_id, validation_unit, fold_id)

if (nrow(top_markers) > 0L) {
  setorder(top_markers, panel_id, validation_unit, fold_id, fold_rank)
}

if (nrow(K_grid) > 0L) {
  setorder(K_grid, panel_id, validation_unit, fold_id, K)
}

# =============================================================================
# Performance summaries
# =============================================================================

fold_performance <- predictions[, {
  evaluable <- fold_evaluable %in% TRUE
  accepted_eval <- accepted[evaluable] %in% TRUE
  correct_eval <- correct_call[evaluable] %in% TRUE
  correct_accepted_eval <- correct_accepted[evaluable] %in% TRUE
  incorrect_accepted_eval <- incorrect_accepted[evaluable] %in% TRUE
  unresolved_eval <- unresolved[evaluable] %in% TRUE

  list(
    panel_label = panel_label[1],
    group_col = group_col[1],
    heldout_country = heldout_country[1],
    heldout_site = heldout_site[1],
    true_groups = collapse_values(true_group),
    fold_evaluable = all(evaluable),
    fold_reason = fold_reason[1],
    n_samples = .N,
    n_evaluable_samples = sum(evaluable),
    n_accepted = sum(accepted_eval),
    n_correct_calls = sum(correct_eval),
    n_correct_accepted = sum(correct_accepted_eval),
    n_incorrect_accepted = sum(incorrect_accepted_eval),
    n_unresolved = sum(unresolved_eval),
    raw_call_accuracy = if (any(evaluable)) safe_mean(correct_eval) else NA_real_,
    coverage = if (any(evaluable)) safe_mean(accepted_eval) else NA_real_,
    accepted_accuracy = if (sum(accepted_eval) > 0L) {
      safe_mean(correct_eval[accepted_eval])
    } else {
      NA_real_
    },
    operational_success_rate = if (any(evaluable)) {
      safe_mean(correct_accepted_eval)
    } else {
      NA_real_
    },
    accepted_error_rate = if (any(evaluable)) {
      safe_mean(incorrect_accepted_eval)
    } else {
      NA_real_
    },
    uncertainty_rate = if (any(evaluable)) {
      safe_mean(unresolved_eval)
    } else {
      NA_real_
    },
    median_maxK_posterior = safe_median(maxK_top_posterior[evaluable]),
    median_maxK_gap = safe_median(maxK_gap[evaluable]),
    median_K_agreement = safe_median(agreement_frac[evaluable])
  )
}, by = .(panel_id, validation_unit, fold_id)]

setorder(fold_performance, panel_id, validation_unit, fold_id)

sample_summary <- predictions[, {
  evaluable <- fold_evaluable %in% TRUE
  accepted_eval <- accepted[evaluable] %in% TRUE
  correct_eval <- correct_call[evaluable] %in% TRUE
  correct_accepted_eval <- correct_accepted[evaluable] %in% TRUE
  incorrect_accepted_eval <- incorrect_accepted[evaluable] %in% TRUE
  unresolved_eval <- unresolved[evaluable] %in% TRUE

  list(
    n_samples = .N,
    n_evaluable_samples = sum(evaluable),
    n_non_evaluable_samples = sum(!evaluable),
    n_accepted = sum(accepted_eval),
    n_correct_calls = sum(correct_eval),
    n_correct_accepted = sum(correct_accepted_eval),
    n_incorrect_accepted = sum(incorrect_accepted_eval),
    n_unresolved = sum(unresolved_eval),
    pooled_raw_call_accuracy = safe_mean(correct_eval),
    pooled_coverage = safe_mean(accepted_eval),
    pooled_accepted_accuracy = if (sum(accepted_eval) > 0L) {
      safe_mean(correct_eval[accepted_eval])
    } else {
      NA_real_
    },
    pooled_operational_success_rate = safe_mean(correct_accepted_eval),
    pooled_accepted_error_rate = safe_mean(incorrect_accepted_eval),
    pooled_uncertainty_rate = safe_mean(unresolved_eval),
    median_maxK_posterior = safe_median(maxK_top_posterior[evaluable]),
    median_maxK_gap = safe_median(maxK_gap[evaluable]),
    median_K_agreement = safe_median(agreement_frac[evaluable])
  )
}, by = .(panel_id, panel_label, group_col, validation_unit)]

fold_summary <- fold_performance[, .(
  n_folds = .N,
  n_evaluable_folds = sum(fold_evaluable %in% TRUE),
  n_non_evaluable_folds = sum(!(fold_evaluable %in% TRUE)),
  fold_balanced_raw_call_accuracy = safe_mean(
    raw_call_accuracy[fold_evaluable %in% TRUE]
  ),
  fold_balanced_coverage = safe_mean(
    coverage[fold_evaluable %in% TRUE]
  ),
  fold_balanced_accepted_accuracy = safe_mean(
    accepted_accuracy[fold_evaluable %in% TRUE]
  ),
  fold_balanced_operational_success_rate = safe_mean(
    operational_success_rate[fold_evaluable %in% TRUE]
  ),
  fold_balanced_accepted_error_rate = safe_mean(
    accepted_error_rate[fold_evaluable %in% TRUE]
  ),
  fold_balanced_uncertainty_rate = safe_mean(
    uncertainty_rate[fold_evaluable %in% TRUE]
  )
), by = .(panel_id, validation_unit)]

panel_summary <- merge(
  sample_summary,
  fold_summary,
  by = c("panel_id", "validation_unit"),
  all = TRUE,
  sort = FALSE
)

setorder(panel_summary, panel_id, validation_unit)

per_class <- predictions[fold_evaluable %in% TRUE, {
  accepted_eval <- accepted %in% TRUE
  correct_eval <- correct_call %in% TRUE

  list(
    n_samples = .N,
    n_folds = uniqueN(fold_id),
    n_accepted = sum(accepted_eval),
    n_correct_calls = sum(correct_eval),
    n_correct_accepted = sum(correct_accepted %in% TRUE),
    n_incorrect_accepted = sum(incorrect_accepted %in% TRUE),
    n_unresolved = sum(unresolved %in% TRUE),
    raw_call_accuracy = safe_mean(correct_eval),
    coverage = safe_mean(accepted_eval),
    accepted_accuracy = if (sum(accepted_eval) > 0L) {
      safe_mean(correct_eval[accepted_eval])
    } else {
      NA_real_
    },
    operational_success_rate = safe_mean(correct_accepted %in% TRUE),
    accepted_error_rate = safe_mean(incorrect_accepted %in% TRUE),
    uncertainty_rate = safe_mean(unresolved %in% TRUE)
  )
}, by = .(panel_id, validation_unit, true_group)]

setorder(per_class, panel_id, validation_unit, true_group)

confusion <- predictions[
  fold_evaluable %in% TRUE & !is.na(predicted_group),
  .N,
  by = .(
    panel_id,
    validation_unit,
    true_group,
    predicted_group,
    accepted
  )
]

if (nrow(confusion) > 0L) {
  confusion[, proportion_within_true := N / sum(N), by = .(
    panel_id,
    validation_unit,
    true_group
  )]

  setorder(
    confusion,
    panel_id,
    validation_unit,
    true_group,
    -N,
    predicted_group
  )
}

reason_counts <- predictions[, .N, by = .(
  panel_id,
  validation_unit,
  fold_evaluable,
  final_reason
)]

setorder(reason_counts, panel_id, validation_unit, -N, final_reason)

non_evaluable_folds <- folds[!(fold_evaluable %in% TRUE)]

# =============================================================================
# Comparable individual LOO summary
# =============================================================================

loo_predictions_path <- file.path(
  opt$loo_dir,
  "step4b_loo_predictions_all_panels.tsv"
)

loo_comparison <- data.table()

if (file.exists(loo_predictions_path)) {
  loo_predictions <- fread(loo_predictions_path)

  required_loo <- c(
    "panel_id",
    "final_confidence",
    "correct"
  )

  if (all(required_loo %in% names(loo_predictions))) {
    loo_predictions[, correct := as_clean_logical(correct)]
    loo_predictions[, accepted := final_confidence %in% c("High", "Moderate")]

    loo_comparison <- loo_predictions[
      panel_id %in% panels_requested,
      {
        accepted_i <- accepted %in% TRUE
        correct_i <- correct %in% TRUE

        list(
          validation_scheme = "individual",
          geographic_holdout = FALSE,
          marker_ranking = "fixed_step02_ranking",
          n_folds = .N,
          n_evaluable_folds = .N,
          n_non_evaluable_folds = 0L,
          n_samples = .N,
          n_evaluable_samples = .N,
          n_non_evaluable_samples = 0L,
          pooled_raw_call_accuracy = safe_mean(correct_i),
          pooled_coverage = safe_mean(accepted_i),
          pooled_accepted_accuracy = if (sum(accepted_i) > 0L) {
            safe_mean(correct_i[accepted_i])
          } else {
            NA_real_
          },
          pooled_operational_success_rate = safe_mean(accepted_i & correct_i),
          pooled_accepted_error_rate = safe_mean(accepted_i & !correct_i),
          pooled_uncertainty_rate = safe_mean(!accepted_i),
          fold_balanced_raw_call_accuracy = safe_mean(correct_i),
          fold_balanced_coverage = safe_mean(accepted_i),
          fold_balanced_accepted_accuracy = if (sum(accepted_i) > 0L) {
            safe_mean(correct_i[accepted_i])
          } else {
            NA_real_
          },
          fold_balanced_operational_success_rate = safe_mean(accepted_i & correct_i),
          fold_balanced_accepted_error_rate = safe_mean(accepted_i & !correct_i),
          fold_balanced_uncertainty_rate = safe_mean(!accepted_i)
        )
      },
      by = panel_id
    ]
  } else {
    warning(
      "Individual LOO comparison skipped because required columns are missing: ",
      paste(setdiff(required_loo, names(loo_predictions)), collapse = ", ")
    )
  }
} else {
  warning(
    "Individual LOO comparison file was not found: ",
    loo_predictions_path
  )
}

grouped_comparison <- panel_summary[, .(
  panel_id,
  validation_scheme = validation_unit,
  geographic_holdout = TRUE,
  marker_ranking = if (isTRUE(settings$foldwise_reranking)) {
    "foldwise_reranking_of_step02_retained_snps"
  } else {
    "fixed_step02_ranking"
  },
  n_folds,
  n_evaluable_folds,
  n_non_evaluable_folds,
  n_samples,
  n_evaluable_samples,
  n_non_evaluable_samples,
  pooled_raw_call_accuracy,
  pooled_coverage,
  pooled_accepted_accuracy,
  pooled_operational_success_rate,
  pooled_accepted_error_rate,
  pooled_uncertainty_rate,
  fold_balanced_raw_call_accuracy,
  fold_balanced_coverage,
  fold_balanced_accepted_accuracy,
  fold_balanced_operational_success_rate,
  fold_balanced_accepted_error_rate,
  fold_balanced_uncertainty_rate
)]

validation_comparison <- rbindlist(
  list(loo_comparison, grouped_comparison),
  use.names = TRUE,
  fill = TRUE
)

scheme_order <- c("individual", "country", "site")
validation_comparison[, scheme_order := match(
  validation_scheme,
  scheme_order
)]
setorder(validation_comparison, panel_id, scheme_order)
validation_comparison[, scheme_order := NULL]

# =============================================================================
# Write outputs
# =============================================================================

write_table_pair <- function(x, stem, compress = FALSE) {
  tsv_path <- file.path(
    opt$out_dir,
    paste0(stem, if (compress) ".tsv.gz" else ".tsv")
  )

  rds_path <- file.path(opt$out_dir, paste0(stem, ".rds"))

  fwrite(
    x,
    tsv_path,
    sep = "\t",
    na = "NA",
    compress = if (compress) "gzip" else "none"
  )

  saveRDS(x, rds_path)

  invisible(list(tsv = tsv_path, rds = rds_path))
}

write_table_pair(predictions, "grouped_cv_predictions")
write_table_pair(raw, "grouped_cv_raw", compress = TRUE)
write_table_pair(folds, "grouped_cv_folds")
write_table_pair(fold_performance, "grouped_cv_fold_performance")
write_table_pair(panel_summary, "grouped_cv_panel_summary")
write_table_pair(per_class, "grouped_cv_per_class")
write_table_pair(confusion, "grouped_cv_confusion")
write_table_pair(reason_counts, "grouped_cv_reason_counts")
write_table_pair(non_evaluable_folds, "grouped_cv_non_evaluable_folds")
write_table_pair(training_counts, "grouped_cv_training_counts")
write_table_pair(marker_stats, "grouped_cv_fold_marker_stats")
write_table_pair(top_markers, "grouped_cv_top_markers")
write_table_pair(K_grid, "grouped_cv_K_grid")
write_table_pair(
  validation_comparison,
  "grouped_cv_comparison_with_individual_loo"
)

params <- data.table(
  parameter = c(
    "panels",
    "fold_units",
    "n_cores",
    "min_ref_per_group",
    "max_site_missing",
    "min_mac",
    "foldwise_reranking",
    "pseudocount",
    "epsilon",
    "min_snps_query_macroregion",
    "min_snps_query_subregion",
    "high_posterior",
    "moderate_posterior",
    "min_gap_macroregion",
    "min_gap_subregion",
    "base_K",
    "extra_K",
    "auto_extend_K",
    "include_all_snps_K",
    "max_K_points",
    "min_agreement",
    "min_K_available",
    "tail_fraction",
    "min_tail_agreement",
    "top_markers_per_fold"
  ),
  value = c(
    paste(panels_requested, collapse = ","),
    paste(fold_units, collapse = ","),
    as.character(opt$n_cores),
    as.character(settings$min_ref_per_group),
    as.character(settings$max_site_missing),
    as.character(settings$min_mac),
    as.character(settings$foldwise_reranking),
    as.character(settings$pseudocount),
    as.character(settings$epsilon),
    as.character(opt$min_snps_query_macroregion),
    as.character(opt$min_snps_query_subregion),
    as.character(settings$high_posterior),
    as.character(settings$moderate_posterior),
    as.character(opt$min_gap_macroregion),
    as.character(opt$min_gap_subregion),
    paste(baseK, collapse = ","),
    paste(extraK, collapse = ","),
    as.character(settings$auto_extend_K),
    as.character(settings$include_all_snps_K),
    as.character(settings$max_K_points),
    as.character(settings$min_agreement),
    as.character(settings$min_K_available),
    as.character(settings$tail_fraction),
    as.character(settings$min_tail_agreement),
    as.character(settings$top_markers_per_fold)
  )
)

write_table_pair(params, "grouped_cv_params")

# =============================================================================
# Workbook with guarded empty table handling
# =============================================================================

add_table_sheet <- function(wb, sheet_name, x) {
  sheet_name <- substr(sheet_name, 1L, 31L)
  addWorksheet(wb, sheet_name)

  if (ncol(x) == 0L || nrow(x) == 0L) {
    writeData(
      wb,
      sheet = sheet_name,
      x = data.frame(note = "No rows for this output"),
      startRow = 1L,
      startCol = 1L
    )
  } else {
    writeDataTable(
      wb,
      sheet = sheet_name,
      x = as.data.frame(x),
      startRow = 1L,
      startCol = 1L,
      tableStyle = "TableStyleMedium2"
    )
  }

  freezePane(wb, sheet = sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet = sheet_name, cols = 1:max(1L, ncol(x)), widths = "auto")
}

workbook_path <- file.path(
  opt$out_dir,
  "grouped_geographic_cv_summary.xlsx"
)

wb <- createWorkbook()

add_table_sheet(wb, "Validation_comparison", validation_comparison)
add_table_sheet(wb, "Panel_summary", panel_summary)
add_table_sheet(wb, "Fold_performance", fold_performance)
add_table_sheet(wb, "Per_class", per_class)
add_table_sheet(wb, "Confusion", confusion)
add_table_sheet(wb, "Non_evaluable_folds", non_evaluable_folds)
add_table_sheet(wb, "Reason_counts", reason_counts)
add_table_sheet(wb, "Marker_stats", marker_stats)
add_table_sheet(wb, "Parameters", params)

saveWorkbook(wb, workbook_path, overwrite = TRUE)

# =============================================================================
# Run record
# =============================================================================

run_info <- list(
  analysis = "grouped_geographic_cross_validation",
  purpose = paste(
    "Evaluate geographic generalisation by withholding complete countries",
    "and collection sites"
  ),
  panels = panels_requested,
  validation_units = fold_units,
  fold_counts = folds[, .N, by = .(panel_id, validation_unit)],
  non_evaluable_folds = non_evaluable_folds,
  candidate_marker_source = paste(
    "Corrected Step 02 retained SNPs, with at most one retained SNP per",
    "orthogroup"
  ),
  marker_ranking_scope = if (isTRUE(settings$foldwise_reranking)) {
    paste(
      "Training fold missingness, minor allele count, and standard Hudson",
      "FST scores are recalculated before foldwise reranking"
    )
  } else {
    "Corrected fixed Step 02 ranking"
  },
  nesting_limitation = paste(
    "The analysis does not return to the original alignments to rediscover",
    "alternative within orthogroup SNPs and is not described as fully nested",
    "de novo marker discovery"
  ),
  hudson_estimator = list(
    numerator = paste(
      "(p1-p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1)"
    ),
    denominator = "p1(1-p2) + p2(1-p1)",
    negative_values = "retained",
    ranking = paste(
      "binary score for P1; maximum pairwise score for multiclass P2 and P3"
    )
  ),
  class_absence_policy = paste(
    "A fold is not evaluable when the true held out class has fewer than",
    settings$min_ref_per_group,
    "training references after withholding"
  ),
  assignment_model = paste(
    "Unchanged Step 04 likelihood, equal class priors, multi K confidence,",
    "posterior gap, global agreement, and tail agreement rules"
  ),
  weighting = paste(
    "Both pooled sample weighted and fold balanced performance estimates",
    "are reported"
  ),
  metadata = metadata_path,
  corrected_step2_dir = opt$step2_dir,
  individual_loo_predictions = if (file.exists(loo_predictions_path)) {
    loo_predictions_path
  } else {
    NA_character_
  },
  output_directory = opt$out_dir,
  workbook = workbook_path,
  parameters = as.list(setNames(params$value, params$parameter)),
  parallel = list(
    requested_workers = opt$n_cores,
    detected_logical_cores = detected_cores,
    platform = .Platform$OS.type,
    windows_backend = "PSOCK"
  ),
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(
  run_info,
  file.path(opt$out_dir, "grouped_cv_run_info.rds")
)

message("\nDone.")
message(
  "Main comparison: ",
  file.path(
    opt$out_dir,
    "grouped_cv_comparison_with_individual_loo.tsv"
  )
)
message(
  "Panel summary:   ",
  file.path(opt$out_dir, "grouped_cv_panel_summary.tsv")
)
message(
  "Fold summary:    ",
  file.path(opt$out_dir, "grouped_cv_fold_performance.tsv")
)
message("Workbook:        ", workbook_path)
message("\nValidation comparison:")
print(validation_comparison)
