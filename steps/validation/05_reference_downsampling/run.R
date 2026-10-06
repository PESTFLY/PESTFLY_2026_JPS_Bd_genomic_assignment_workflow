#!/usr/bin/env Rscript

# PESTFLY supplementary analysis: Balanced reference downsampling
#
# Purpose
# Repeat balanced sampling of class references at multiple training sizes, using independently
# held out reference test sets and thirty replicates per design.
#
# Interpretation
# Training only filtering, Hudson ranking and allele frequency estimation are repeated within
# the retained candidate SNP set. The analysis does not rediscover alternative ortholog sites.
# Raw accuracy and uncertainty rate should be read together; high accuracy among calls does not
# imply that all specimens were reportable.
#
# Technical notes
# Base seed is 20260930. The seed rule adds panel index times 1000000, training size index times
# 10000 and replicate. Default requested sizes are 3, 5, 10, 20, 40, 80 and 100 references per
# class, plus each panel maximum evaluable balanced size. Class size constraints limit the
# actual designs, recorded in downsampling_design.tsv.
#
# Run from the repository root:
#   Rscript steps/validation/05_reference_downsampling/run.R
# See the adjacent README.md for inputs, outputs and complete CLI defaults.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(parallel)
})

# =============================================================================
# Repository and general helpers
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

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::sd(x)
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::quantile(
    x,
    probs = probability,
    na.rm = TRUE,
    names = FALSE,
    type = 7
  ))
}

safe_cor <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)
  if (length(unique(x[keep])) < 2L || length(unique(y[keep])) < 2L) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

mode_value <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)

  tab <- sort(table(x), decreasing = TRUE)
  candidates <- names(tab)[tab == max(tab)]
  sort(candidates)[1]
}

write_table_pair <- function(x, stem, out_dir, compress = FALSE) {
  x <- as.data.table(x)
  extension <- if (isTRUE(compress)) ".tsv.gz" else ".tsv"
  tsv_path <- file.path(out_dir, paste0(stem, extension))
  rds_path <- file.path(out_dir, paste0(stem, ".rds"))

  fwrite(
    x,
    tsv_path,
    sep = "\t",
    na = "NA",
    quote = FALSE,
    compress = if (isTRUE(compress)) "gzip" else "none"
  )

  saveRDS(x, rds_path)
  invisible(x)
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
    "--assignment_dir",
    type = "character",
    default = "results/04_origin_assignment",
    help = "Corrected Step 04 assignment directory [default %default]"
  ),
  make_option(
    "--loo_dir",
    type = "character",
    default = "results/05_loo_validation",
    help = "Corrected individual LOO directory [default %default]"
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/validation/05_reference_downsampling",
    help = "Supplementary validation output directory [default %default]"
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
    "--n_replicates",
    type = "integer",
    default = 30L,
    help = "Random replicates per panel and training size [default %default]"
  ),
  make_option(
    "--n_cores",
    type = "integer",
    default = 0L,
    help = "Parallel workers; zero uses up to eight detected cores minus one [default %default]"
  ),
  make_option(
    "--seed",
    type = "integer",
    default = 20260930L,
    help = "Base random seed [default %default]"
  ),
  make_option(
    "--training_sizes",
    type = "character",
    default = "3,5,10,20,40,80,100",
    help = "Candidate references per class; each panel also includes its maximum evaluable balanced size [default %default]"
  ),
  make_option(
    "--test_per_group",
    type = "integer",
    default = 5L,
    help = "Maximum balanced held out references per class and replicate [default %default]"
  ),
  make_option(
    "--min_ref_per_group",
    type = "integer",
    default = 3L,
    help = "Minimum training references per class [default %default]"
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
    help = "Minimum training minor allele count [default %default]"
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
    help = "Minimum global agreement across usable K values [default %default]"
  ),
  make_option(
    "--min_K_available",
    type = "integer",
    default = 6L,
    help = "Minimum usable K values for a reliable call [default %default]"
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
    "--marker_overlap_K",
    type = "character",
    default = "20,100,500,2000,5000",
    help = "K values for overlap with the complete Step 02 ranking [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$qc_dir <- resolve_path(opt$qc_dir)
opt$step2_dir <- resolve_path(opt$step2_dir)
opt$assignment_dir <- resolve_path(opt$assignment_dir)
opt$loo_dir <- resolve_path(opt$loo_dir)
opt$out_dir <- resolve_path(opt$out_dir)
opt$n_replicates <- max(1L, as.integer(opt$n_replicates))
opt$test_per_group <- max(1L, as.integer(opt$test_per_group))
opt$min_ref_per_group <- max(2L, as.integer(opt$min_ref_per_group))
opt$seed <- as.integer(opt$seed)

panels_requested <- clean_text(strsplit(opt$panels, ",")[[1]])
panels_requested <- panels_requested[!is.na(panels_requested)]

candidate_training_sizes <- parse_integer_vector(opt$training_sizes)
baseK <- parse_integer_vector(opt$base_K)
extraK <- parse_integer_vector(opt$extra_K)
marker_overlap_K <- parse_integer_vector(opt$marker_overlap_K)

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

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

if (is.na(opt$n_cores) || opt$n_cores <= 0L) {
  n_workers <- min(8L, max(1L, as.integer(detected_cores) - 1L))
} else {
  n_workers <- max(1L, as.integer(opt$n_cores))
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
return_dir <- file.path(opt$out_dir, "return_bundle")
dir.create(return_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Panel definitions
# =============================================================================

panel_config <- function(panel_id) {
  if (identical(panel_id, "P1_macroregion_africa_vs_asia")) {
    return(list(
      panel_id = panel_id,
      panel_label = "P1 Africa versus Asia",
      panel_type = "binary",
      group_col = "macroregion_3",
      macroregion_filter = NA_character_,
      allowed_groups = c("Africa", "Asia"),
      query_panel = "P1_macroregion",
      min_snps_query = as.integer(opt$min_snps_query_macroregion),
      min_gap = as.numeric(opt$min_gap_macroregion)
    ))
  }

  if (identical(panel_id, "P2_subregion_within_africa")) {
    return(list(
      panel_id = panel_id,
      panel_label = "P2 Africa subregion",
      panel_type = "multiclass",
      group_col = "subregion",
      macroregion_filter = "Africa",
      allowed_groups = NULL,
      query_panel = "P2_Africa",
      min_snps_query = as.integer(opt$min_snps_query_subregion),
      min_gap = as.numeric(opt$min_gap_subregion)
    ))
  }

  if (identical(panel_id, "P3_subregion_within_asia")) {
    return(list(
      panel_id = panel_id,
      panel_label = "P3 Asia subregion",
      panel_type = "multiclass",
      group_col = "subregion",
      macroregion_filter = "Asia",
      allowed_groups = NULL,
      query_panel = "P3_Asia",
      min_snps_query = as.integer(opt$min_snps_query_subregion),
      min_gap = as.numeric(opt$min_gap_subregion)
    ))
  }

  stop("Unknown panel identifier: ", panel_id)
}

# =============================================================================
# Corrected Hudson scoring and training only marker ranking
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

  # Negative finite estimates are deliberately retained.
  out
}

score_training_snps <- function(
    X,
    train_meta,
    group_col,
    panel_type,
    snp_map,
    max_site_missing,
    min_mac
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
  groups <- groups[!is.na(groups)]
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

  training_score <- rep(NA_real_, nrow(X))
  best_pair <- rep(NA_character_, nrow(X))

  if (identical(panel_type, "binary")) {
    g1 <- groups[1]
    g2 <- groups[2]

    training_score <- hudson_fst_vector(
      pmat[, g1],
      pmat[, g2],
      nmat[, g1],
      nmat[, g2]
    )

    best_pair[is.finite(training_score)] <- paste(g1, g2, sep = "__vs__")
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
          (!is.finite(training_score) | score_ab > training_score)

        training_score[take] <- score_ab[take]
        best_pair[take] <- paste(g1, g2, sep = "__vs__")
      }
    }
  }

  keep <- eligible_basic & is.finite(training_score)
  if (!any(keep)) return(data.table())

  out <- data.table(
    snp_id = rownames(X)[keep],
    training_score = training_score[keep],
    training_best_pair = best_pair[keep],
    training_missing_frac = missing_frac[keep],
    training_mac = as.integer(mac[keep]),
    training_n_nonmissing = as.integer(n_nonmissing[keep])
  )

  idx_map <- match(out$snp_id, snp_map$snp_id)
  out[, `:=`(
    og = if ("og" %in% names(snp_map)) snp_map$og[idx_map] else NA_character_,
    original_score = if ("score" %in% names(snp_map)) {
      as.numeric(snp_map$score[idx_map])
    } else {
      NA_real_
    },
    original_global_rank = if ("global_rank" %in% names(snp_map)) {
      as.integer(snp_map$global_rank[idx_map])
    } else {
      NA_integer_
    }
  )]

  setorder(
    out,
    -training_score,
    training_missing_frac,
    -training_mac,
    snp_id
  )

  out[, training_rank := seq_len(.N)]
  out
}

# =============================================================================
# Likelihood and multi K helpers
# =============================================================================

conf_label <- function(pmax, high = 0.95, mod = 0.85) {
  if (is.na(pmax)) return("NA")
  if (pmax >= high) return("High")
  if (pmax >= mod) return("Moderate")
  "Low"
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

build_ref_freqs <- function(snp_mat, meta_ref, group_col, pseudocount = 0.5) {
  ref_ids <- intersect(meta_ref$sample_id, colnames(snp_mat))
  if (length(ref_ids) < 2L) {
    return(list(freq = NULL, groups = character(), ref_ids = ref_ids))
  }

  meta_ref2 <- meta_ref[sample_id %in% ref_ids]
  meta_ref2 <- meta_ref2[!is.na(get(group_col))]

  groups <- sort(unique(meta_ref2[[group_col]]))
  groups <- groups[!is.na(groups)]

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

stable_from_threshold <- function(dt, final_group) {
  dt2 <- dt[status %in% c("ok", "low_gap") & !is.na(top_group)]
  if (nrow(dt2) == 0L || is.na(final_group)) return(NA_real_)

  dt2 <- dt2[order(K)]
  for (k in sort(unique(dt2$K))) {
    later <- dt2[K >= k]
    if (nrow(later) > 0L && all(later$top_group == final_group)) {
      return(as.numeric(k))
    }
  }

  NA_real_
}

tail_agreement_with_final <- function(dt, final_group, tail_fraction = 0.50) {
  dt2 <- dt[status %in% c("ok", "low_gap") & !is.na(top_group)]
  if (nrow(dt2) == 0L || is.na(final_group)) return(NA_real_)

  dt2 <- dt2[order(K)]
  n_tail <- max(1L, ceiling(nrow(dt2) * tail_fraction))
  mean(tail(dt2, n_tail)$top_group == final_group)
}

summarise_multiK_one_sample <- function(raw_dt, tail_fraction) {
  dt_all <- raw_dt
  dt_ok <- raw_dt[status %in% c("ok", "low_gap") & !is.na(top_group)]

  nK_total <- length(unique(dt_all$K))
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

  max_row <- if (!is.na(maxK_used)) dt_ok[K == maxK_used][1] else NULL

  data.table(
    final_group = final_group,
    stable_from_K = stable_from_threshold(dt_ok, final_group),
    agreement_frac = as.numeric(agreement),
    tail_agreement_frac = tail_agreement_with_final(
      dt_ok,
      final_group,
      tail_fraction
    ),
    tail_fraction = as.numeric(tail_fraction),
    n_K_total_tested = as.integer(nK_total),
    n_K_usable_for_stability = as.integer(nK),
    n_K_available = as.integer(nK),
    maxK_used = as.numeric(maxK_used),
    maxK_top_posterior = if (!is.null(max_row)) {
      as.numeric(max_row$top_posterior)
    } else {
      NA_real_
    },
    maxK_second_posterior = if (!is.null(max_row)) {
      as.numeric(max_row$second_posterior)
    } else {
      NA_real_
    },
    maxK_gap = if (!is.null(max_row)) {
      as.numeric(max_row$posterior_gap)
    } else {
      NA_real_
    },
    maxK_confidence = if (!is.null(max_row)) {
      as.character(max_row$confidence)
    } else {
      NA_character_
    },
    maxK_status = if (!is.null(max_row)) {
      as.character(max_row$status)
    } else {
      NA_character_
    },
    maxK_n_snps_used = if (!is.null(max_row)) {
      as.integer(max_row$n_snps_used)
    } else {
      NA_integer_
    },
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
    return(list(call = NA_character_, label = "Uncertain", reason = "posterior_NA"))
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

assign_one_sample_fast <- function(
    sample_id,
    X_ranked,
    freq,
    groups,
    Ks,
    cfg,
    settings
) {
  qvec <- X_ranked[, sample_id]
  n_markers <- length(qvec)

  cumulative_ll <- matrix(
    0,
    nrow = n_markers,
    ncol = length(groups),
    dimnames = list(NULL, groups)
  )

  cumulative_group_observed <- matrix(
    0L,
    nrow = n_markers,
    ncol = length(groups),
    dimnames = list(NULL, groups)
  )

  for (g in groups) {
    pvec <- freq[, g]
    ok <- !is.na(qvec) & !is.na(pvec)

    p1 <- pvec * (1 - settings$epsilon) +
      (1 - pvec) * settings$epsilon
    p0 <- (1 - pvec) * (1 - settings$epsilon) +
      pvec * settings$epsilon

    p1 <- pmin(pmax(p1, 1e-12), 1 - 1e-12)
    p0 <- pmin(pmax(p0, 1e-12), 1 - 1e-12)

    contribution <- numeric(n_markers)
    contribution[ok] <- ifelse(
      qvec[ok] == 1L,
      log(p1[ok]),
      log(p0[ok])
    )

    cumulative_ll[, g] <- cumsum(contribution)
    cumulative_group_observed[, g] <- cumsum(ok)
  }

  any_frequency <- rowSums(!is.na(freq)) > 0L
  cumulative_used <- cumsum(!is.na(qvec) & any_frequency)
  raw_rows <- vector("list", length(Ks))

  for (i in seq_along(Ks)) {
    K <- as.integer(Ks[i])
    idx <- min(K, n_markers)
    n_used <- as.integer(cumulative_used[idx])

    if (is.na(n_used) || n_used < cfg$min_snps_query) {
      raw_rows[[i]] <- data.table(
        K = K,
        top_group = NA_character_,
        top_posterior = NA_real_,
        second_group = NA_character_,
        second_posterior = NA_real_,
        posterior_gap = NA_real_,
        n_snps_used = n_used,
        confidence = "Low",
        status = "too_few_snps"
      )
      next
    }

    lls <- cumulative_ll[idx, ]
    lls[cumulative_group_observed[idx, ] == 0L] <- NA_real_
    post <- post_from_ll(lls)

    if (all(is.na(post))) {
      raw_rows[[i]] <- data.table(
        K = K,
        top_group = NA_character_,
        top_posterior = NA_real_,
        second_group = NA_character_,
        second_posterior = NA_real_,
        posterior_gap = NA_real_,
        n_snps_used = n_used,
        confidence = "Low",
        status = "no_likelihood"
      )
      next
    }

    ord <- order(post, decreasing = TRUE)
    top_index <- ord[1]
    second_index <- if (length(ord) >= 2L) ord[2] else NA_integer_

    top_group <- names(post)[top_index]
    top_posterior <- unname(post[top_index])
    second_group <- if (!is.na(second_index)) {
      names(post)[second_index]
    } else {
      NA_character_
    }
    second_posterior <- if (!is.na(second_index)) {
      unname(post[second_index])
    } else {
      NA_real_
    }
    posterior_gap <- if (!is.na(second_posterior)) {
      top_posterior - second_posterior
    } else {
      NA_real_
    }

    raw_rows[[i]] <- data.table(
      K = K,
      top_group = top_group,
      top_posterior = top_posterior,
      second_group = second_group,
      second_posterior = second_posterior,
      posterior_gap = posterior_gap,
      n_snps_used = n_used,
      confidence = conf_label(
        top_posterior,
        high = settings$high_posterior,
        mod = settings$moderate_posterior
      ),
      status = if (!is.na(posterior_gap) && posterior_gap < cfg$min_gap) {
        "low_gap"
      } else {
        "ok"
      }
    )
  }

  raw_dt <- rbindlist(raw_rows, use.names = TRUE, fill = TRUE)
  stability <- summarise_multiK_one_sample(
    raw_dt,
    tail_fraction = settings$tail_fraction
  )

  call <- final_call(
    stability,
    high_thr = settings$high_posterior,
    mod_thr = settings$moderate_posterior,
    min_gap = cfg$min_gap,
    min_agreement = settings$min_agreement,
    min_K = settings$min_K_available,
    min_tail_agreement = settings$min_tail_agreement
  )

  stability[, `:=`(
    sample_id = sample_id,
    predicted_group = call$call,
    final_confidence = call$label,
    final_reason = call$reason,
    accepted = call$label %in% c("High", "Moderate") &&
      identical(call$reason, "ok")
  )]

  stability
}

assign_sample_set <- function(
    sample_ids,
    X_ranked,
    freq,
    groups,
    Ks,
    cfg,
    settings
) {
  sample_ids <- intersect(sample_ids, colnames(X_ranked))
  if (length(sample_ids) == 0L) return(data.table())

  rows <- lapply(
    sample_ids,
    assign_one_sample_fast,
    X_ranked = X_ranked,
    freq = freq,
    groups = groups,
    Ks = Ks,
    cfg = cfg,
    settings = settings
  )

  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

# =============================================================================
# One downsampling replicate
# =============================================================================

make_failed_replicate <- function(task, cfg, reason) {
  replicate_row <- data.table(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label,
    n_per_group = as.integer(task$n_per_group),
    replicate = as.integer(task$replicate),
    seed = as.integer(task$seed),
    status = "not_evaluable",
    reason = reason,
    n_training = NA_integer_,
    n_test = NA_integer_,
    n_eligible_snps = NA_integer_,
    accuracy_all = NA_real_,
    accepted_rate = NA_real_,
    accuracy_accepted = NA_real_,
    correct_accepted_rate = NA_real_,
    uncertainty_rate = NA_real_,
    mean_maxK_posterior = NA_real_,
    mean_maxK_gap = NA_real_,
    mean_agreement = NA_real_,
    mean_tail_agreement = NA_real_
  )

  list(
    replicate_performance = replicate_row,
    test_predictions = data.table(),
    per_class = data.table(),
    query_predictions = data.table(),
    query_replicate = data.table(),
    marker_stats = data.table(),
    marker_overlap = data.table(),
    training_manifest = data.table(),
    test_manifest = data.table(),
    K_grid = data.table()
  )
}

process_one_replicate <- function(
    task,
    X,
    snp_map,
    ref_meta,
    query_baseline,
    cfg,
    settings,
    baseK,
    extraK,
    marker_overlap_K,
    full_snp_order
) {
  set.seed(task$seed)

  groups <- sort(unique(ref_meta[[cfg$group_col]]))
  training_ids <- character()
  test_ids <- character()

  for (g in groups) {
    pool <- ref_meta[get(cfg$group_col) == g, sample_id]

    if (length(pool) < task$n_per_group + task$n_test_per_group) {
      return(make_failed_replicate(
        task,
        cfg,
        paste0("insufficient_samples_in_", g)
      ))
    }

    train_g <- sample(pool, size = task$n_per_group, replace = FALSE)
    remaining <- setdiff(pool, train_g)
    test_g <- sample(
      remaining,
      size = task$n_test_per_group,
      replace = FALSE
    )

    training_ids <- c(training_ids, train_g)
    test_ids <- c(test_ids, test_g)
  }

  train_meta <- ref_meta[sample_id %in% training_ids]
  test_meta <- ref_meta[sample_id %in% test_ids]

  training_manifest <- train_meta[, .(
    panel_id = cfg$panel_id,
    n_per_group = as.integer(task$n_per_group),
    replicate = as.integer(task$replicate),
    seed = as.integer(task$seed),
    sample_id,
    analytical_group = get(cfg$group_col)
  )]

  test_manifest <- test_meta[, .(
    panel_id = cfg$panel_id,
    n_per_group = as.integer(task$n_per_group),
    replicate = as.integer(task$replicate),
    seed = as.integer(task$seed),
    sample_id,
    analytical_group = get(cfg$group_col)
  )]

  ranked_snps <- score_training_snps(
    X = X,
    train_meta = train_meta,
    group_col = cfg$group_col,
    panel_type = cfg$panel_type,
    snp_map = snp_map,
    max_site_missing = settings$max_site_missing,
    min_mac = settings$min_mac
  )

  if (nrow(ranked_snps) == 0L) {
    failed <- make_failed_replicate(task, cfg, "no_training_eligible_snps")
    failed$training_manifest <- training_manifest
    failed$test_manifest <- test_manifest
    return(failed)
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
    failed <- make_failed_replicate(task, cfg, "empty_K_grid")
    failed$training_manifest <- training_manifest
    failed$test_manifest <- test_manifest
    return(failed)
  }

  X_ranked <- X[snp_order, , drop = FALSE]
  ref_freqs <- build_ref_freqs(
    snp_mat = X_ranked,
    meta_ref = train_meta,
    group_col = cfg$group_col,
    pseudocount = settings$pseudocount
  )

  if (is.null(ref_freqs$freq) || length(ref_freqs$groups) < 2L) {
    failed <- make_failed_replicate(task, cfg, "frequency_estimation_failed")
    failed$training_manifest <- training_manifest
    failed$test_manifest <- test_manifest
    return(failed)
  }

  test_predictions <- assign_sample_set(
    sample_ids = test_meta$sample_id,
    X_ranked = X_ranked,
    freq = ref_freqs$freq,
    groups = ref_freqs$groups,
    Ks = Ks,
    cfg = cfg,
    settings = settings
  )

  test_truth <- test_meta[, .(
    sample_id,
    true_group = get(cfg$group_col)
  )]

  test_predictions <- merge(
    test_predictions,
    test_truth,
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  test_predictions[, `:=`(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label,
    n_per_group = as.integer(task$n_per_group),
    replicate = as.integer(task$replicate),
    seed = as.integer(task$seed),
    correct_call = !is.na(predicted_group) & predicted_group == true_group,
    correct_accepted = accepted %in% TRUE &
      !is.na(predicted_group) & predicted_group == true_group,
    incorrect_accepted = accepted %in% TRUE &
      !is.na(predicted_group) & predicted_group != true_group,
    unresolved = !(accepted %in% TRUE)
  )]

  per_class <- test_predictions[, .(
    n_test = .N,
    correct_n = sum(correct_call, na.rm = TRUE),
    accuracy_all = mean(correct_call, na.rm = TRUE),
    accepted_n = sum(accepted, na.rm = TRUE),
    accepted_rate = mean(accepted, na.rm = TRUE),
    accuracy_accepted = if (sum(accepted, na.rm = TRUE) > 0L) {
      mean(correct_call[accepted %in% TRUE], na.rm = TRUE)
    } else {
      NA_real_
    },
    uncertainty_rate = mean(final_confidence == "Uncertain", na.rm = TRUE)
  ), by = .(
    panel_id,
    panel_label,
    n_per_group,
    replicate,
    seed,
    true_group
  )]

  replicate_performance <- test_predictions[, .(
    status = "completed",
    reason = "ok",
    n_training = as.integer(nrow(train_meta)),
    n_test = .N,
    n_eligible_snps = as.integer(length(snp_order)),
    accuracy_all = mean(correct_call, na.rm = TRUE),
    accepted_rate = mean(accepted, na.rm = TRUE),
    accuracy_accepted = if (sum(accepted, na.rm = TRUE) > 0L) {
      mean(correct_call[accepted %in% TRUE], na.rm = TRUE)
    } else {
      NA_real_
    },
    correct_accepted_rate = mean(correct_accepted, na.rm = TRUE),
    uncertainty_rate = mean(final_confidence == "Uncertain", na.rm = TRUE),
    mean_maxK_posterior = safe_mean(maxK_top_posterior),
    mean_maxK_gap = safe_mean(maxK_gap),
    mean_agreement = safe_mean(agreement_frac),
    mean_tail_agreement = safe_mean(tail_agreement_frac)
  ), by = .(
    panel_id,
    panel_label,
    n_per_group,
    replicate,
    seed
  )]

  query_predictions <- assign_sample_set(
    sample_ids = query_baseline$sample_id,
    X_ranked = X_ranked,
    freq = ref_freqs$freq,
    groups = ref_freqs$groups,
    Ks = Ks,
    cfg = cfg,
    settings = settings
  )

  if (nrow(query_predictions) > 0L) {
    query_predictions <- merge(
      query_predictions,
      query_baseline,
      by = "sample_id",
      all.x = TRUE,
      sort = FALSE
    )

    query_predictions[, `:=`(
      panel_id = cfg$panel_id,
      panel_label = cfg$panel_label,
      n_per_group = as.integer(task$n_per_group),
      replicate = as.integer(task$replicate),
      seed = as.integer(task$seed),
      call_matches_baseline = !is.na(predicted_group) &
        !is.na(baseline_call) & predicted_group == baseline_call,
      confidence_matches_baseline = !is.na(final_confidence) &
        !is.na(baseline_confidence) &
        final_confidence == baseline_confidence
    )]

    query_replicate <- query_predictions[, .(
      n_queries = .N,
      calls_matching_baseline = sum(call_matches_baseline, na.rm = TRUE),
      fraction_calls_matching_baseline = mean(
        call_matches_baseline,
        na.rm = TRUE
      ),
      all_calls_match_baseline = all(call_matches_baseline %in% TRUE),
      accepted_fraction = mean(accepted, na.rm = TRUE),
      uncertainty_rate = mean(
        final_confidence == "Uncertain",
        na.rm = TRUE
      )
    ), by = .(
      panel_id,
      panel_label,
      n_per_group,
      replicate,
      seed
    )]
  } else {
    query_replicate <- data.table()
  }

  common_markers <- intersect(full_snp_order, snp_order)
  full_rank <- match(common_markers, full_snp_order)
  training_rank <- match(common_markers, snp_order)

  marker_stats <- data.table(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label,
    n_per_group = as.integer(task$n_per_group),
    replicate = as.integer(task$replicate),
    seed = as.integer(task$seed),
    n_complete_panel_snps = as.integer(length(full_snp_order)),
    n_training_eligible_snps = as.integer(length(snp_order)),
    n_common_snps = as.integer(length(common_markers)),
    common_fraction_of_complete = length(common_markers) /
      length(full_snp_order),
    spearman_rank_correlation_common = safe_cor(
      full_rank,
      training_rank,
      method = "spearman"
    ),
    n_positive_scores = as.integer(sum(
      ranked_snps$training_score > 0,
      na.rm = TRUE
    )),
    n_zero_scores = as.integer(sum(
      ranked_snps$training_score == 0,
      na.rm = TRUE
    )),
    n_negative_scores = as.integer(sum(
      ranked_snps$training_score < 0,
      na.rm = TRUE
    ))
  )

  marker_overlap <- rbindlist(lapply(marker_overlap_K, function(K_requested) {
    K_complete <- min(as.integer(K_requested), length(full_snp_order))
    K_training <- min(as.integer(K_requested), length(snp_order))

    complete_top <- head(full_snp_order, K_complete)
    training_top <- head(snp_order, K_training)
    overlap_n <- length(intersect(complete_top, training_top))
    union_n <- length(union(complete_top, training_top))

    data.table(
      panel_id = cfg$panel_id,
      panel_label = cfg$panel_label,
      n_per_group = as.integer(task$n_per_group),
      replicate = as.integer(task$replicate),
      seed = as.integer(task$seed),
      K_requested = as.integer(K_requested),
      K_complete_used = as.integer(K_complete),
      K_training_used = as.integer(K_training),
      overlap_n = as.integer(overlap_n),
      overlap_fraction_of_complete = if (K_complete > 0L) {
        overlap_n / K_complete
      } else {
        NA_real_
      },
      overlap_fraction_of_training = if (K_training > 0L) {
        overlap_n / K_training
      } else {
        NA_real_
      },
      jaccard = if (union_n > 0L) overlap_n / union_n else NA_real_
    )
  }), use.names = TRUE, fill = TRUE)

  K_grid <- data.table(
    panel_id = cfg$panel_id,
    n_per_group = as.integer(task$n_per_group),
    replicate = as.integer(task$replicate),
    seed = as.integer(task$seed),
    n_training_eligible_snps = as.integer(length(snp_order)),
    K = as.integer(Ks)
  )

  list(
    replicate_performance = replicate_performance,
    test_predictions = test_predictions,
    per_class = per_class,
    query_predictions = query_predictions,
    query_replicate = query_replicate,
    marker_stats = marker_stats,
    marker_overlap = marker_overlap,
    training_manifest = training_manifest,
    test_manifest = test_manifest,
    K_grid = K_grid
  )
}

# =============================================================================
# Load and validate primary inputs
# =============================================================================

metadata_path <- file.path(opt$qc_dir, "metadata_clean.tsv")
assignment_path <- file.path(opt$assignment_dir, "final_assignment.tsv")
baseline_loo_path <- file.path(
  opt$loo_dir,
  "step4b_loo_summary_all_panels.tsv"
)

required_paths <- c(metadata_path, assignment_path, baseline_loo_path)
missing_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_paths) > 0L) {
  stop("Missing required input file(s):\n", paste(missing_paths, collapse = "\n"))
}

meta <- fread(metadata_path)
setnames(meta, tolower(names(meta)))

required_meta <- c("sample_id", "is_reference", "macroregion_3", "subregion")
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta) > 0L) {
  stop("metadata_clean.tsv is missing: ", paste(missing_meta, collapse = ", "))
}

for (nm in names(meta)) {
  if (!inherits(meta[[nm]], c("numeric", "integer", "logical", "Date", "POSIXct"))) {
    meta[[nm]] <- clean_text(meta[[nm]])
  }
}

meta[, `:=`(
  sample_id = clean_text(sample_id),
  is_reference = as_clean_logical(is_reference),
  macroregion_3 = clean_text(macroregion_3),
  subregion = clean_text(subregion)
)]

if (anyDuplicated(meta$sample_id)) {
  stop("metadata_clean.tsv contains duplicated sample_id values.")
}

baseline_assignment <- fread(assignment_path)
setnames(baseline_assignment, tolower(names(baseline_assignment)))

required_assignment <- c(
  "sample_id",
  "macroregion_call",
  "macroregion_confidence",
  "macroregion_reason",
  "subregion_call",
  "subregion_confidence",
  "subregion_reason",
  "subregion_panel"
)

missing_assignment <- setdiff(required_assignment, names(baseline_assignment))
if (length(missing_assignment) > 0L) {
  stop("final_assignment.tsv is missing: ", paste(missing_assignment, collapse = ", "))
}

baseline_loo <- fread(baseline_loo_path)

settings <- list(
  max_site_missing = as.numeric(opt$max_site_missing),
  min_mac = as.integer(opt$min_mac),
  pseudocount = as.numeric(opt$pseudocount),
  epsilon = as.numeric(opt$epsilon),
  high_posterior = as.numeric(opt$high_posterior),
  moderate_posterior = as.numeric(opt$moderate_posterior),
  auto_extend_K = auto_extend_K,
  include_all_snps_K = include_all_snps_K,
  max_K_points = as.integer(opt$max_K_points),
  min_agreement = as.numeric(opt$min_agreement),
  min_K_available = as.integer(opt$min_K_available),
  tail_fraction = as.numeric(opt$tail_fraction),
  min_tail_agreement = as.numeric(opt$min_tail_agreement)
)

message("Repository root: ", repo_root)
message("Corrected Step 02: ", opt$step2_dir)
message("Output directory: ", opt$out_dir)
message("Replicates per training size: ", opt$n_replicates)
message("Parallel workers: ", n_workers)

# =============================================================================
# Panel runs
# =============================================================================

result_names <- c(
  "replicate_performance",
  "test_predictions",
  "per_class",
  "query_predictions",
  "query_replicate",
  "marker_stats",
  "marker_overlap",
  "training_manifest",
  "test_manifest",
  "K_grid"
)

all_results <- setNames(lapply(result_names, function(x) list()), result_names)
design_rows <- list()
reference_count_rows <- list()
baseline_query_rows <- list()
panel_index_rows <- list()

worker_functions <- c(
  "clean_text",
  "safe_mean",
  "safe_cor",
  "hudson_fst_vector",
  "score_training_snps",
  "conf_label",
  "post_from_ll",
  "make_Ks",
  "build_ref_freqs",
  "stable_from_threshold",
  "tail_agreement_with_final",
  "summarise_multiK_one_sample",
  "final_call",
  "assign_one_sample_fast",
  "assign_sample_set",
  "make_failed_replicate",
  "process_one_replicate"
)

for (panel_index in seq_along(panels_requested)) {
  panel_id_value <- panels_requested[panel_index]
  cfg <- panel_config(panel_id_value)

  message("\n============================================================")
  message("Panel: ", cfg$panel_label)
  message("============================================================")

  panel_dir <- file.path(opt$step2_dir, "panels", cfg$panel_id)
  matrix_path <- file.path(panel_dir, "snp_matrix.rds")
  map_path <- file.path(panel_dir, "snp_map.tsv")

  if (!file.exists(matrix_path) || !file.exists(map_path)) {
    stop("Missing corrected Step 02 files for panel: ", cfg$panel_id)
  }

  X <- readRDS(matrix_path)
  if (!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "integer"

  if (is.null(rownames(X)) || is.null(colnames(X))) {
    stop("snp_matrix.rds lacks row or column names for ", cfg$panel_id)
  }

  snp_map <- fread(map_path)
  if (!("snp_id" %in% names(snp_map))) {
    stop("snp_map.tsv lacks snp_id for ", cfg$panel_id)
  }

  if ("global_rank" %in% names(snp_map)) {
    setorder(snp_map, global_rank)
  } else if ("score" %in% names(snp_map)) {
    setorder(snp_map, -score)
  }

  snp_map <- snp_map[snp_id %in% rownames(X)]
  full_snp_order <- snp_map$snp_id

  if (length(full_snp_order) == 0L) {
    stop("No mapped SNPs remain for ", cfg$panel_id)
  }

  if (anyDuplicated(full_snp_order)) {
    stop("Duplicated snp_id values in corrected Step 02 map for ", cfg$panel_id)
  }

  X <- X[full_snp_order, , drop = FALSE]

  ref_meta <- meta[is_reference %in% TRUE]

  if (!is.null(cfg$allowed_groups)) {
    ref_meta <- ref_meta[get(cfg$group_col) %in% cfg$allowed_groups]
  }

  if (!is.na(cfg$macroregion_filter)) {
    ref_meta <- ref_meta[macroregion_3 == cfg$macroregion_filter]
  }

  ref_meta <- ref_meta[!is.na(get(cfg$group_col))]
  ref_meta <- ref_meta[sample_id %in% colnames(X)]

  # c() makes the character column name explicit for recent data.table
  # versions, which no longer accept a bare character scalar in by.
  group_counts <- ref_meta[, .N, by = c(cfg$group_col)]
  setnames(group_counts, cfg$group_col, "analytical_group")
  setorder(group_counts, analytical_group)

  if (nrow(group_counts) < 2L) {
    stop("Fewer than two reference classes for ", cfg$panel_id)
  }

  minimum_class_n <- min(group_counts$N)
  maximum_balanced_training_n <- minimum_class_n - 1L

  if (maximum_balanced_training_n < opt$min_ref_per_group) {
    stop(
      "Panel ", cfg$panel_id,
      " cannot retain an independent test reference per class with at least ",
      opt$min_ref_per_group,
      " training references."
    )
  }

  training_sizes <- sort(unique(c(
    candidate_training_sizes[
      candidate_training_sizes >= opt$min_ref_per_group &
        candidate_training_sizes <= maximum_balanced_training_n
    ],
    maximum_balanced_training_n
  )))

  if (identical(cfg$query_panel, "P1_macroregion")) {
    query_baseline <- baseline_assignment[, .(
      sample_id,
      baseline_call = clean_text(macroregion_call),
      baseline_confidence = clean_text(macroregion_confidence),
      baseline_reason = clean_text(macroregion_reason)
    )]
  } else {
    query_baseline <- baseline_assignment[
      subregion_panel == cfg$query_panel,
      .(
        sample_id,
        baseline_call = clean_text(subregion_call),
        baseline_confidence = clean_text(subregion_confidence),
        baseline_reason = clean_text(subregion_reason)
      )
    ]
  }

  query_baseline <- query_baseline[sample_id %in% colnames(X)]
  if (nrow(query_baseline) == 0L) {
    stop("No applicable baseline queries found for ", cfg$panel_id)
  }

  reference_count_rows[[panel_index]] <- group_counts[, `:=`(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label
  )]

  baseline_query_rows[[panel_index]] <- copy(query_baseline)[, `:=`(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label
  )]

  panel_index_rows[[panel_index]] <- data.table(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label,
    group_col = cfg$group_col,
    n_reference_groups = nrow(group_counts),
    n_reference_samples = nrow(ref_meta),
    minimum_class_n = as.integer(minimum_class_n),
    maximum_balanced_training_n = as.integer(maximum_balanced_training_n),
    training_sizes = paste(training_sizes, collapse = ","),
    n_complete_panel_snps = length(full_snp_order),
    n_applicable_queries = nrow(query_baseline)
  )

  panel_workers <- min(n_workers, opt$n_replicates)
  cl <- NULL

  if (panel_workers > 1L) {
    cl <- parallel::makeCluster(panel_workers)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(data.table))
      NULL
    })

    parallel::clusterExport(
      cl,
      varlist = c(
        worker_functions,
        "X",
        "snp_map",
        "ref_meta",
        "query_baseline",
        "cfg",
        "settings",
        "baseK",
        "extraK",
        "marker_overlap_K",
        "full_snp_order"
      ),
      envir = environment()
    )
  }

  for (size_index in seq_along(training_sizes)) {
    n_per_group_value <- training_sizes[size_index]
    n_test_per_group <- min(
      opt$test_per_group,
      minimum_class_n - n_per_group_value
    )

    design_rows[[length(design_rows) + 1L]] <- data.table(
      panel_id = cfg$panel_id,
      panel_label = cfg$panel_label,
      n_per_group = as.integer(n_per_group_value),
      n_groups = nrow(group_counts),
      n_training_total = as.integer(n_per_group_value * nrow(group_counts)),
      n_test_per_group = as.integer(n_test_per_group),
      n_test_total = as.integer(n_test_per_group * nrow(group_counts)),
      n_replicates = as.integer(opt$n_replicates)
    )

    tasks <- lapply(seq_len(opt$n_replicates), function(replicate_value) {
      list(
        n_per_group = as.integer(n_per_group_value),
        n_test_per_group = as.integer(n_test_per_group),
        replicate = as.integer(replicate_value),
        seed = as.integer(
          opt$seed + panel_index * 1000000L +
            size_index * 10000L + replicate_value
        )
      )
    })

    message(
      "Training references per class: ", n_per_group_value,
      " | held out per class: ", n_test_per_group,
      " | replicates: ", opt$n_replicates
    )

    if (!is.null(cl)) {
      size_results <- parallel::parLapplyLB(
        cl,
        tasks,
        function(task) {
          process_one_replicate(
            task = task,
            X = X,
            snp_map = snp_map,
            ref_meta = ref_meta,
            query_baseline = query_baseline,
            cfg = cfg,
            settings = settings,
            baseK = baseK,
            extraK = extraK,
            marker_overlap_K = marker_overlap_K,
            full_snp_order = full_snp_order
          )
        }
      )
    } else {
      size_results <- lapply(tasks, function(task) {
        process_one_replicate(
          task = task,
          X = X,
          snp_map = snp_map,
          ref_meta = ref_meta,
          query_baseline = query_baseline,
          cfg = cfg,
          settings = settings,
          baseK = baseK,
          extraK = extraK,
          marker_overlap_K = marker_overlap_K,
          full_snp_order = full_snp_order
        )
      })
    }

    for (result_name in result_names) {
      all_results[[result_name]][[length(all_results[[result_name]]) + 1L]] <-
        rbindlist(
          lapply(size_results, `[[`, result_name),
          use.names = TRUE,
          fill = TRUE
        )
    }
  }

  if (!is.null(cl)) parallel::stopCluster(cl)
}

# =============================================================================
# Combine and summarise
# =============================================================================

combined <- lapply(all_results, function(parts) {
  rbindlist(parts, use.names = TRUE, fill = TRUE)
})

replicate_performance <- combined$replicate_performance
test_predictions <- combined$test_predictions
per_class_replicate <- combined$per_class
query_predictions <- combined$query_predictions
query_replicate <- combined$query_replicate
marker_stats <- combined$marker_stats
marker_overlap <- combined$marker_overlap
training_manifest <- combined$training_manifest
test_manifest <- combined$test_manifest
K_grid <- combined$K_grid

design <- rbindlist(design_rows, use.names = TRUE, fill = TRUE)
reference_counts <- rbindlist(reference_count_rows, use.names = TRUE, fill = TRUE)
baseline_queries <- rbindlist(baseline_query_rows, use.names = TRUE, fill = TRUE)
panel_index_table <- rbindlist(panel_index_rows, use.names = TRUE, fill = TRUE)

if (nrow(replicate_performance) == 0L) {
  stop("No downsampling replicate results were produced.")
}

performance_summary <- replicate_performance[, .(
  n_replicates_requested = as.integer(opt$n_replicates),
  n_replicates_completed = sum(status == "completed"),
  n_replicates_not_evaluable = sum(status != "completed"),
  n_training_total = as.integer(round(safe_mean(n_training))),
  n_test_per_replicate = as.integer(round(safe_mean(n_test))),
  mean_eligible_snps = safe_mean(n_eligible_snps),
  accuracy_mean = safe_mean(accuracy_all),
  accuracy_sd = safe_sd(accuracy_all),
  accuracy_q025 = safe_quantile(accuracy_all, 0.025),
  accuracy_median = safe_quantile(accuracy_all, 0.50),
  accuracy_q975 = safe_quantile(accuracy_all, 0.975),
  accepted_rate_mean = safe_mean(accepted_rate),
  accepted_accuracy_mean = safe_mean(accuracy_accepted),
  correct_accepted_rate_mean = safe_mean(correct_accepted_rate),
  uncertainty_rate_mean = safe_mean(uncertainty_rate),
  uncertainty_rate_q025 = safe_quantile(uncertainty_rate, 0.025),
  uncertainty_rate_q975 = safe_quantile(uncertainty_rate, 0.975),
  mean_maxK_posterior = safe_mean(mean_maxK_posterior),
  mean_maxK_gap = safe_mean(mean_maxK_gap),
  mean_global_agreement = safe_mean(mean_agreement),
  mean_tail_agreement = safe_mean(mean_tail_agreement)
), by = .(
  panel_id,
  panel_label,
  n_per_group
)][order(panel_id, n_per_group)]

pooled_performance <- test_predictions[, .(
  n_predictions = .N,
  correct_n = sum(correct_call, na.rm = TRUE),
  accuracy_all = mean(correct_call, na.rm = TRUE),
  accepted_n = sum(accepted, na.rm = TRUE),
  accepted_rate = mean(accepted, na.rm = TRUE),
  accuracy_accepted = if (sum(accepted, na.rm = TRUE) > 0L) {
    mean(correct_call[accepted %in% TRUE], na.rm = TRUE)
  } else {
    NA_real_
  },
  correct_accepted_rate = mean(correct_accepted, na.rm = TRUE),
  uncertainty_rate = mean(final_confidence == "Uncertain", na.rm = TRUE)
), by = .(
  panel_id,
  panel_label,
  n_per_group
)][order(panel_id, n_per_group)]

per_class_summary <- per_class_replicate[, .(
  n_replicates = .N,
  total_test_predictions = sum(n_test, na.rm = TRUE),
  accuracy_mean = safe_mean(accuracy_all),
  accuracy_sd = safe_sd(accuracy_all),
  accuracy_q025 = safe_quantile(accuracy_all, 0.025),
  accuracy_median = safe_quantile(accuracy_all, 0.50),
  accuracy_q975 = safe_quantile(accuracy_all, 0.975),
  accepted_rate_mean = safe_mean(accepted_rate),
  accepted_accuracy_mean = safe_mean(accuracy_accepted),
  uncertainty_rate_mean = safe_mean(uncertainty_rate)
), by = .(
  panel_id,
  panel_label,
  n_per_group,
  true_group
)][order(panel_id, n_per_group, true_group)]

query_stability <- query_predictions[, {
  modal <- mode_value(predicted_group)
  list(
    n_replicates = .N,
    baseline_call = baseline_call[1],
    baseline_confidence = baseline_confidence[1],
    modal_call = modal,
    modal_call_fraction = mean(predicted_group == modal, na.rm = TRUE),
    baseline_call_agreement = mean(call_matches_baseline, na.rm = TRUE),
    accepted_fraction = mean(accepted, na.rm = TRUE),
    high_fraction = mean(final_confidence == "High", na.rm = TRUE),
    moderate_fraction = mean(final_confidence == "Moderate", na.rm = TRUE),
    uncertainty_fraction = mean(final_confidence == "Uncertain", na.rm = TRUE),
    mean_maxK_posterior = safe_mean(maxK_top_posterior),
    minimum_maxK_posterior = if (any(is.finite(maxK_top_posterior))) {
      min(maxK_top_posterior, na.rm = TRUE)
    } else {
      NA_real_
    },
    mean_maxK_gap = safe_mean(maxK_gap),
    mean_global_agreement = safe_mean(agreement_frac),
    mean_tail_agreement = safe_mean(tail_agreement_frac)
  )
}, by = .(
  panel_id,
  panel_label,
  n_per_group,
  sample_id
)][order(panel_id, n_per_group, sample_id)]

query_overall_summary <- query_replicate[, .(
  n_replicates = .N,
  n_queries = as.integer(round(safe_mean(n_queries))),
  mean_fraction_calls_matching_baseline = safe_mean(
    fraction_calls_matching_baseline
  ),
  q025_fraction_calls_matching_baseline = safe_quantile(
    fraction_calls_matching_baseline,
    0.025
  ),
  q975_fraction_calls_matching_baseline = safe_quantile(
    fraction_calls_matching_baseline,
    0.975
  ),
  fraction_replicates_all_calls_match = mean(
    all_calls_match_baseline,
    na.rm = TRUE
  ),
  mean_accepted_fraction = safe_mean(accepted_fraction),
  mean_uncertainty_rate = safe_mean(uncertainty_rate)
), by = .(
  panel_id,
  panel_label,
  n_per_group
)][order(panel_id, n_per_group)]

marker_summary <- marker_stats[, .(
  n_replicates = .N,
  eligible_snps_mean = safe_mean(n_training_eligible_snps),
  eligible_snps_q025 = safe_quantile(n_training_eligible_snps, 0.025),
  eligible_snps_q975 = safe_quantile(n_training_eligible_snps, 0.975),
  common_fraction_mean = safe_mean(common_fraction_of_complete),
  rank_correlation_mean = safe_mean(spearman_rank_correlation_common),
  rank_correlation_q025 = safe_quantile(
    spearman_rank_correlation_common,
    0.025
  ),
  rank_correlation_q975 = safe_quantile(
    spearman_rank_correlation_common,
    0.975
  )
), by = .(
  panel_id,
  panel_label,
  n_per_group
)][order(panel_id, n_per_group)]

marker_overlap_summary <- marker_overlap[, .(
  n_replicates = .N,
  overlap_fraction_mean = safe_mean(overlap_fraction_of_complete),
  overlap_fraction_q025 = safe_quantile(overlap_fraction_of_complete, 0.025),
  overlap_fraction_q975 = safe_quantile(overlap_fraction_of_complete, 0.975),
  jaccard_mean = safe_mean(jaccard)
), by = .(
  panel_id,
  panel_label,
  n_per_group,
  K_requested
)][order(panel_id, n_per_group, K_requested)]

baseline_loo_selected <- baseline_loo[panel_id %in% panels_requested]

# =============================================================================
# Write tables
# =============================================================================

write_table_pair(panel_index_table, "downsampling_panel_index", opt$out_dir)
write_table_pair(reference_counts, "downsampling_reference_counts", opt$out_dir)
write_table_pair(design, "downsampling_design", opt$out_dir)
write_table_pair(baseline_loo_selected, "baseline_individual_loo_summary", opt$out_dir)
write_table_pair(baseline_queries, "baseline_query_assignments_by_panel", opt$out_dir)
write_table_pair(replicate_performance, "downsampling_replicate_performance", opt$out_dir)
write_table_pair(performance_summary, "downsampling_performance_summary", opt$out_dir)
write_table_pair(pooled_performance, "downsampling_pooled_performance", opt$out_dir)
write_table_pair(per_class_replicate, "downsampling_per_class_replicate", opt$out_dir)
write_table_pair(per_class_summary, "downsampling_per_class_summary", opt$out_dir)
write_table_pair(test_predictions, "downsampling_test_predictions", opt$out_dir, TRUE)
write_table_pair(query_predictions, "downsampling_query_predictions", opt$out_dir, TRUE)
write_table_pair(query_replicate, "downsampling_query_replicate", opt$out_dir)
write_table_pair(query_stability, "downsampling_query_stability", opt$out_dir)
write_table_pair(query_overall_summary, "downsampling_query_overall_summary", opt$out_dir)
write_table_pair(marker_stats, "downsampling_marker_stats", opt$out_dir)
write_table_pair(marker_summary, "downsampling_marker_summary", opt$out_dir)
write_table_pair(marker_overlap, "downsampling_marker_overlap", opt$out_dir, TRUE)
write_table_pair(marker_overlap_summary, "downsampling_marker_overlap_summary", opt$out_dir)
write_table_pair(training_manifest, "downsampling_training_manifest", opt$out_dir, TRUE)
write_table_pair(test_manifest, "downsampling_test_manifest", opt$out_dir, TRUE)
write_table_pair(K_grid, "downsampling_K_grid", opt$out_dir, TRUE)

params <- data.table(
  parameter = c(
    "n_replicates",
    "n_workers",
    "base_seed",
    "candidate_training_sizes",
    "test_per_group",
    "min_ref_per_group",
    "max_site_missing",
    "min_mac",
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
    "marker_overlap_K",
    "ranking_scope"
  ),
  value = c(
    as.character(opt$n_replicates),
    as.character(n_workers),
    as.character(opt$seed),
    paste(candidate_training_sizes, collapse = ","),
    as.character(opt$test_per_group),
    as.character(opt$min_ref_per_group),
    as.character(opt$max_site_missing),
    as.character(opt$min_mac),
    as.character(opt$pseudocount),
    as.character(opt$epsilon),
    as.character(opt$min_snps_query_macroregion),
    as.character(opt$min_snps_query_subregion),
    as.character(opt$high_posterior),
    as.character(opt$moderate_posterior),
    as.character(opt$min_gap_macroregion),
    as.character(opt$min_gap_subregion),
    paste(baseK, collapse = ","),
    paste(extraK, collapse = ","),
    as.character(auto_extend_K),
    as.character(include_all_snps_K),
    as.character(opt$max_K_points),
    as.character(opt$min_agreement),
    as.character(opt$min_K_available),
    as.character(opt$tail_fraction),
    as.character(opt$min_tail_agreement),
    paste(marker_overlap_K, collapse = ","),
    "training_only_reranking_of_corrected_step02_retained_snps"
  )
)

write_table_pair(params, "downsampling_params", opt$out_dir)

# =============================================================================
# Base R sensitivity figure
# =============================================================================

figure_path <- file.path(
  opt$out_dir,
  "FIG_reference_downsampling_sensitivity.pdf"
)

grDevices::pdf(figure_path, width = 8.27, height = 11.69)

for (panel_id_value in panels_requested) {
  perf <- performance_summary[panel_id == panel_id_value][order(n_per_group)]
  qry <- query_overall_summary[panel_id == panel_id_value][order(n_per_group)]

  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(3, 1), mar = c(4.5, 4.5, 3.5, 1.5))

  if (nrow(perf) > 0L) {
    plot(
      perf$n_per_group,
      perf$accuracy_mean,
      type = "b",
      log = "x",
      ylim = c(0, 1),
      xlab = "Training references per class",
      ylab = "Held out accuracy",
      main = paste0(perf$panel_label[1], ": balanced downsampling")
    )
    arrows(
      perf$n_per_group,
      perf$accuracy_q025,
      perf$n_per_group,
      perf$accuracy_q975,
      angle = 90,
      code = 3,
      length = 0.03
    )
    abline(h = seq(0, 1, 0.25), col = "grey85", lty = 3)

    plot(
      perf$n_per_group,
      perf$uncertainty_rate_mean,
      type = "b",
      log = "x",
      ylim = c(0, 1),
      xlab = "Training references per class",
      ylab = "Uncertainty rate",
      main = "Held out reference uncertainty"
    )
    arrows(
      perf$n_per_group,
      perf$uncertainty_rate_q025,
      perf$n_per_group,
      perf$uncertainty_rate_q975,
      angle = 90,
      code = 3,
      length = 0.03
    )
    abline(h = seq(0, 1, 0.25), col = "grey85", lty = 3)
  } else {
    plot.new()
    title("No held out performance results")
    plot.new()
    title("No uncertainty results")
  }

  if (nrow(qry) > 0L) {
    plot(
      qry$n_per_group,
      qry$mean_fraction_calls_matching_baseline,
      type = "b",
      log = "x",
      ylim = c(0, 1),
      xlab = "Training references per class",
      ylab = "Query call agreement",
      main = "Agreement with complete reference assignment"
    )
    arrows(
      qry$n_per_group,
      qry$q025_fraction_calls_matching_baseline,
      qry$n_per_group,
      qry$q975_fraction_calls_matching_baseline,
      angle = 90,
      code = 3,
      length = 0.03
    )
    abline(h = seq(0, 1, 0.25), col = "grey85", lty = 3)
  } else {
    plot.new()
    title("No query stability results")
  }

  par(old_par)
}

dev.off()

# =============================================================================
# Plain text summary and run information
# =============================================================================

summary_lines <- c(
  "PESTFLY BALANCED REFERENCE DOWNSAMPLING SENSITIVITY",
  "",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Repository root: ", repo_root),
  paste0("Replicates per training size: ", opt$n_replicates),
  paste0("Parallel workers: ", n_workers),
  "",
  "METHOD",
  paste0(
    "Equal numbers of training references were sampled per class. ",
    "A balanced independent test set was sampled from the remaining references."
  ),
  paste0(
    "SNP missingness, minor allele count, corrected Hudson FST ranking, and ",
    "allele frequencies were recalculated from training references only."
  ),
  paste0(
    "The analysis reranked the corrected Step 02 retained SNP set and did not ",
    "rediscover alternative SNPs from the original alignments."
  ),
  "",
  "PANEL RESULTS"
)

for (panel_id_value in panels_requested) {
  perf <- performance_summary[panel_id == panel_id_value][order(n_per_group)]
  qry <- query_overall_summary[panel_id == panel_id_value][order(n_per_group)]
  counts <- reference_counts[panel_id == panel_id_value]

  summary_lines <- c(
    summary_lines,
    "",
    if (nrow(perf) > 0L) perf$panel_label[1] else panel_id_value,
    paste0(
      "Reference counts: ",
      paste(
        paste0(counts$analytical_group, "=", counts$N),
        collapse = ", "
      )
    ),
    paste0(
      "Balanced training sizes: ",
      paste(perf$n_per_group, collapse = ", ")
    )
  )

  if (nrow(perf) > 0L) {
    first_perf <- perf[1]
    last_perf <- perf[nrow(perf)]
    last_qry <- qry[n_per_group == last_perf$n_per_group]

    summary_lines <- c(
      summary_lines,
      paste0(
        "Smallest size accuracy mean: ",
        format(round(first_perf$accuracy_mean, 4), nsmall = 4),
        " | uncertainty mean: ",
        format(round(first_perf$uncertainty_rate_mean, 4), nsmall = 4)
      ),
      paste0(
        "Largest size accuracy mean: ",
        format(round(last_perf$accuracy_mean, 4), nsmall = 4),
        " | uncertainty mean: ",
        format(round(last_perf$uncertainty_rate_mean, 4), nsmall = 4)
      )
    )

    if (nrow(last_qry) > 0L) {
      summary_lines <- c(
        summary_lines,
        paste0(
          "Largest size query call agreement with complete reference analysis: ",
          format(
            round(last_qry$mean_fraction_calls_matching_baseline, 4),
            nsmall = 4
          )
        )
      )
    }
  }
}

summary_lines <- c(
  summary_lines,
  "",
  "INTERPRETATION",
  paste0(
    "Use the performance, per class, query stability, and marker stability ",
    "tables for interpretation."
  ),
  paste0(
    "The individual LOO results from the complete reference set are supplied ",
    "for comparison but are not treated as downsampling replicates."
  ),
  "",
  "SUMMARY EXPORTS",
  "Zip the complete return_bundle directory."
)

summary_path <- file.path(opt$out_dir, "REFERENCE_DOWNSAMPLING_SUMMARY.txt")
writeLines(summary_lines, summary_path, useBytes = TRUE)

run_info <- list(
  analysis = "balanced_reference_downsampling_sensitivity",
  timestamp = Sys.time(),
  repo_root = repo_root,
  inputs = list(
    metadata = metadata_path,
    step2_dir = opt$step2_dir,
    assignment = assignment_path,
    baseline_loo = baseline_loo_path
  ),
  panels = panel_index_table,
  design = design,
  params = params,
  hudson_estimator = list(
    numerator = paste0(
      "(p1-p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1)"
    ),
    denominator = "p1(1-p2) + p2(1-p1)",
    negative_estimates_retained = TRUE,
    multiclass_score = "maximum pairwise Hudson FST"
  ),
  randomisation = list(
    base_seed = opt$seed,
    seed_rule = paste0(
      "base seed + panel index * 1000000 + training size index * 10000 + replicate"
    )
  ),
  session_info = utils::sessionInfo()
)

saveRDS(
  run_info,
  file.path(opt$out_dir, "reference_downsampling_run_info.rds")
)

# =============================================================================
# Compact return bundle
# =============================================================================

return_files <- c(
  "REFERENCE_DOWNSAMPLING_SUMMARY.txt",
  "FIG_reference_downsampling_sensitivity.pdf",
  "downsampling_panel_index.tsv",
  "downsampling_reference_counts.tsv",
  "downsampling_design.tsv",
  "baseline_individual_loo_summary.tsv",
  "baseline_query_assignments_by_panel.tsv",
  "downsampling_replicate_performance.tsv",
  "downsampling_performance_summary.tsv",
  "downsampling_pooled_performance.tsv",
  "downsampling_per_class_replicate.tsv",
  "downsampling_per_class_summary.tsv",
  "downsampling_test_predictions.tsv.gz",
  "downsampling_query_predictions.tsv.gz",
  "downsampling_query_replicate.tsv",
  "downsampling_query_stability.tsv",
  "downsampling_query_overall_summary.tsv",
  "downsampling_marker_stats.tsv",
  "downsampling_marker_summary.tsv",
  "downsampling_marker_overlap.tsv.gz",
  "downsampling_marker_overlap_summary.tsv",
  "downsampling_training_manifest.tsv.gz",
  "downsampling_test_manifest.tsv.gz",
  "downsampling_K_grid.tsv.gz",
  "downsampling_params.tsv",
  "reference_downsampling_run_info.rds"
)

for (filename in return_files) {
  source_path <- file.path(opt$out_dir, filename)
  if (!file.exists(source_path)) {
    stop("Expected output was not created: ", source_path)
  }

  copied <- file.copy(
    source_path,
    file.path(return_dir, filename),
    overwrite = TRUE,
    copy.date = TRUE
  )

  if (!isTRUE(copied)) stop("Could not copy into return bundle: ", filename)
}

message("\nBalanced reference downsampling sensitivity completed.")
message("Summary exports: ", return_dir)
