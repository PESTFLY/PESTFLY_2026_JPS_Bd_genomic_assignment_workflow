#!/usr/bin/env Rscript

# PESTFLY: Alternative marker set resampling
#
# Partition the top 2K Hudson ranked markers into rank matched A and B sets of K markers,
# giving thirty panels from fifteen partitions.
# A/B sets are disjoint within a partition; different partitions may overlap.
#
# Run from the repository root:
#   Rscript steps/validation/06_marker_resampling/run.R
# Inputs, outputs and options: adjacent README.md. Full methods: associated manuscript.

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

  y <- tolower(clean_text(x))
  out <- rep(NA, length(y))

  out[y %in% c("true", "t", "1", "yes", "y", "reference", "ref")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n", "query", "intercept")] <- FALSE

  out
}

parse_integer_vector <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) {
    return(integer())
  }

  values <- suppressWarnings(as.integer(clean_text(strsplit(x, ",")[[1]])))
  values <- values[!is.na(values) & values > 0L]
  sort(unique(values))
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

mode_value <- function(x) {
  x <- clean_text(x)
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA_character_)

  counts <- sort(table(x), decreasing = TRUE)
  candidates <- names(counts)[counts == max(counts)]
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
    default = "results/validation/06_marker_resampling",
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
    "--panel_sizes",
    type = "character",
    default = "500,1000,2000,5000",
    help = "Comma separated fixed marker panel sizes [default %default]"
  ),
  make_option(
    "--n_partitions",
    type = "integer",
    default = 15L,
    help = "Complementary A/B partitions per panel and K [default %default]"
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
    default = 20261001L,
    help = "Base random seed [default %default]"
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
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$qc_dir <- resolve_path(opt$qc_dir)
opt$step2_dir <- resolve_path(opt$step2_dir)
opt$assignment_dir <- resolve_path(opt$assignment_dir)
opt$loo_dir <- resolve_path(opt$loo_dir)
opt$out_dir <- resolve_path(opt$out_dir)
opt$n_partitions <- max(1L, as.integer(opt$n_partitions))
opt$seed <- as.integer(opt$seed)

panels_requested <- clean_text(strsplit(opt$panels, ",")[[1]])
panels_requested <- panels_requested[!is.na(panels_requested)]
panel_sizes_requested <- parse_integer_vector(opt$panel_sizes)

if (length(panels_requested) == 0L) stop("No panels were requested.")
if (length(panel_sizes_requested) == 0L) stop("No marker panel sizes were requested.")

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
# Likelihood and fixed panel prediction helpers
# =============================================================================

conf_label <- function(pmax, high = 0.95, moderate = 0.85) {
  if (is.na(pmax)) return("NA")
  if (pmax >= high) return("High")
  if (pmax >= moderate) return("Moderate")
  "Low"
}

softmax_one <- function(log_likelihoods) {
  if (all(is.na(log_likelihoods))) {
    return(rep(NA_real_, length(log_likelihoods)))
  }

  maximum <- max(log_likelihoods, na.rm = TRUE)
  weights <- exp(log_likelihoods - maximum)
  weights[is.na(weights)] <- 0

  if (sum(weights) <= 0) {
    return(rep(NA_real_, length(log_likelihoods)))
  }

  weights / sum(weights)
}

frequency_from_counts <- function(alt_counts, observed_counts, pseudocount) {
  frequency <- matrix(
    NA_real_,
    nrow = nrow(alt_counts),
    ncol = ncol(alt_counts),
    dimnames = dimnames(alt_counts)
  )

  available <- observed_counts > 0L
  frequency[available] <-
    (alt_counts[available] + pseudocount) /
    (observed_counts[available] + 2 * pseudocount)

  frequency
}

loglik_matrix <- function(X_samples, frequency, epsilon) {
  sample_ids <- colnames(X_samples)
  groups <- colnames(frequency)

  log_likelihood <- matrix(
    NA_real_,
    nrow = ncol(X_samples),
    ncol = ncol(frequency),
    dimnames = list(sample_ids, groups)
  )

  group_observed <- matrix(
    0L,
    nrow = ncol(X_samples),
    ncol = ncol(frequency),
    dimnames = list(sample_ids, groups)
  )

  query_is_one <- X_samples == 1L
  query_is_zero <- X_samples == 0L
  query_is_one[is.na(query_is_one)] <- FALSE
  query_is_zero[is.na(query_is_zero)] <- FALSE

  for (group_index in seq_along(groups)) {
    p <- frequency[, group_index]
    p1 <- p * (1 - epsilon) + (1 - p) * epsilon
    p0 <- (1 - p) * (1 - epsilon) + p * epsilon

    p1 <- pmin(pmax(p1, 1e-12), 1 - 1e-12)
    p0 <- pmin(pmax(p0, 1e-12), 1 - 1e-12)

    log_p1 <- log(p1)
    log_p0 <- log(p0)
    frequency_available <- is.finite(p)

    log_p1[!is.finite(log_p1)] <- 0
    log_p0[!is.finite(log_p0)] <- 0

    log_likelihood[, group_index] <- as.numeric(
      crossprod(log_p1, query_is_one) +
        crossprod(log_p0, query_is_zero)
    )

    group_observed[, group_index] <- as.integer(crossprod(
      as.integer(frequency_available),
      !is.na(X_samples)
    ))
  }

  log_likelihood[group_observed == 0L] <- NA_real_

  list(
    log_likelihood = log_likelihood,
    group_observed = group_observed
  )
}

summarise_prediction_matrix <- function(
    log_likelihood,
    n_snps_used,
    sample_ids,
    cfg,
    settings,
    true_groups = NULL
) {
  groups <- colnames(log_likelihood)
  rows <- vector("list", length(sample_ids))

  if (is.null(true_groups)) {
    true_groups <- rep(NA_character_, length(sample_ids))
  }

  for (i in seq_along(sample_ids)) {
    sid <- sample_ids[i]
    true_group <- clean_text(true_groups[i])
    n_used <- as.integer(n_snps_used[i])

    if (is.na(n_used) || n_used < cfg$min_snps_query) {
      rows[[i]] <- data.table(
        sample_id = sid,
        true_group = true_group,
        predicted_group = NA_character_,
        top_posterior = NA_real_,
        second_group = NA_character_,
        second_posterior = NA_real_,
        posterior_gap = NA_real_,
        n_snps_used = n_used,
        fixed_panel_confidence = "Uncertain",
        fixed_panel_reason = "too_few_snps",
        accepted = FALSE,
        correct = if (!is.na(true_group)) FALSE else NA
      )
      next
    }

    posterior <- softmax_one(log_likelihood[i, ])
    names(posterior) <- groups

    if (all(is.na(posterior))) {
      rows[[i]] <- data.table(
        sample_id = sid,
        true_group = true_group,
        predicted_group = NA_character_,
        top_posterior = NA_real_,
        second_group = NA_character_,
        second_posterior = NA_real_,
        posterior_gap = NA_real_,
        n_snps_used = n_used,
        fixed_panel_confidence = "Uncertain",
        fixed_panel_reason = "no_likelihood",
        accepted = FALSE,
        correct = if (!is.na(true_group)) FALSE else NA
      )
      next
    }

    order_index <- order(posterior, decreasing = TRUE)
    top_index <- order_index[1]
    second_index <- if (length(order_index) >= 2L) order_index[2] else NA_integer_
    predicted_group <- groups[top_index]
    top_posterior <- unname(posterior[top_index])
    second_group <- if (!is.na(second_index)) groups[second_index] else NA_character_
    second_posterior <- if (!is.na(second_index)) {
      unname(posterior[second_index])
    } else {
      NA_real_
    }
    posterior_gap <- if (!is.na(second_posterior)) {
      top_posterior - second_posterior
    } else {
      NA_real_
    }

    posterior_label <- conf_label(
      top_posterior,
      high = settings$high_posterior,
      moderate = settings$moderate_posterior
    )

    if (!(posterior_label %in% c("High", "Moderate"))) {
      fixed_confidence <- "Uncertain"
      fixed_reason <- "low_posterior"
    } else if (!is.na(posterior_gap) && posterior_gap < cfg$min_gap) {
      fixed_confidence <- "Uncertain"
      fixed_reason <- "low_posterior_gap"
    } else {
      fixed_confidence <- posterior_label
      fixed_reason <- "ok"
    }

    accepted <- fixed_confidence %in% c("High", "Moderate") &&
      identical(fixed_reason, "ok")

    rows[[i]] <- data.table(
      sample_id = sid,
      true_group = true_group,
      predicted_group = predicted_group,
      top_posterior = top_posterior,
      second_group = second_group,
      second_posterior = second_posterior,
      posterior_gap = posterior_gap,
      n_snps_used = n_used,
      fixed_panel_confidence = fixed_confidence,
      fixed_panel_reason = fixed_reason,
      accepted = accepted,
      correct = if (!is.na(true_group)) predicted_group == true_group else NA
    )
  }

  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

predict_reference_loo <- function(
    X_reference,
    reference_groups,
    alt_counts,
    observed_counts,
    cfg,
    settings
) {
  reference_ids <- colnames(X_reference)
  groups <- colnames(alt_counts)
  full_frequency <- frequency_from_counts(
    alt_counts,
    observed_counts,
    settings$pseudocount
  )

  likelihood <- loglik_matrix(
    X_reference,
    full_frequency,
    settings$epsilon
  )$log_likelihood

  for (group_index in seq_along(groups)) {
    group_value <- groups[group_index]
    in_group <- which(reference_groups == group_value)
    if (length(in_group) == 0L) next

    X_group <- X_reference[, in_group, drop = FALSE]
    query_observed <- !is.na(X_group)
    query_is_one <- X_group == 1L
    query_is_zero <- X_group == 0L
    query_is_one[is.na(query_is_one)] <- FALSE
    query_is_zero[is.na(query_is_zero)] <- FALSE

    alt_excluding <- matrix(
      alt_counts[, group_index],
      nrow = nrow(X_group),
      ncol = ncol(X_group)
    ) - query_is_one

    observed_excluding <- matrix(
      observed_counts[, group_index],
      nrow = nrow(X_group),
      ncol = ncol(X_group)
    ) - query_observed

    frequency_excluding <- matrix(
      NA_real_,
      nrow = nrow(X_group),
      ncol = ncol(X_group)
    )

    available <- observed_excluding > 0L
    frequency_excluding[available] <-
      (alt_excluding[available] + settings$pseudocount) /
      (observed_excluding[available] + 2 * settings$pseudocount)

    p1 <- frequency_excluding * (1 - settings$epsilon) +
      (1 - frequency_excluding) * settings$epsilon
    p0 <- (1 - frequency_excluding) * (1 - settings$epsilon) +
      frequency_excluding * settings$epsilon

    p1 <- pmin(pmax(p1, 1e-12), 1 - 1e-12)
    p0 <- pmin(pmax(p0, 1e-12), 1 - 1e-12)

    log_p1 <- log(p1)
    log_p0 <- log(p0)
    log_p1[!is.finite(log_p1)] <- 0
    log_p0[!is.finite(log_p0)] <- 0

    likelihood[in_group, group_index] <- colSums(
      query_is_one * log_p1 + query_is_zero * log_p0
    )

    no_group_data <- colSums(query_observed & available) == 0L
    likelihood[in_group[no_group_data], group_index] <- NA_real_
  }

  total_observed <- rowSums(observed_counts)
  n_used <- integer(length(reference_ids))

  for (i in seq_along(reference_ids)) {
    query_observed <- !is.na(X_reference[, i])
    n_used[i] <- sum(query_observed & (total_observed - query_observed) > 0L)
  }

  summarise_prediction_matrix(
    log_likelihood = likelihood,
    n_snps_used = n_used,
    sample_ids = reference_ids,
    cfg = cfg,
    settings = settings,
    true_groups = reference_groups
  )
}

predict_queries <- function(
    X_query,
    alt_counts,
    observed_counts,
    cfg,
    settings
) {
  if (ncol(X_query) == 0L) return(data.table())

  full_frequency <- frequency_from_counts(
    alt_counts,
    observed_counts,
    settings$pseudocount
  )

  likelihood_result <- loglik_matrix(
    X_query,
    full_frequency,
    settings$epsilon
  )

  any_frequency <- rowSums(!is.na(full_frequency)) > 0L
  n_used <- colSums(!is.na(X_query) & any_frequency)

  summarise_prediction_matrix(
    log_likelihood = likelihood_result$log_likelihood,
    n_snps_used = n_used,
    sample_ids = colnames(X_query),
    cfg = cfg,
    settings = settings,
    true_groups = NULL
  )
}

# =============================================================================
# Marker partition and evaluation helpers
# =============================================================================

make_rank_matched_partition <- function(K, seed) {
  set.seed(seed)
  pairs <- matrix(seq_len(2L * K), ncol = 2L, byrow = TRUE)
  choose_first <- sample(c(TRUE, FALSE), K, replace = TRUE)

  set_a <- ifelse(choose_first, pairs[, 1], pairs[, 2])
  set_b <- ifelse(choose_first, pairs[, 2], pairs[, 1])

  list(
    A = sort(as.integer(set_a)),
    B = sort(as.integer(set_b))
  )
}

evaluate_marker_set <- function(
    selected_index,
    K,
    set_id,
    set_type,
    partition,
    complement_side,
    seed,
    X,
    snp_map,
    reference_ids,
    reference_groups,
    query_ids,
    query_baseline,
    alt_all,
    observed_all,
    cfg,
    settings
) {
  selected_index <- sort(as.integer(selected_index))
  selected_ids <- snp_map$snp_id[selected_index]

  X_reference <- X[selected_ids, reference_ids, drop = FALSE]
  alt_counts <- alt_all[selected_index, , drop = FALSE]
  observed_counts <- observed_all[selected_index, , drop = FALSE]

  loo_predictions <- predict_reference_loo(
    X_reference = X_reference,
    reference_groups = reference_groups,
    alt_counts = alt_counts,
    observed_counts = observed_counts,
    cfg = cfg,
    settings = settings
  )

  query_predictions <- if (length(query_ids) > 0L) {
    predict_queries(
      X_query = X[selected_ids, query_ids, drop = FALSE],
      alt_counts = alt_counts,
      observed_counts = observed_counts,
      cfg = cfg,
      settings = settings
    )
  } else {
    data.table()
  }

  annotation <- list(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label,
    K = as.integer(K),
    set_id = set_id,
    set_type = set_type,
    partition = as.integer(partition),
    complement_side = complement_side,
    seed = as.integer(seed)
  )

  for (name in names(annotation)) {
    set(loo_predictions, j = name, value = annotation[[name]])
    if (nrow(query_predictions) > 0L) {
      set(query_predictions, j = name, value = annotation[[name]])
    }
  }

  if (nrow(query_predictions) > 0L) {
    query_predictions <- merge(
      query_predictions,
      query_baseline,
      by = "sample_id",
      all.x = TRUE,
      sort = FALSE
    )

    query_predictions[, matches_primary_multiK :=
      !is.na(predicted_group) & predicted_group == primary_baseline_call]
  }

  loo_performance <- loo_predictions[, .(
    n_references = .N,
    n_attempted = sum(!is.na(predicted_group)),
    n_accepted = sum(accepted),
    accuracy_all = mean(correct %in% TRUE),
    accuracy_attempted = if (sum(!is.na(predicted_group)) > 0L) {
      mean(correct[!is.na(predicted_group)] %in% TRUE)
    } else {
      NA_real_
    },
    accepted_rate = mean(accepted),
    accepted_accuracy = if (sum(accepted) > 0L) {
      mean(correct[accepted] %in% TRUE)
    } else {
      NA_real_
    },
    correct_accepted_rate = mean(correct %in% TRUE & accepted),
    uncertainty_rate = mean(!accepted),
    mean_top_posterior = safe_mean(top_posterior),
    mean_posterior_gap = safe_mean(posterior_gap),
    mean_snps_used = safe_mean(n_snps_used)
  )]

  for (name in names(annotation)) {
    set(loo_performance, j = name, value = annotation[[name]])
  }

  loo_per_class <- loo_predictions[, .(
    n_references = .N,
    n_attempted = sum(!is.na(predicted_group)),
    n_accepted = sum(accepted),
    accuracy_all = mean(correct %in% TRUE),
    accepted_rate = mean(accepted),
    accepted_accuracy = if (sum(accepted) > 0L) {
      mean(correct[accepted] %in% TRUE)
    } else {
      NA_real_
    },
    uncertainty_rate = mean(!accepted)
  ), by = .(analytical_group = true_group)]

  for (name in names(annotation)) {
    set(loo_per_class, j = name, value = annotation[[name]])
  }

  query_performance <- if (nrow(query_predictions) > 0L) {
    query_predictions[, .(
      n_queries = .N,
      n_attempted = sum(!is.na(predicted_group)),
      n_accepted = sum(accepted),
      primary_call_agreement = mean(matches_primary_multiK %in% TRUE),
      all_primary_calls_match = all(matches_primary_multiK %in% TRUE),
      accepted_rate = mean(accepted),
      uncertainty_rate = mean(!accepted),
      mean_top_posterior = safe_mean(top_posterior),
      mean_posterior_gap = safe_mean(posterior_gap),
      mean_snps_used = safe_mean(n_snps_used)
    )]
  } else {
    data.table()
  }

  if (nrow(query_performance) > 0L) {
    for (name in names(annotation)) {
      set(query_performance, j = name, value = annotation[[name]])
    }
  }

  marker_manifest <- data.table(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label,
    K = as.integer(K),
    set_id = set_id,
    set_type = set_type,
    partition = as.integer(partition),
    complement_side = complement_side,
    seed = as.integer(seed),
    snp_id = selected_ids,
    orthogroup = clean_text(snp_map$og[selected_index]),
    corrected_hudson_rank = as.integer(selected_index),
    corrected_hudson_score = as.numeric(snp_map$score[selected_index]),
    in_conventional_topK = selected_index <= K
  )

  marker_summary <- marker_manifest[, .(
    n_markers = .N,
    n_unique_snps = uniqueN(snp_id),
    n_unique_orthogroups = uniqueN(orthogroup),
    n_from_conventional_topK = sum(in_conventional_topK),
    fraction_from_conventional_topK = mean(in_conventional_topK),
    minimum_rank = min(corrected_hudson_rank),
    median_rank = stats::median(corrected_hudson_rank),
    maximum_rank = max(corrected_hudson_rank),
    mean_hudson_score = safe_mean(corrected_hudson_score)
  ), by = .(
    panel_id,
    panel_label,
    K,
    set_id,
    set_type,
    partition,
    complement_side,
    seed
  )]

  list(
    loo_predictions = loo_predictions,
    loo_performance = loo_performance,
    loo_per_class = loo_per_class,
    query_predictions = query_predictions,
    query_performance = query_performance,
    marker_manifest = marker_manifest,
    marker_summary = marker_summary
  )
}

evaluate_partition <- function(
    task,
    X,
    snp_map,
    reference_ids,
    reference_groups,
    query_ids,
    query_baseline,
    alt_all,
    observed_all,
    cfg,
    settings
) {
  partition_sets <- make_rank_matched_partition(task$K, task$seed)
  outputs <- vector("list", 2L)
  sides <- c("A", "B")

  for (i in seq_along(sides)) {
    side <- sides[i]
    set_id <- paste0("K", task$K, "_P", task$partition, "_", side)

    outputs[[i]] <- evaluate_marker_set(
      selected_index = partition_sets[[side]],
      K = task$K,
      set_id = set_id,
      set_type = "rank_matched_alternative",
      partition = task$partition,
      complement_side = side,
      seed = task$seed,
      X = X,
      snp_map = snp_map,
      reference_ids = reference_ids,
      reference_groups = reference_groups,
      query_ids = query_ids,
      query_baseline = query_baseline,
      alt_all = alt_all,
      observed_all = observed_all,
      cfg = cfg,
      settings = settings
    )
  }

  combined <- lapply(names(outputs[[1]]), function(result_name) {
    rbindlist(
      lapply(outputs, `[[`, result_name),
      use.names = TRUE,
      fill = TRUE
    )
  })
  names(combined) <- names(outputs[[1]])

  set_a <- partition_sets$A
  set_b <- partition_sets$B

  combined$complement_check <- data.table(
    panel_id = cfg$panel_id,
    panel_label = cfg$panel_label,
    K = as.integer(task$K),
    partition = as.integer(task$partition),
    seed = as.integer(task$seed),
    set_A_n = length(set_a),
    set_B_n = length(set_b),
    intersection_n = length(intersect(set_a, set_b)),
    union_n = length(union(set_a, set_b)),
    expected_union_n = as.integer(2L * task$K),
    passed = length(set_a) == task$K &&
      length(set_b) == task$K &&
      length(intersect(set_a, set_b)) == 0L &&
      length(union(set_a, set_b)) == 2L * task$K
  )

  combined
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

metadata <- fread(metadata_path)
setnames(metadata, tolower(names(metadata)))

required_metadata <- c("sample_id", "is_reference", "macroregion_3", "subregion")
missing_metadata <- setdiff(required_metadata, names(metadata))

if (length(missing_metadata) > 0L) {
  stop("metadata_clean.tsv is missing: ", paste(missing_metadata, collapse = ", "))
}

for (name in names(metadata)) {
  if (!inherits(metadata[[name]], c("numeric", "integer", "logical", "Date", "POSIXct"))) {
    metadata[[name]] <- clean_text(metadata[[name]])
  }
}

metadata[, `:=`(
  sample_id = clean_text(sample_id),
  is_reference = as_clean_logical(is_reference),
  macroregion_3 = clean_text(macroregion_3),
  subregion = clean_text(subregion)
)]

if (anyDuplicated(metadata$sample_id)) {
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
  pseudocount = as.numeric(opt$pseudocount),
  epsilon = as.numeric(opt$epsilon),
  high_posterior = as.numeric(opt$high_posterior),
  moderate_posterior = as.numeric(opt$moderate_posterior)
)

message("Repository root: ", repo_root)
message("Corrected Step 02: ", opt$step2_dir)
message("Output directory: ", opt$out_dir)
message("Complementary partitions per K: ", opt$n_partitions)
message("Alternative panels per K: ", 2L * opt$n_partitions)
message("Parallel workers: ", n_workers)

# =============================================================================
# Run panels
# =============================================================================

result_names <- c(
  "loo_predictions",
  "loo_performance",
  "loo_per_class",
  "query_predictions",
  "query_performance",
  "marker_manifest",
  "marker_summary",
  "complement_check"
)

alternative_results <- setNames(lapply(result_names, function(x) list()), result_names)
baseline_results <- setNames(lapply(result_names[1:7], function(x) list()), result_names[1:7])
panel_index_rows <- list()
design_rows <- list()
reference_count_rows <- list()
baseline_query_rows <- list()

worker_functions <- c(
  "clean_text",
  "safe_mean",
  "conf_label",
  "softmax_one",
  "frequency_from_counts",
  "loglik_matrix",
  "summarise_prediction_matrix",
  "predict_reference_loo",
  "predict_queries",
  "make_rank_matched_partition",
  "evaluate_marker_set",
  "evaluate_partition"
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
  required_map <- c("snp_id", "og", "score")
  missing_map <- setdiff(required_map, names(snp_map))

  if (length(missing_map) > 0L) {
    stop("snp_map.tsv is missing: ", paste(missing_map, collapse = ", "))
  }

  if ("global_rank" %in% names(snp_map)) {
    setorder(snp_map, global_rank)
  } else {
    setorder(snp_map, -score)
  }

  snp_map <- snp_map[snp_id %in% rownames(X)]

  if (nrow(snp_map) == 0L) stop("No mapped SNPs remain for ", cfg$panel_id)
  if (anyDuplicated(snp_map$snp_id)) stop("Duplicated snp_id values in ", cfg$panel_id)
  if (anyDuplicated(snp_map$og)) {
    stop(
      "The corrected panel contains more than one retained SNP per orthogroup: ",
      cfg$panel_id
    )
  }

  X <- X[snp_map$snp_id, , drop = FALSE]

  valid_sizes <- panel_sizes_requested[2L * panel_sizes_requested <= nrow(snp_map)]
  omitted_sizes <- setdiff(panel_sizes_requested, valid_sizes)

  if (length(omitted_sizes) > 0L) {
    warning(
      "Omitting K values without at least 2K ranked markers for ",
      cfg$panel_id,
      ": ",
      paste(omitted_sizes, collapse = ", ")
    )
  }

  if (length(valid_sizes) == 0L) {
    stop("No requested marker panel size is evaluable for ", cfg$panel_id)
  }

  reference_metadata <- metadata[is_reference %in% TRUE]

  if (!is.null(cfg$allowed_groups)) {
    reference_metadata <- reference_metadata[
      get(cfg$group_col) %in% cfg$allowed_groups
    ]
  }

  if (!is.na(cfg$macroregion_filter)) {
    reference_metadata <- reference_metadata[
      macroregion_3 == cfg$macroregion_filter
    ]
  }

  reference_metadata <- reference_metadata[!is.na(get(cfg$group_col))]
  reference_metadata <- reference_metadata[sample_id %in% colnames(X)]
  setorder(reference_metadata, sample_id)

  group_counts <- reference_metadata[, .N, by = c(cfg$group_col)]
  setnames(group_counts, cfg$group_col, "analytical_group")
  setorder(group_counts, analytical_group)

  if (nrow(group_counts) < 2L) {
    stop("Fewer than two reference classes for ", cfg$panel_id)
  }

  reference_ids <- reference_metadata$sample_id
  reference_groups <- clean_text(reference_metadata[[cfg$group_col]])
  groups <- group_counts$analytical_group

  if (identical(cfg$query_panel, "P1_macroregion")) {
    query_baseline <- baseline_assignment[, .(
      sample_id,
      primary_baseline_call = clean_text(macroregion_call),
      primary_baseline_confidence = clean_text(macroregion_confidence),
      primary_baseline_reason = clean_text(macroregion_reason)
    )]
  } else {
    query_baseline <- baseline_assignment[
      subregion_panel == cfg$query_panel,
      .(
        sample_id,
        primary_baseline_call = clean_text(subregion_call),
        primary_baseline_confidence = clean_text(subregion_confidence),
        primary_baseline_reason = clean_text(subregion_reason)
      )
    ]
  }

  query_baseline <- query_baseline[sample_id %in% colnames(X)]
  setorder(query_baseline, sample_id)
  query_ids <- query_baseline$sample_id

  if (nrow(query_baseline) == 0L) {
    stop("No applicable baseline queries found for ", cfg$panel_id)
  }

  X_reference_all <- X[, reference_ids, drop = FALSE]
  alt_all <- matrix(
    0L,
    nrow = nrow(X),
    ncol = length(groups),
    dimnames = list(rownames(X), groups)
  )
  observed_all <- alt_all

  for (group_index in seq_along(groups)) {
    ids <- reference_ids[reference_groups == groups[group_index]]
    X_group <- X[, ids, drop = FALSE]
    alt_all[, group_index] <- rowSums(X_group == 1L, na.rm = TRUE)
    observed_all[, group_index] <- rowSums(!is.na(X_group))
  }

  reference_count_rows[[panel_index]] <- copy(group_counts)[, `:=`(
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
    n_reference_samples = nrow(reference_metadata),
    n_complete_panel_snps = nrow(snp_map),
    one_snp_per_orthogroup = !anyDuplicated(snp_map$og),
    marker_panel_sizes = paste(valid_sizes, collapse = ","),
    n_partitions_per_K = as.integer(opt$n_partitions),
    n_alternative_sets_per_K = as.integer(2L * opt$n_partitions),
    n_applicable_queries = nrow(query_baseline)
  )

  for (K in valid_sizes) {
    design_rows[[length(design_rows) + 1L]] <- data.table(
      panel_id = cfg$panel_id,
      panel_label = cfg$panel_label,
      K = as.integer(K),
      candidate_pool_size = as.integer(2L * K),
      n_partitions = as.integer(opt$n_partitions),
      n_alternative_sets = as.integer(2L * opt$n_partitions),
      expected_markers_per_set = as.integer(K),
      expected_complement_overlap = 0L,
      expected_fraction_from_conventional_topK = 0.5
    )

    message(
      "K: ", K,
      " | candidate pool: top ", 2L * K,
      " | alternative sets: ", 2L * opt$n_partitions
    )

    baseline_result <- evaluate_marker_set(
      selected_index = seq_len(K),
      K = K,
      set_id = paste0("K", K, "_baseline_topK"),
      set_type = "conventional_topK",
      partition = 0L,
      complement_side = "baseline",
      seed = NA_integer_,
      X = X,
      snp_map = snp_map,
      reference_ids = reference_ids,
      reference_groups = reference_groups,
      query_ids = query_ids,
      query_baseline = query_baseline,
      alt_all = alt_all,
      observed_all = observed_all,
      cfg = cfg,
      settings = settings
    )

    for (result_name in names(baseline_results)) {
      baseline_results[[result_name]][[length(baseline_results[[result_name]]) + 1L]] <-
        baseline_result[[result_name]]
    }
  }

  tasks <- list()
  for (size_index in seq_along(valid_sizes)) {
    K <- valid_sizes[size_index]

    for (partition in seq_len(opt$n_partitions)) {
      tasks[[length(tasks) + 1L]] <- list(
        K = as.integer(K),
        partition = as.integer(partition),
        seed = as.integer(
          opt$seed + panel_index * 1000000L +
            size_index * 10000L + partition
        )
      )
    }
  }

  panel_workers <- min(n_workers, length(tasks))
  cluster <- NULL

  if (panel_workers > 1L) {
    cluster <- parallel::makeCluster(panel_workers)
    parallel::clusterEvalQ(cluster, {
      suppressPackageStartupMessages(library(data.table))
      NULL
    })

    parallel::clusterExport(
      cluster,
      varlist = c(
        worker_functions,
        "X",
        "snp_map",
        "reference_ids",
        "reference_groups",
        "query_ids",
        "query_baseline",
        "alt_all",
        "observed_all",
        "cfg",
        "settings"
      ),
      envir = environment()
    )

    partition_results <- parallel::parLapplyLB(
      cluster,
      tasks,
      function(task) {
        evaluate_partition(
          task = task,
          X = X,
          snp_map = snp_map,
          reference_ids = reference_ids,
          reference_groups = reference_groups,
          query_ids = query_ids,
          query_baseline = query_baseline,
          alt_all = alt_all,
          observed_all = observed_all,
          cfg = cfg,
          settings = settings
        )
      }
    )

    parallel::stopCluster(cluster)
  } else {
    partition_results <- lapply(tasks, function(task) {
      evaluate_partition(
        task = task,
        X = X,
        snp_map = snp_map,
        reference_ids = reference_ids,
        reference_groups = reference_groups,
        query_ids = query_ids,
        query_baseline = query_baseline,
        alt_all = alt_all,
        observed_all = observed_all,
        cfg = cfg,
        settings = settings
      )
    })
  }

  for (result_name in result_names) {
    alternative_results[[result_name]][[length(alternative_results[[result_name]]) + 1L]] <-
      rbindlist(
        lapply(partition_results, `[[`, result_name),
        use.names = TRUE,
        fill = TRUE
      )
  }
}

# =============================================================================
# Combine results and attach fixed Top K comparisons
# =============================================================================

alternatives <- lapply(alternative_results, function(parts) {
  rbindlist(parts, use.names = TRUE, fill = TRUE)
})

baselines <- lapply(baseline_results, function(parts) {
  rbindlist(parts, use.names = TRUE, fill = TRUE)
})

panel_index_table <- rbindlist(panel_index_rows, use.names = TRUE, fill = TRUE)
design <- rbindlist(design_rows, use.names = TRUE, fill = TRUE)
reference_counts <- rbindlist(reference_count_rows, use.names = TRUE, fill = TRUE)
primary_query_baseline <- rbindlist(
  baseline_query_rows,
  use.names = TRUE,
  fill = TRUE
)

baseline_loo_key <- baselines$loo_predictions[, .(
  panel_id,
  K,
  sample_id,
  baseline_fixed_topK_group = predicted_group,
  baseline_fixed_topK_confidence = fixed_panel_confidence,
  baseline_fixed_topK_reason = fixed_panel_reason,
  baseline_fixed_topK_accepted = accepted
)]

alternatives$loo_predictions <- merge(
  alternatives$loo_predictions,
  baseline_loo_key,
  by = c("panel_id", "K", "sample_id"),
  all.x = TRUE,
  sort = FALSE
)

alternatives$loo_predictions[, `:=`(
  matches_baseline_fixed_topK =
    !is.na(predicted_group) & predicted_group == baseline_fixed_topK_group,
  confidence_matches_baseline_fixed_topK =
    fixed_panel_confidence == baseline_fixed_topK_confidence
)]

baseline_query_key <- baselines$query_predictions[, .(
  panel_id,
  K,
  sample_id,
  baseline_fixed_topK_group = predicted_group,
  baseline_fixed_topK_confidence = fixed_panel_confidence,
  baseline_fixed_topK_reason = fixed_panel_reason,
  baseline_fixed_topK_accepted = accepted
)]

alternatives$query_predictions <- merge(
  alternatives$query_predictions,
  baseline_query_key,
  by = c("panel_id", "K", "sample_id"),
  all.x = TRUE,
  sort = FALSE
)

alternatives$query_predictions[, `:=`(
  matches_baseline_fixed_topK =
    !is.na(predicted_group) & predicted_group == baseline_fixed_topK_group,
  confidence_matches_baseline_fixed_topK =
    fixed_panel_confidence == baseline_fixed_topK_confidence
)]

# =============================================================================
# Marker overlap summaries
# =============================================================================

pairwise_overlap_rows <- list()
manifest_split <- split(
  alternatives$marker_manifest,
  by = c("panel_id", "K"),
  keep.by = TRUE,
  drop = TRUE
)

for (block in manifest_split) {
  set_ids <- sort(unique(block$set_id))
  marker_sets <- lapply(set_ids, function(set_id_value) {
    block[set_id == set_id_value, snp_id]
  })
  names(marker_sets) <- set_ids

  set_metadata <- unique(block[, .(
    set_id,
    partition,
    complement_side
  )])

  if (length(set_ids) >= 2L) {
    combinations <- utils::combn(set_ids, 2L, simplify = FALSE)

    for (combination in combinations) {
      set_1 <- combination[1]
      set_2 <- combination[2]
      markers_1 <- marker_sets[[set_1]]
      markers_2 <- marker_sets[[set_2]]
      intersection_n <- length(intersect(markers_1, markers_2))
      union_n <- length(union(markers_1, markers_2))
      metadata_1 <- set_metadata[set_id == set_1][1]
      metadata_2 <- set_metadata[set_id == set_2][1]
      complementary_pair <-
        metadata_1$partition == metadata_2$partition &&
        metadata_1$complement_side != metadata_2$complement_side

      pairwise_overlap_rows[[length(pairwise_overlap_rows) + 1L]] <- data.table(
        panel_id = block$panel_id[1],
        panel_label = block$panel_label[1],
        K = block$K[1],
        set_id_1 = set_1,
        set_id_2 = set_2,
        pair_relation = if (complementary_pair) {
          "complementary_same_partition"
        } else {
          "different_partition"
        },
        intersection_n = intersection_n,
        overlap_fraction_of_K = intersection_n / block$K[1],
        union_n = union_n,
        jaccard = if (union_n > 0L) intersection_n / union_n else NA_real_
      )
    }
  }
}

pairwise_overlap <- rbindlist(
  pairwise_overlap_rows,
  use.names = TRUE,
  fill = TRUE
)

pairwise_overlap_summary <- pairwise_overlap[, .(
  n_set_pairs = .N,
  overlap_fraction_mean = safe_mean(overlap_fraction_of_K),
  overlap_fraction_q025 = safe_quantile(overlap_fraction_of_K, 0.025),
  overlap_fraction_median = safe_quantile(overlap_fraction_of_K, 0.5),
  overlap_fraction_q975 = safe_quantile(overlap_fraction_of_K, 0.975),
  jaccard_mean = safe_mean(jaccard),
  jaccard_q025 = safe_quantile(jaccard, 0.025),
  jaccard_median = safe_quantile(jaccard, 0.5),
  jaccard_q975 = safe_quantile(jaccard, 0.975)
), by = .(panel_id, panel_label, K, pair_relation)]

# =============================================================================
# Performance and query stability summaries
# =============================================================================

baseline_fixed_performance <- copy(baselines$loo_performance)
setorder(baseline_fixed_performance, panel_id, K)

performance_summary <- alternatives$loo_performance[, .(
  n_alternative_sets = .N,
  accuracy_mean = safe_mean(accuracy_all),
  accuracy_sd = safe_sd(accuracy_all),
  accuracy_q025 = safe_quantile(accuracy_all, 0.025),
  accuracy_median = safe_quantile(accuracy_all, 0.5),
  accuracy_q975 = safe_quantile(accuracy_all, 0.975),
  accepted_rate_mean = safe_mean(accepted_rate),
  accepted_rate_q025 = safe_quantile(accepted_rate, 0.025),
  accepted_rate_q975 = safe_quantile(accepted_rate, 0.975),
  accepted_accuracy_mean = safe_mean(accepted_accuracy),
  correct_accepted_rate_mean = safe_mean(correct_accepted_rate),
  uncertainty_rate_mean = safe_mean(uncertainty_rate),
  mean_top_posterior = safe_mean(mean_top_posterior),
  mean_posterior_gap = safe_mean(mean_posterior_gap),
  mean_snps_used = safe_mean(mean_snps_used)
), by = .(panel_id, panel_label, K)]

performance_summary <- merge(
  performance_summary,
  baseline_fixed_performance[, .(
    panel_id,
    K,
    baseline_fixed_topK_accuracy = accuracy_all,
    baseline_fixed_topK_accepted_rate = accepted_rate,
    baseline_fixed_topK_accepted_accuracy = accepted_accuracy,
    baseline_fixed_topK_uncertainty_rate = uncertainty_rate
  )],
  by = c("panel_id", "K"),
  all.x = TRUE,
  sort = FALSE
)

performance_summary[, mean_accuracy_difference_from_fixed_topK :=
  accuracy_mean - baseline_fixed_topK_accuracy]

per_class_summary <- alternatives$loo_per_class[, .(
  n_alternative_sets = .N,
  accuracy_mean = safe_mean(accuracy_all),
  accuracy_sd = safe_sd(accuracy_all),
  accuracy_q025 = safe_quantile(accuracy_all, 0.025),
  accuracy_median = safe_quantile(accuracy_all, 0.5),
  accuracy_q975 = safe_quantile(accuracy_all, 0.975),
  accepted_rate_mean = safe_mean(accepted_rate),
  accepted_accuracy_mean = safe_mean(accepted_accuracy),
  uncertainty_rate_mean = safe_mean(uncertainty_rate)
), by = .(panel_id, panel_label, K, analytical_group)]

query_stability <- alternatives$query_predictions[, .(
  n_alternative_sets = .N,
  primary_baseline_call = primary_baseline_call[1],
  primary_baseline_confidence = primary_baseline_confidence[1],
  baseline_fixed_topK_group = baseline_fixed_topK_group[1],
  baseline_fixed_topK_confidence = baseline_fixed_topK_confidence[1],
  modal_alternative_call = mode_value(predicted_group),
  modal_call_fraction = {
    modal <- mode_value(predicted_group)
    if (is.na(modal)) NA_real_ else mean(predicted_group == modal, na.rm = TRUE)
  },
  primary_multiK_agreement = mean(matches_primary_multiK %in% TRUE),
  fixed_topK_agreement = mean(matches_baseline_fixed_topK %in% TRUE),
  accepted_fraction = mean(accepted),
  uncertainty_fraction = mean(!accepted),
  mean_top_posterior = safe_mean(top_posterior),
  minimum_top_posterior = {
    values <- top_posterior[is.finite(top_posterior)]
    if (length(values) > 0L) min(values) else NA_real_
  },
  mean_posterior_gap = safe_mean(posterior_gap),
  mean_snps_used = safe_mean(n_snps_used)
), by = .(panel_id, panel_label, K, sample_id)]

query_overall_summary <- alternatives$query_performance[, .(
  n_alternative_sets = .N,
  n_queries = n_queries[1],
  primary_call_agreement_mean = safe_mean(primary_call_agreement),
  primary_call_agreement_q025 = safe_quantile(primary_call_agreement, 0.025),
  primary_call_agreement_median = safe_quantile(primary_call_agreement, 0.5),
  primary_call_agreement_q975 = safe_quantile(primary_call_agreement, 0.975),
  all_primary_calls_match_fraction = mean(all_primary_calls_match),
  accepted_rate_mean = safe_mean(accepted_rate),
  uncertainty_rate_mean = safe_mean(uncertainty_rate),
  mean_top_posterior = safe_mean(mean_top_posterior),
  mean_posterior_gap = safe_mean(mean_posterior_gap),
  mean_snps_used = safe_mean(mean_snps_used)
), by = .(panel_id, panel_label, K)]

fixed_query_agreement <- alternatives$query_predictions[, .(
  fixed_topK_call_agreement = mean(matches_baseline_fixed_topK %in% TRUE),
  all_fixed_topK_calls_match = all(matches_baseline_fixed_topK %in% TRUE)
), by = .(panel_id, panel_label, K, set_id)]

fixed_query_agreement_summary <- fixed_query_agreement[, .(
  fixed_topK_call_agreement_mean = safe_mean(fixed_topK_call_agreement),
  fixed_topK_call_agreement_q025 = safe_quantile(fixed_topK_call_agreement, 0.025),
  fixed_topK_call_agreement_q975 = safe_quantile(fixed_topK_call_agreement, 0.975),
  all_fixed_topK_calls_match_fraction = mean(all_fixed_topK_calls_match)
), by = .(panel_id, panel_label, K)]

query_overall_summary <- merge(
  query_overall_summary,
  fixed_query_agreement_summary,
  by = c("panel_id", "panel_label", "K"),
  all.x = TRUE,
  sort = FALSE
)

loo_call_stability <- alternatives$loo_predictions[, .(
  call_agreement_with_fixed_topK = mean(matches_baseline_fixed_topK %in% TRUE),
  confidence_agreement_with_fixed_topK =
    mean(confidence_matches_baseline_fixed_topK %in% TRUE)
), by = .(panel_id, panel_label, K, set_id)]

loo_call_stability_summary <- loo_call_stability[, .(
  n_alternative_sets = .N,
  call_agreement_mean = safe_mean(call_agreement_with_fixed_topK),
  call_agreement_q025 = safe_quantile(call_agreement_with_fixed_topK, 0.025),
  call_agreement_median = safe_quantile(call_agreement_with_fixed_topK, 0.5),
  call_agreement_q975 = safe_quantile(call_agreement_with_fixed_topK, 0.975),
  confidence_agreement_mean = safe_mean(confidence_agreement_with_fixed_topK)
), by = .(panel_id, panel_label, K)]

setorder(performance_summary, panel_id, K)
setorder(per_class_summary, panel_id, K, analytical_group)
setorder(query_stability, panel_id, K, sample_id)
setorder(query_overall_summary, panel_id, K)
setorder(pairwise_overlap_summary, panel_id, K, pair_relation)

# =============================================================================
# Write tables
# =============================================================================

write_table_pair(panel_index_table, "marker_resampling_panel_index", opt$out_dir)
write_table_pair(design, "marker_resampling_design", opt$out_dir)
write_table_pair(reference_counts, "marker_resampling_reference_counts", opt$out_dir)
write_table_pair(primary_query_baseline, "primary_multiK_query_baseline", opt$out_dir)
write_table_pair(baseline_loo, "baseline_individual_loo_summary", opt$out_dir)

write_table_pair(
  baselines$loo_predictions,
  "baseline_fixed_topK_loo_predictions",
  opt$out_dir,
  TRUE
)
write_table_pair(
  baseline_fixed_performance,
  "baseline_fixed_topK_loo_performance",
  opt$out_dir
)
write_table_pair(
  baselines$query_predictions,
  "baseline_fixed_topK_query_predictions",
  opt$out_dir
)

write_table_pair(
  alternatives$loo_predictions,
  "marker_resampling_loo_predictions",
  opt$out_dir,
  TRUE
)
write_table_pair(
  alternatives$loo_performance,
  "marker_resampling_loo_set_performance",
  opt$out_dir
)
write_table_pair(
  performance_summary,
  "marker_resampling_loo_performance_summary",
  opt$out_dir
)
write_table_pair(
  alternatives$loo_per_class,
  "marker_resampling_loo_per_class_set",
  opt$out_dir,
  TRUE
)
write_table_pair(
  per_class_summary,
  "marker_resampling_loo_per_class_summary",
  opt$out_dir
)
write_table_pair(
  loo_call_stability,
  "marker_resampling_loo_call_stability",
  opt$out_dir
)
write_table_pair(
  loo_call_stability_summary,
  "marker_resampling_loo_call_stability_summary",
  opt$out_dir
)

write_table_pair(
  alternatives$query_predictions,
  "marker_resampling_query_predictions",
  opt$out_dir,
  TRUE
)
write_table_pair(
  alternatives$query_performance,
  "marker_resampling_query_set_summary",
  opt$out_dir
)
write_table_pair(
  query_stability,
  "marker_resampling_query_stability",
  opt$out_dir
)
write_table_pair(
  query_overall_summary,
  "marker_resampling_query_overall_summary",
  opt$out_dir
)

write_table_pair(
  alternatives$marker_manifest,
  "marker_resampling_marker_manifest",
  opt$out_dir,
  TRUE
)
write_table_pair(
  alternatives$marker_summary,
  "marker_resampling_marker_set_summary",
  opt$out_dir
)
write_table_pair(
  alternatives$complement_check,
  "marker_resampling_complement_checks",
  opt$out_dir
)
write_table_pair(
  pairwise_overlap,
  "marker_resampling_pairwise_overlap",
  opt$out_dir,
  TRUE
)
write_table_pair(
  pairwise_overlap_summary,
  "marker_resampling_pairwise_overlap_summary",
  opt$out_dir
)

params <- data.table(
  parameter = c(
    "panels",
    "panel_sizes_requested",
    "n_partitions",
    "n_alternative_sets_per_K",
    "n_workers",
    "seed",
    "pseudocount",
    "epsilon",
    "min_snps_query_macroregion",
    "min_snps_query_subregion",
    "high_posterior",
    "moderate_posterior",
    "min_gap_macroregion",
    "min_gap_subregion",
    "marker_resampling_design",
    "fixed_panel_acceptance_scope"
  ),
  value = c(
    paste(panels_requested, collapse = ","),
    paste(panel_sizes_requested, collapse = ","),
    as.character(opt$n_partitions),
    as.character(2L * opt$n_partitions),
    as.character(n_workers),
    as.character(opt$seed),
    as.character(opt$pseudocount),
    as.character(opt$epsilon),
    as.character(opt$min_snps_query_macroregion),
    as.character(opt$min_snps_query_subregion),
    as.character(opt$high_posterior),
    as.character(opt$moderate_posterior),
    as.character(opt$min_gap_macroregion),
    as.character(opt$min_gap_subregion),
    "paired_disjoint_adjacent_rank_partition_of_top_2K",
    paste(
      "posterior, posterior gap, and minimum usable SNP thresholds;",
      "stability is estimated across independent fixed size marker panels"
    )
  )
)

write_table_pair(params, "marker_resampling_params", opt$out_dir)

# =============================================================================
# Base R figure
# =============================================================================

figure_path <- file.path(opt$out_dir, "FIG_marker_resampling_sensitivity.pdf")
grDevices::pdf(figure_path, width = 11.69, height = 8.27, onefile = TRUE)
graphics::par(mfrow = c(length(panels_requested), 2L), mar = c(4.2, 4.4, 2.5, 1.0))

for (panel_id_value in panels_requested) {
  panel_perf <- alternatives$loo_performance[panel_id == panel_id_value]
  panel_baseline <- baseline_fixed_performance[panel_id == panel_id_value]
  panel_query <- alternatives$query_performance[panel_id == panel_id_value]
  panel_label <- unique(panel_perf$panel_label)[1]
  K_values <- sort(unique(panel_perf$K))

  accuracy_list <- lapply(K_values, function(K_value) {
    panel_perf[K == K_value, accuracy_all]
  })

  graphics::boxplot(
    accuracy_list,
    names = K_values,
    ylim = c(0, 1),
    xlab = "Markers per alternative panel",
    ylab = "Reference LOO accuracy",
    main = panel_label,
    col = "#D7E9F7",
    border = "#33658A",
    outline = FALSE
  )

  graphics::points(
    seq_along(K_values),
    panel_baseline[match(K_values, K), accuracy_all],
    pch = 19,
    col = "#D1495B"
  )

  query_means <- vapply(K_values, function(K_value) {
    safe_mean(panel_query[K == K_value, primary_call_agreement])
  }, numeric(1))

  query_low <- vapply(K_values, function(K_value) {
    safe_quantile(panel_query[K == K_value, primary_call_agreement], 0.025)
  }, numeric(1))

  query_high <- vapply(K_values, function(K_value) {
    safe_quantile(panel_query[K == K_value, primary_call_agreement], 0.975)
  }, numeric(1))

  graphics::plot(
    K_values,
    query_means,
    type = "b",
    log = "x",
    ylim = c(0, 1),
    xlab = "Markers per alternative panel",
    ylab = "Agreement with primary query calls",
    main = paste(panel_label, "queries"),
    pch = 19,
    col = "#33658A"
  )
  graphics::segments(K_values, query_low, K_values, query_high, col = "#33658A")
  graphics::abline(h = 1, lty = 2, col = "grey50")
}

grDevices::dev.off()

# =============================================================================
# Plain text summary and run information
# =============================================================================

all_complement_checks_pass <- all(alternatives$complement_check$passed %in% TRUE)
all_marker_sets_have_unique_orthogroups <- all(
  alternatives$marker_summary$n_markers ==
    alternatives$marker_summary$n_unique_orthogroups
)

summary_lines <- c(
  "PESTFLY supplementary validation 06: alternative marker set resampling",
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Design",
  paste0("Complementary partitions per panel and K: ", opt$n_partitions),
  paste0("Alternative fixed size marker panels per K: ", 2L * opt$n_partitions),
  paste0("Requested K values: ", paste(panel_sizes_requested, collapse = ", ")),
  paste0("Parallel workers: ", n_workers),
  paste0("All complementary pair checks passed: ", all_complement_checks_pass),
  paste0(
    "Every alternative panel retained one SNP per orthogroup: ",
    all_marker_sets_have_unique_orthogroups
  ),
  "",
  paste(
    "Each rank matched A/B pair is disjoint, partitions the top 2K corrected",
    "Hudson ranked orthogroups, and gives each set exactly half of the",
    "conventional Top K markers."
  ),
  "",
  "Results by panel and K"
)

for (i in seq_len(nrow(performance_summary))) {
  row <- performance_summary[i]
  query_row <- query_overall_summary[
    panel_id == row$panel_id & K == row$K
  ][1]

  summary_lines <- c(
    summary_lines,
    paste0(
      row$panel_label,
      ", K=", row$K,
      ": alternative LOO accuracy mean=", sprintf("%.4f", row$accuracy_mean),
      " (2.5% to 97.5%=", sprintf("%.4f", row$accuracy_q025),
      " to ", sprintf("%.4f", row$accuracy_q975), ")",
      "; conventional Top K accuracy=",
      sprintf("%.4f", row$baseline_fixed_topK_accuracy),
      "; query agreement with primary calls=",
      sprintf("%.4f", query_row$primary_call_agreement_mean),
      "; query agreement with fixed Top K=",
      sprintf("%.4f", query_row$fixed_topK_call_agreement_mean)
    )
  )
}

summary_lines <- c(
  summary_lines,
  "",
  "Interpretation boundary",
  paste(
    "This analysis tests robustness to nonnested, rank matched alternative",
    "marker panels conditional on the corrected Step 02 ranking. Grouped",
    "geographic CV separately evaluates feature selection and assignment in",
    "held out geographic groups."
  ),
  "",
  "Primary tables",
  "marker_resampling_loo_performance_summary.tsv",
  "marker_resampling_loo_per_class_summary.tsv",
  "marker_resampling_query_stability.tsv",
  "marker_resampling_query_overall_summary.tsv",
  "marker_resampling_pairwise_overlap_summary.tsv",
  "marker_resampling_complement_checks.tsv",
  "FIG_marker_resampling_sensitivity.pdf"
)

summary_path <- file.path(opt$out_dir, "MARKER_RESAMPLING_SUMMARY.txt")
writeLines(summary_lines, summary_path)

run_info <- list(
  analysis = "validation_06_marker_resampling",
  timestamp = Sys.time(),
  repo_root = repo_root,
  inputs = list(
    metadata = metadata_path,
    step2_dir = opt$step2_dir,
    assignment = assignment_path,
    baseline_loo = baseline_loo_path
  ),
  output_dir = opt$out_dir,
  parameters = as.list(opt),
  panels = panel_index_table,
  design = design,
  checks = list(
    all_complement_checks_pass = all_complement_checks_pass,
    all_marker_sets_have_unique_orthogroups =
      all_marker_sets_have_unique_orthogroups
  ),
  session_info = utils::sessionInfo()
)

run_info_path <- file.path(opt$out_dir, "marker_resampling_run_info.rds")
saveRDS(run_info, run_info_path)

# =============================================================================
# Compact return bundle
# =============================================================================

return_files <- c(
  "marker_resampling_panel_index.tsv",
  "marker_resampling_design.tsv",
  "marker_resampling_reference_counts.tsv",
  "primary_multiK_query_baseline.tsv",
  "baseline_individual_loo_summary.tsv",
  "baseline_fixed_topK_loo_performance.tsv",
  "baseline_fixed_topK_query_predictions.tsv",
  "marker_resampling_loo_set_performance.tsv",
  "marker_resampling_loo_performance_summary.tsv",
  "marker_resampling_loo_per_class_summary.tsv",
  "marker_resampling_loo_call_stability_summary.tsv",
  "marker_resampling_query_set_summary.tsv",
  "marker_resampling_query_stability.tsv",
  "marker_resampling_query_overall_summary.tsv",
  "marker_resampling_marker_set_summary.tsv",
  "marker_resampling_complement_checks.tsv",
  "marker_resampling_pairwise_overlap_summary.tsv",
  "marker_resampling_loo_predictions.tsv.gz",
  "marker_resampling_query_predictions.tsv.gz",
  "marker_resampling_marker_manifest.tsv.gz",
  "marker_resampling_pairwise_overlap.tsv.gz",
  "marker_resampling_params.tsv",
  "FIG_marker_resampling_sensitivity.pdf",
  "MARKER_RESAMPLING_SUMMARY.txt",
  "marker_resampling_run_info.rds"
)

for (file_name in return_files) {
  source_path <- file.path(opt$out_dir, file_name)
  if (!file.exists(source_path)) stop("Expected output is missing: ", source_path)

  copied <- file.copy(
    source_path,
    file.path(return_dir, file_name),
    overwrite = TRUE
  )

  if (!isTRUE(copied)) stop("Could not copy output into return_bundle: ", file_name)
}

message("\nMarker resampling validation complete.")
message("Summary: ", summary_path)
message("Figure:  ", figure_path)
message("Return:  ", return_dir)
