# PESTFLY: pure reporting rules for saved Step 04 outputs.
# This file uses base R only. It does not calculate genomic assignments.

read_reporting_parameters <- function(path) {
  if (!file.exists(path)) stop("Missing Step 04 parameter table: ", path)
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = "character", na.strings = c("", "NA"))
  if (!all(c("parameter", "value") %in% names(x))) {
    stop("Step 04 parameter table must contain parameter and value columns")
  }
  if (anyDuplicated(x$parameter)) stop("Duplicate Step 04 parameter names")
  needed <- c("min_snps_query_macroregion", "min_snps_query_subregion",
              "moderate_posterior", "high_posterior", "min_gap_macroregion",
              "min_gap_subregion", "min_agreement", "min_tail_agreement",
              "min_K_available", "tail_fraction")
  if (!all(needed %in% x$parameter)) {
    stop("Missing reporting parameters: ", paste(setdiff(needed, x$parameter), collapse = ", "))
  }
  p <- as.list(setNames(suppressWarnings(as.numeric(x$value[match(needed, x$parameter)])), needed))
  if (any(!is.finite(unlist(p)))) stop("Reporting parameters must be finite numbers")
  integer_params <- c("min_snps_query_macroregion", "min_snps_query_subregion", "min_K_available")
  for (key in integer_params) {
    if (p[[key]] < 1 || p[[key]] != floor(p[[key]])) stop("Invalid count parameter: ", key)
  }
  probability_params <- setdiff(needed, integer_params)
  for (key in probability_params) {
    if (p[[key]] < 0 || p[[key]] > 1) stop("Invalid proportion parameter: ", key)
  }
  if (p$tail_fraction <= 0 || p$high_posterior < p$moderate_posterior) {
    stop("Invalid tail fraction or posterior category thresholds")
  }
  p
}

reporting_threshold_table <- function(p) {
  ids <- c("usable_snps", "posterior", "posterior_gap", "global_K_agreement",
           "tail_K_agreement", "usable_K_count")
  suffixes <- c("maxk_nsnps", "maxk_post", "maxk_gap", "agreement",
                "tail_agreement", "nk_usable_for_stability")
  make <- function(resolution) {
    macro <- resolution == "macroregion"
    data.frame(resolution = resolution, criterion = ids,
               observed_field = paste0(resolution, "_", suffixes),
               operator = ">=", threshold = c(
                 if (macro) p$min_snps_query_macroregion else p$min_snps_query_subregion,
                 p$moderate_posterior,
                 if (macro) p$min_gap_macroregion else p$min_gap_subregion,
                 p$min_agreement, p$min_tail_agreement, p$min_K_available),
               stringsAsFactors = FALSE)
  }
  rbind(make("macroregion"), make("subregion"))
}

assess_reporting_criteria <- function(saved, p) {
  saved <- as.data.frame(saved, stringsAsFactors = FALSE)
  names(saved) <- tolower(names(saved))
  rules <- reporting_threshold_table(p)
  required <- unique(c("sample_id", rules$observed_field,
                       "macroregion_call", "subregion_call", "subregion_panel",
                       "macroregion_confidence", "subregion_confidence",
                       "macroregion_tail_fraction", "subregion_tail_fraction"))
  if (!all(required %in% names(saved))) {
    stop("Saved Step 04 table is missing reporting fields: ",
         paste(setdiff(required, names(saved)), collapse = ", "))
  }
  if (!nrow(saved)) stop("Saved Step 04 table has no samples")
  valid_text <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))
  if (any(!valid_text(saved$sample_id)) || anyDuplicated(saved$sample_id)) {
    stop("Step 04 sample_id values must be populated and unique")
  }
  if (any(!saved$macroregion_call %in% c("Africa", "Asia"))) {
    stop("Macroregion candidates must be Africa or Asia")
  }
  details <- list()
  summaries <- list()
  for (i in seq_len(nrow(saved))) {
    macro_pass <- FALSE
    for (resolution in c("macroregion", "subregion")) {
      macro <- resolution == "macroregion"
      evaluated <- macro || macro_pass
      candidate <- as.character(saved[[paste0(resolution, "_call")]][i])
      panel <- if (macro) "P1_macroregion_africa_vs_asia" else as.character(saved$subregion_panel[i])
      rule <- rules[rules$resolution == resolution, , drop = FALSE]
      if (evaluated) {
        if (!valid_text(candidate) || !valid_text(panel)) {
          stop("Evaluated branch lacks a candidate or panel: ", saved$sample_id[i], " ", resolution)
        }
        expected_panel <- if (saved$macroregion_call[i] == "Africa") "P2_Africa" else "P3_Asia"
        if (!macro && panel != expected_panel) {
          stop("Conditional subregion panel does not match macroregion: ", saved$sample_id[i])
        }
        observed_tail <- suppressWarnings(as.numeric(saved[[paste0(resolution, "_tail_fraction")]][i]))
        if (!is.finite(observed_tail) || observed_tail != p$tail_fraction) {
          stop("Saved tail fraction differs from the parameter table: ", saved$sample_id[i], " ", resolution)
        }
      } else {
        if (valid_text(panel) || valid_text(candidate)) {
          stop("A subregion branch was supplied despite a failed macroregion: ", saved$sample_id[i])
        }
      }
      observed <- vapply(rule$observed_field, function(key) {
        suppressWarnings(as.numeric(saved[[key]][i]))
      }, numeric(1))
      # Proportions outside [0, 1] and fractional/negative counts indicate a bad input.
      finite <- is.finite(observed)
      if (evaluated && any(finite[2:5] & (observed[2:5] < 0 | observed[2:5] > 1))) {
        stop("Invalid reporting proportion: ", saved$sample_id[i], " ", resolution)
      }
      if (evaluated && any(finite[c(1, 6)] &
                          (observed[c(1, 6)] < 0 | observed[c(1, 6)] != floor(observed[c(1, 6)])))) {
        stop("Invalid reporting count: ", saved$sample_id[i], " ", resolution)
      }
      met <- finite & observed >= rule$threshold
      result <- if (evaluated) ifelse(met, "PASS", "FAIL") else rep("NOT_EVALUATED", 6L)
      reason <- if (evaluated) ifelse(!finite, "missing_value",
                                     ifelse(met, "criterion_met", "below_threshold")) else
        rep("macroregion_failed", 6L)
      if (!evaluated) observed[] <- NA_real_
      n_met <- if (evaluated) as.integer(sum(met)) else NA_integer_
      status <- if (!evaluated) "NOT_EVALUATED" else if (n_met == 6L) "PASS" else "FAIL"
      unmet <- if (!evaluated) NA_character_ else if (n_met == 6L) "None" else
        paste(rule$criterion[!met], collapse = "; ")
      confidence <- as.character(saved[[paste0(resolution, "_confidence")]][i])
      accepted <- !is.na(confidence) && tolower(confidence) %in% c("high", "moderate")
      if (accepted != (status == "PASS")) {
        stop("Criterion decision differs from Step 04 confidence eligibility: ",
             saved$sample_id[i], " ", resolution)
      }
      if (status == "PASS") {
        expected_confidence <- if (observed[2] >= p$high_posterior) "High" else "Moderate"
        if (tolower(confidence) != tolower(expected_confidence)) {
          stop("Posterior category differs from Step 04 confidence: ", saved$sample_id[i], " ", resolution)
        }
      }
      details[[length(details) + 1L]] <- data.frame(
        sample_id = saved$sample_id[i], resolution = resolution, panel_id = panel,
        candidate_group = candidate, evaluated = evaluated, criterion = rule$criterion,
        observed_field = rule$observed_field, observed_value = observed,
        operator = rule$operator, threshold = rule$threshold, result = result, reason = reason,
        stringsAsFactors = FALSE)
      summaries[[length(summaries) + 1L]] <- data.frame(
        sample_id = saved$sample_id[i], resolution = resolution, status = status,
        criteria_met = n_met, criteria_unmet = if (evaluated) 6L - n_met else NA_integer_,
        unmet_criteria = unmet, step4_confidence_consistent = TRUE,
        stringsAsFactors = FALSE)
      if (macro) macro_pass <- status == "PASS"
    }
  }
  list(detail = do.call(rbind, details), summary = do.call(rbind, summaries), thresholds = rules)
}
