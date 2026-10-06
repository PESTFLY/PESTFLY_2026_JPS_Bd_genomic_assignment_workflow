#!/usr/bin/env Rscript

# PESTFLY: Reference set audit
#
# Check identifiers, metadata, class representation, ortholog completeness, consensus
# profile similarity and published specimen flags.
# Annotations identify predefined sensitivity sets; the primary snapshot is retained.
#
# Run from the repository root:
#   Rscript steps/validation/03_reference_set_audit/run.R
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

collapse_values <- function(x) {
  x <- sort(unique(clean_text(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(NA_character_)
  paste(x, collapse = ";")
}

safe_fraction <- function(num, den) {
  ifelse(is.finite(den) & den > 0, num / den, NA_real_)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
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

make_flag_rows <- function(
    sample_ids,
    flag_category,
    evidence_domain,
    source,
    source_detail,
    recommended_handling = "retain_primary_and_test_in_sensitivity"
) {
  data.table(
    sample_id = as.character(sample_ids),
    flag_category = flag_category,
    evidence_domain = evidence_domain,
    source = source,
    source_detail = source_detail,
    recommended_handling = recommended_handling
  )
}

# =============================================================================
# Command line options
# =============================================================================

option_list <- list(
  make_option(
    "--qc_dir",
    type = "character",
    default = "results/01_qc",
    help = "Step 01 QC directory [default %default]"
  ),
  make_option(
    "--step2_dir",
    type = "character",
    default = "results/02_snp_panels",
    help = "Corrected Step 02 panel directory [default %default]"
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/validation/03_reference_set_audit",
    help = "Audit output directory [default %default]"
  ),
  make_option(
    "--duplicate_panel",
    type = "character",
    default = "P1_macroregion_africa_vs_asia",
    help = "Panel matrix used for genotype duplicate checks [default %default]"
  ),
  make_option(
    "--min_pairwise_snps",
    type = "integer",
    default = 1000L,
    help = "Minimum pairwise nonmissing SNPs for similarity checks [default %default]"
  ),
  make_option(
    "--near_duplicate_concordance",
    type = "double",
    default = 0.995,
    help = "Concordance threshold for near duplicate flags [default %default]"
  ),
  make_option(
    "--top_pair_rows",
    type = "integer",
    default = 200L,
    help = "Most similar reference pairs retained for audit [default %default]"
  ),
  make_option(
    "--min_fraction_ogs_seen",
    type = "double",
    default = 0.95,
    help = "Ortholog completeness review threshold [default %default]"
  ),
  make_option(
    "--min_ref_per_group",
    type = "integer",
    default = 3L,
    help = "Minimum training references per class for holdout feasibility [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$qc_dir <- resolve_path(opt$qc_dir)
opt$step2_dir <- resolve_path(opt$step2_dir)
opt$out_dir <- resolve_path(opt$out_dir)

opt$min_pairwise_snps <- max(1L, as.integer(opt$min_pairwise_snps))
opt$near_duplicate_concordance <- as.numeric(opt$near_duplicate_concordance)
opt$top_pair_rows <- max(1L, as.integer(opt$top_pair_rows))
opt$min_fraction_ogs_seen <- as.numeric(opt$min_fraction_ogs_seen)
opt$min_ref_per_group <- max(1L, as.integer(opt$min_ref_per_group))

if (!is.finite(opt$near_duplicate_concordance) ||
    opt$near_duplicate_concordance < 0 ||
    opt$near_duplicate_concordance > 1) {
  stop("--near_duplicate_concordance must be between 0 and 1.")
}

if (!is.finite(opt$min_fraction_ogs_seen) ||
    opt$min_fraction_ogs_seen < 0 ||
    opt$min_fraction_ogs_seen > 1) {
  stop("--min_fraction_ogs_seen must be between 0 and 1.")
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Repository root: ", repo_root)
message("QC directory:    ", opt$qc_dir)
message("Step 02 panels:  ", opt$step2_dir)
message("Audit output:    ", opt$out_dir)

# =============================================================================
# Load and clean metadata
# =============================================================================

metadata_path <- file.path(opt$qc_dir, "metadata_clean.tsv")

if (!file.exists(metadata_path)) {
  stop("Missing cleaned metadata: ", metadata_path)
}

meta <- fread(
  metadata_path,
  na.strings = c("", "NA", "NaN", "NULL", "null")
)

setnames(meta, tolower(names(meta)))

required_meta <- c(
  "sample_id",
  "is_reference",
  "macroregion_3",
  "subregion",
  "country",
  "site",
  "collection",
  "sample_accession",
  "organism"
)

missing_meta_columns <- setdiff(required_meta, names(meta))

if (length(missing_meta_columns) > 0L) {
  stop(
    "metadata_clean.tsv is missing required column(s): ",
    paste(missing_meta_columns, collapse = ", ")
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
meta[, collection := clean_text(collection)]
meta[, sample_accession := clean_text(sample_accession)]
meta[, organism := clean_text(organism)]

if (!("collection_date" %in% names(meta))) {
  meta[, collection_date := NA_character_]
}

if (!("population" %in% names(meta))) {
  meta[, population := NA_character_]
}

if (!("geo_loc_name" %in% names(meta))) {
  meta[, geo_loc_name := NA_character_]
}

missing_sample_id_rows <- meta[is.na(sample_id)]

if (nrow(missing_sample_id_rows) > 0L) {
  fwrite(
    missing_sample_id_rows,
    file.path(opt$out_dir, "metadata_rows_missing_sample_id.tsv"),
    sep = "\t",
    na = "NA"
  )
  stop("Metadata contains missing sample identifiers. See audit output.")
}

duplicate_sample_ids <- meta[
  duplicated(sample_id) | duplicated(sample_id, fromLast = TRUE)
][order(sample_id)]

write_table_pair(
  duplicate_sample_ids,
  "metadata_duplicate_sample_ids",
  opt$out_dir
)

if (nrow(duplicate_sample_ids) > 0L) {
  stop("Metadata contains duplicated sample identifiers. See audit output.")
}

refs <- copy(meta[is_reference %in% TRUE])
queries <- copy(meta[is_reference %in% FALSE])

if (nrow(refs) == 0L) {
  stop("No reference samples were found in metadata_clean.tsv.")
}

message("Reference samples: ", nrow(refs))
message("Query samples:     ", nrow(queries))

# =============================================================================
# Metadata identity and completeness audit
# =============================================================================

duplicate_accessions <- refs[
  !is.na(sample_accession) &
    (duplicated(sample_accession) | duplicated(sample_accession, fromLast = TRUE))
][order(sample_accession, sample_id)]

write_table_pair(
  duplicate_accessions,
  "metadata_duplicate_accessions",
  opt$out_dir
)

issue_list <- list()

add_missing_issue <- function(field, severity, issue_type = "missing_value") {
  idx <- which(is.na(refs[[field]]))
  if (length(idx) == 0L) return(invisible(NULL))

  issue_list[[length(issue_list) + 1L]] <<- data.table(
    sample_id = refs$sample_id[idx],
    field = field,
    issue_type = issue_type,
    severity = severity,
    observed_value = NA_character_,
    detail = paste0("Reference metadata field '", field, "' is missing")
  )

  invisible(NULL)
}

add_missing_issue("sample_accession", "critical")
add_missing_issue("organism", "high")
add_missing_issue("macroregion_3", "high")
add_missing_issue("subregion", "high")
add_missing_issue("country", "high")
add_missing_issue("collection", "high")
add_missing_issue("site", "documentation")
add_missing_issue("collection_date", "documentation")

bad_organism <- refs[
  !is.na(organism) & tolower(organism) != "bactrocera dorsalis"
]

if (nrow(bad_organism) > 0L) {
  issue_list[[length(issue_list) + 1L]] <- bad_organism[, .(
    sample_id,
    field = "organism",
    issue_type = "unexpected_organism_label",
    severity = "critical",
    observed_value = organism,
    detail = "Reference organism label is not Bactrocera dorsalis"
  )]
}

valid_subregion_map <- list(
  Africa = c("C_Africa", "E_Africa", "S_Africa", "W_Africa"),
  Asia = c("E_Asia", "SE_Asia", "S_Asia"),
  Others = "Other"
)

for (region in names(valid_subregion_map)) {
  bad <- refs[
    macroregion_3 == region &
      !is.na(subregion) &
      !(subregion %in% valid_subregion_map[[region]])
  ]

  if (nrow(bad) > 0L) {
    issue_list[[length(issue_list) + 1L]] <- bad[, .(
      sample_id,
      field = "macroregion_3/subregion",
      issue_type = "incompatible_class_mapping",
      severity = "critical",
      observed_value = paste(macroregion_3, subregion, sep = "/"),
      detail = "Subregion is incompatible with macroregion_3"
    )]
  }
}

unexpected_regions <- refs[
  !is.na(macroregion_3) & !(macroregion_3 %in% names(valid_subregion_map))
]

if (nrow(unexpected_regions) > 0L) {
  issue_list[[length(issue_list) + 1L]] <- unexpected_regions[, .(
    sample_id,
    field = "macroregion_3",
    issue_type = "unexpected_class_label",
    severity = "critical",
    observed_value = macroregion_3,
    detail = "macroregion_3 is not Africa, Asia, or Others"
  )]
}

if (nrow(duplicate_accessions) > 0L) {
  issue_list[[length(issue_list) + 1L]] <- duplicate_accessions[, .(
    sample_id,
    field = "sample_accession",
    issue_type = "duplicated_accession",
    severity = "critical",
    observed_value = sample_accession,
    detail = "sample_accession occurs more than once among references"
  )]
}

metadata_issues <- if (length(issue_list) > 0L) {
  unique(rbindlist(issue_list, use.names = TRUE, fill = TRUE))
} else {
  data.table(
    sample_id = character(),
    field = character(),
    issue_type = character(),
    severity = character(),
    observed_value = character(),
    detail = character()
  )
}

setorder(metadata_issues, severity, field, sample_id)

write_table_pair(metadata_issues, "metadata_issues", opt$out_dir)

metadata_issue_summary <- metadata_issues[, .N, by = .(
  severity,
  field,
  issue_type
)][order(severity, field, issue_type)]

write_table_pair(
  metadata_issue_summary,
  "metadata_issue_summary",
  opt$out_dir
)

# =============================================================================
# Reference composition and panel membership
# =============================================================================

reference_class_counts <- refs[, .N, by = .(
  macroregion_3,
  subregion
)][order(macroregion_3, subregion)]

reference_collection_counts <- refs[, .N, by = .(
  collection,
  macroregion_3,
  subregion
)][order(collection, macroregion_3, subregion)]

reference_geographic_counts <- refs[, .N, by = .(
  macroregion_3,
  subregion,
  country,
  site,
  collection
)][order(macroregion_3, subregion, country, site, collection)]

write_table_pair(
  reference_class_counts,
  "reference_class_counts",
  opt$out_dir
)

write_table_pair(
  reference_collection_counts,
  "reference_collection_counts",
  opt$out_dir
)

write_table_pair(
  reference_geographic_counts,
  "reference_geographic_counts",
  opt$out_dir
)

panel_definitions <- list(
  P1_macroregion_africa_vs_asia = list(
    panel_label = "P1 macroregion",
    group_col = "macroregion_3",
    macroregion_filter = c("Africa", "Asia")
  ),
  P2_subregion_within_africa = list(
    panel_label = "P2 Africa subregion",
    group_col = "subregion",
    macroregion_filter = "Africa"
  ),
  P3_subregion_within_asia = list(
    panel_label = "P3 Asia subregion",
    group_col = "subregion",
    macroregion_filter = "Asia"
  )
)

panel_membership_list <- list()

for (panel_id in names(panel_definitions)) {
  cfg <- panel_definitions[[panel_id]]
  d <- refs[
    macroregion_3 %in% cfg$macroregion_filter &
      !is.na(get(cfg$group_col))
  ]

  panel_membership_list[[panel_id]] <- d[, .(
    panel_id = panel_id,
    panel_label = cfg$panel_label,
    group_col = cfg$group_col,
    sample_id,
    analytical_group = get(cfg$group_col),
    macroregion_3,
    subregion,
    country,
    site,
    collection,
    sample_accession
  )]
}

panel_membership <- rbindlist(
  panel_membership_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(panel_membership, panel_id, analytical_group, country, site, sample_id)

panel_reference_counts <- panel_membership[, .N, by = .(
  panel_id,
  panel_label,
  group_col,
  analytical_group
)][order(panel_id, analytical_group)]

write_table_pair(
  panel_membership,
  "panel_reference_membership",
  opt$out_dir
)

write_table_pair(
  panel_reference_counts,
  "panel_reference_counts",
  opt$out_dir
)

# =============================================================================
# Published specimen flags
# =============================================================================

# Exact list from Vanbergen et al. 2025, Supplementary Table S9.
vanbergen_pca_outliers <- c(
  "GXNN08",
  "YNJH01", "YNJH02", "YNJH03", "YNJH05",
  "YNJH07", "YNJH08", "YNJH09", "YNJH10",
  "INTG06",
  "IDJI04", "IDJI10",
  "MMYG01", "MMYG02", "MMYG04", "MMYG05", "MMYG06",
  "MMYG07", "MMYG08", "MMYG10",
  "PHDM09",
  "LKAD07"
)

# Hybrid candidates identified by Zhang et al. and reconciled to the corrected
# identifier note in Vanbergen et al. Supplementary Table S9.  The Zhang main
# text used IDJI06, while its Figure 2 and the later audit identify IDJI10.
zhang_hybrid_candidates <- c(
  "IDJI04", "IDJI10",
  "MMYG01", "MMYG02", "MMYG04", "MMYG05", "MMYG06", "MMYG08", "MMYG10"
)

# WGS specimens in the COI cluster that produced B. dorsalis and
# B. kandiensis matches within the reporting threshold.
coi_taxonomic_caution <- c(
  "CDKB02", "LKAD01", "LKAD04", "LKAD06", "LKAD07", "LKAD08"
)

# Other explicit sample level cautions reported by Vanbergen et al.
high_missingness_sample <- "BIBJ10"
africa_asia_discordant_samples <- c("ZAMP06", "ZAMP07", "ZAMP09")

published_flags <- rbindlist(
  list(
    make_flag_rows(
      vanbergen_pca_outliers,
      "nuclear_PCA_outlier",
      "nuclear_genome",
      "Vanbergen et al. 2025, Supplementary Table S9",
      "Did not cluster with the main B. dorsalis nuclear PCA groups"
    ),
    make_flag_rows(
      zhang_hybrid_candidates,
      "published_hybrid_candidate",
      "nuclear_genome",
      "Zhang et al. 2023; identifier reconciliation in Vanbergen et al. 2025 Table S9",
      "Reported as showing B. dorsalis and B. carambolae hybrid ancestry"
    ),
    make_flag_rows(
      coi_taxonomic_caution,
      "COI_dorsalis_kandiensis_ambiguity",
      "mitochondrial_COI",
      "Vanbergen et al. 2025, Supplementary Table S15",
      "COI similarity results included both B. dorsalis and B. kandiensis",
      "retain_primary_and_treat_as_extended_sensitivity_only"
    ),
    make_flag_rows(
      high_missingness_sample,
      "published_high_nuclear_missingness",
      "nuclear_genome",
      "Vanbergen et al. 2025, Supplementary Table S10",
      "Reported nuclear missingness of 96.7 percent",
      "exclude_if_present_due_to_data_quality"
    ),
    make_flag_rows(
      africa_asia_discordant_samples,
      "published_geography_discordance",
      "nuclear_and_mitochondrial",
      "Vanbergen et al. 2025, Figure 2 and mitogenome results",
      "South African specimens clustered with predominantly Asian profiles"
    ),
    make_flag_rows(
      "LKAD08",
      "published_identifier_discrepancy",
      "publication_audit",
      "Vanbergen et al. 2025 main text versus Supplementary Table S9",
      "Main text names LKAD08 as a PCA outlier, while Table S9 lists LKAD07",
      "document_discrepancy_no_automatic_exclusion"
    ),
    make_flag_rows(
      "IDJI10",
      "published_identifier_discrepancy",
      "publication_audit",
      "Vanbergen et al. 2025, Supplementary Table S9",
      "Table S9 states that Zhang et al. main text used IDJI06 but Figure 2 supports IDJI10",
      "use_corrected_IDJI10_annotation"
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

published_flags <- unique(published_flags)
setorder(published_flags, sample_id, flag_category)

flagged_specimen_status <- merge(
  published_flags,
  meta[, .(
    sample_id,
    present_in_snapshot = TRUE,
    is_reference,
    macroregion_3,
    subregion,
    country,
    site,
    collection,
    sample_accession,
    organism
  )],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

flagged_specimen_status[
  is.na(present_in_snapshot),
  present_in_snapshot := FALSE
]

flagged_specimen_status[, snapshot_status := fifelse(
  !present_in_snapshot,
  "absent_from_frozen_snapshot",
  fifelse(
    is_reference %in% TRUE,
    "retained_reference",
    "present_not_reference"
  )
)]

setorder(flagged_specimen_status, snapshot_status, sample_id, flag_category)

write_table_pair(
  published_flags,
  "published_specimen_flags",
  opt$out_dir
)

write_table_pair(
  flagged_specimen_status,
  "published_flag_status_in_snapshot",
  opt$out_dir
)

retained_published_flags <- flagged_specimen_status[
  snapshot_status == "retained_reference"
]

write_table_pair(
  retained_published_flags,
  "retained_reference_published_flags",
  opt$out_dir
)

reference_flag_summary <- retained_published_flags[, .(
  n_published_flags = uniqueN(flag_category),
  flag_categories = collapse_values(flag_category),
  evidence_domains = collapse_values(evidence_domain),
  sources = collapse_values(source)
), by = .(
  sample_id,
  macroregion_3,
  subregion,
  country,
  site,
  collection,
  sample_accession
)][order(macroregion_3, subregion, country, sample_id)]

write_table_pair(
  reference_flag_summary,
  "retained_reference_flag_summary",
  opt$out_dir
)

# =============================================================================
# Predeclared sensitivity tiers and remaining class counts
# =============================================================================

retained_ids <- refs$sample_id

scenario_sets <- list(
  primary_all_references = character(),
  exclude_published_hybrid_candidates = intersect(
    retained_ids,
    zhang_hybrid_candidates
  ),
  exclude_all_nuclear_PCA_outliers = intersect(
    retained_ids,
    vanbergen_pca_outliers
  ),
  exclude_extended_taxonomic_caution = intersect(
    retained_ids,
    union(vanbergen_pca_outliers, coi_taxonomic_caution)
  )
)

sensitivity_scenario_membership <- rbindlist(
  lapply(names(scenario_sets), function(scenario_name) {
    ids <- scenario_sets[[scenario_name]]

    if (length(ids) == 0L) {
      return(data.table(
        scenario = scenario_name,
        sample_id = NA_character_,
        exclusion_order = NA_integer_
      ))
    }

    data.table(
      scenario = scenario_name,
      sample_id = sort(ids),
      exclusion_order = seq_along(sort(ids))
    )
  }),
  use.names = TRUE,
  fill = TRUE
)

sensitivity_scenario_membership <- merge(
  sensitivity_scenario_membership,
  refs[, .(
    sample_id,
    macroregion_3,
    subregion,
    country,
    site,
    collection,
    sample_accession
  )],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(sensitivity_scenario_membership, scenario, exclusion_order)

write_table_pair(
  sensitivity_scenario_membership,
  "sensitivity_scenario_membership",
  opt$out_dir
)

scenario_count_list <- list()

for (scenario_name in names(scenario_sets)) {
  excluded_ids <- scenario_sets[[scenario_name]]

  for (panel_id in names(panel_definitions)) {
    cfg <- panel_definitions[[panel_id]]

    d <- refs[
      macroregion_3 %in% cfg$macroregion_filter &
        !is.na(get(cfg$group_col))
    ]

    before <- d[, .(n_before = .N), by = .(
      analytical_group = get(cfg$group_col)
    )]

    after <- d[!(sample_id %in% excluded_ids), .(
      n_after = .N
    ), by = .(
      analytical_group = get(cfg$group_col)
    )]

    counts <- merge(
      before,
      after,
      by = "analytical_group",
      all = TRUE,
      sort = TRUE
    )

    counts[is.na(n_before), n_before := 0L]
    counts[is.na(n_after), n_after := 0L]
    counts[, n_excluded := n_before - n_after]
    counts[, `:=`(
      scenario = scenario_name,
      panel_id = panel_id,
      panel_label = cfg$panel_label,
      group_col = cfg$group_col,
      class_supported_after_exclusion = n_after >= opt$min_ref_per_group
    )]

    scenario_count_list[[length(scenario_count_list) + 1L]] <- counts
  }
}

sensitivity_scenario_counts <- rbindlist(
  scenario_count_list,
  use.names = TRUE,
  fill = TRUE
)

setcolorder(
  sensitivity_scenario_counts,
  c(
    "scenario",
    "panel_id",
    "panel_label",
    "group_col",
    "analytical_group",
    "n_before",
    "n_excluded",
    "n_after",
    "class_supported_after_exclusion"
  )
)

setorder(
  sensitivity_scenario_counts,
  scenario,
  panel_id,
  analytical_group
)

write_table_pair(
  sensitivity_scenario_counts,
  "sensitivity_scenario_counts",
  opt$out_dir
)

# =============================================================================
# Ortholog coverage audit
# =============================================================================

ortholog_summary_path <- file.path(
  opt$qc_dir,
  "label_consistency_sample_summary_phy.tsv"
)

if (!file.exists(ortholog_summary_path)) {
  stop("Missing sample ortholog summary: ", ortholog_summary_path)
}

ortholog_summary <- fread(ortholog_summary_path)
setnames(ortholog_summary, tolower(names(ortholog_summary)))

required_ortholog_columns <- c(
  "sample_id",
  "n_ogs_total",
  "n_ogs_missing",
  "n_ogs_seen",
  "frac_ogs_seen"
)

missing_ortholog_columns <- setdiff(
  required_ortholog_columns,
  names(ortholog_summary)
)

if (length(missing_ortholog_columns) > 0L) {
  stop(
    "Ortholog summary is missing column(s): ",
    paste(missing_ortholog_columns, collapse = ", ")
  )
}

ortholog_summary[, sample_id := clean_text(sample_id)]

reference_ortholog_coverage <- merge(
  refs[, .(
    sample_id,
    macroregion_3,
    subregion,
    country,
    site,
    collection,
    sample_accession
  )],
  ortholog_summary[, .(
    sample_id,
    n_ogs_total,
    n_ogs_missing,
    n_ogs_seen,
    frac_ogs_seen
  )],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

reference_ortholog_coverage[, ortholog_summary_present := !is.na(n_ogs_total)]
reference_ortholog_coverage[, below_completeness_threshold :=
  !is.na(frac_ogs_seen) & frac_ogs_seen < opt$min_fraction_ogs_seen]

setorder(
  reference_ortholog_coverage,
  frac_ogs_seen,
  macroregion_3,
  subregion,
  sample_id,
  na.last = FALSE
)

ortholog_coverage_by_collection <- reference_ortholog_coverage[, .(
  n_references = .N,
  n_missing_summary = sum(!ortholog_summary_present),
  n_below_threshold = sum(below_completeness_threshold, na.rm = TRUE),
  minimum_fraction_ogs_seen = safe_min(frac_ogs_seen),
  median_fraction_ogs_seen = safe_median(frac_ogs_seen),
  mean_fraction_ogs_seen = safe_mean(frac_ogs_seen),
  maximum_fraction_ogs_seen = safe_max(frac_ogs_seen)
), by = .(
  collection,
  macroregion_3,
  subregion
)][order(collection, macroregion_3, subregion)]

write_table_pair(
  reference_ortholog_coverage,
  "reference_ortholog_coverage",
  opt$out_dir
)

write_table_pair(
  ortholog_coverage_by_collection,
  "ortholog_coverage_by_collection",
  opt$out_dir
)

# =============================================================================
# Study and geographic support audit
# =============================================================================

study_class_count_list <- list()
study_class_summary_list <- list()
geographic_support_list <- list()

for (panel_id in names(panel_definitions)) {
  cfg <- panel_definitions[[panel_id]]

  # Store the loop values under names that cannot be confused with columns.
  # data.table does not recycle one value expressions supplied directly in a
  # by list on all supported versions, so add these labels to d first.
  panel_id_value <- panel_id
  panel_label_value <- cfg$panel_label
  group_col_value <- cfg$group_col

  d <- refs[
    macroregion_3 %in% cfg$macroregion_filter &
      !is.na(get(group_col_value))
  ]

  d[, `:=`(
    panel_id = panel_id_value,
    panel_label = panel_label_value,
    group_col = group_col_value,
    analytical_group = get(group_col_value)
  )]

  counts <- d[, .N, by = .(
    panel_id,
    panel_label,
    group_col,
    analytical_group,
    collection
  )]

  totals <- counts[, .(class_total = sum(N)), by = .(
    panel_id,
    panel_label,
    group_col,
    analytical_group
  )]

  counts <- merge(
    counts,
    totals,
    by = c("panel_id", "panel_label", "group_col", "analytical_group"),
    all.x = TRUE,
    sort = FALSE
  )

  counts[, class_fraction_from_collection := safe_fraction(N, class_total)]

  summary <- counts[order(-N, collection), .(
    class_total = class_total[1],
    n_source_collections = uniqueN(collection),
    dominant_collection = collection[1],
    dominant_collection_n = N[1],
    dominant_collection_fraction = class_fraction_from_collection[1],
    single_collection_class = uniqueN(collection) == 1L
  ), by = .(
    panel_id,
    panel_label,
    group_col,
    analytical_group
  )]

  geographic <- d[, .(
    n_references = .N,
    n_countries = uniqueN(country, na.rm = TRUE),
    n_sites = uniqueN(paste(country, site, sep = "::"), na.rm = TRUE),
    countries = collapse_values(country),
    source_collections = collapse_values(collection)
  ), by = .(
    panel_id,
    panel_label,
    group_col,
    analytical_group
  )]

  study_class_count_list[[panel_id]] <- counts
  study_class_summary_list[[panel_id]] <- summary
  geographic_support_list[[panel_id]] <- geographic
}

study_class_counts <- rbindlist(
  study_class_count_list,
  use.names = TRUE,
  fill = TRUE
)

study_class_summary <- rbindlist(
  study_class_summary_list,
  use.names = TRUE,
  fill = TRUE
)

geographic_support_summary <- rbindlist(
  geographic_support_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(study_class_counts, panel_id, analytical_group, collection)
setorder(study_class_summary, panel_id, analytical_group)
setorder(geographic_support_summary, panel_id, analytical_group)

write_table_pair(
  study_class_counts,
  "study_class_counts",
  opt$out_dir
)

write_table_pair(
  study_class_summary,
  "study_class_summary",
  opt$out_dir
)

write_table_pair(
  geographic_support_summary,
  "geographic_support_summary",
  opt$out_dir
)

# =============================================================================
# Leave one source collection out feasibility
# =============================================================================

collection_holdout_fold_list <- list()
collection_holdout_count_list <- list()

for (panel_id in names(panel_definitions)) {
  cfg <- panel_definitions[[panel_id]]

  d <- refs[
    macroregion_3 %in% cfg$macroregion_filter &
      !is.na(get(cfg$group_col)) &
      !is.na(collection)
  ]

  d[, analytical_group := get(cfg$group_col)]
  full_groups <- sort(unique(d$analytical_group))

  for (heldout_collection in sort(unique(d$collection))) {
    test <- d[collection == heldout_collection]
    train <- d[collection != heldout_collection]

    observed <- train[, .N, by = analytical_group]

    counts <- merge(
      data.table(analytical_group = full_groups),
      observed,
      by = "analytical_group",
      all.x = TRUE,
      sort = TRUE
    )

    counts[is.na(N), N := 0L]
    counts[, supported := N >= opt$min_ref_per_group]
    counts[, `:=`(
      panel_id = panel_id,
      panel_label = cfg$panel_label,
      group_col = cfg$group_col,
      heldout_collection = heldout_collection
    )]

    supported_groups <- counts[supported %in% TRUE, analytical_group]
    true_groups <- sort(unique(test$analytical_group))

    fold_reason <- "ok"

    if (length(supported_groups) < 2L) {
      fold_reason <- "fewer_than_two_supported_training_classes"
    } else if (!all(true_groups %in% supported_groups)) {
      fold_reason <- "true_class_absent_or_undersupported_after_holdout"
    }

    fold <- data.table(
      panel_id = panel_id,
      panel_label = cfg$panel_label,
      group_col = cfg$group_col,
      heldout_collection = heldout_collection,
      heldout_true_groups = collapse_values(true_groups),
      n_test = nrow(test),
      n_train = nrow(train),
      n_supported_training_classes = length(supported_groups),
      supported_training_classes = collapse_values(supported_groups),
      fold_evaluable = identical(fold_reason, "ok"),
      fold_reason = fold_reason
    )

    collection_holdout_fold_list[[length(collection_holdout_fold_list) + 1L]] <- fold
    collection_holdout_count_list[[length(collection_holdout_count_list) + 1L]] <- counts
  }
}

collection_holdout_folds <- rbindlist(
  collection_holdout_fold_list,
  use.names = TRUE,
  fill = TRUE
)

collection_holdout_training_counts <- rbindlist(
  collection_holdout_count_list,
  use.names = TRUE,
  fill = TRUE
)

setcolorder(
  collection_holdout_training_counts,
  c(
    "panel_id",
    "panel_label",
    "group_col",
    "heldout_collection",
    "analytical_group",
    "N",
    "supported"
  )
)

setorder(collection_holdout_folds, panel_id, heldout_collection)
setorder(
  collection_holdout_training_counts,
  panel_id,
  heldout_collection,
  analytical_group
)

write_table_pair(
  collection_holdout_folds,
  "collection_holdout_feasibility",
  opt$out_dir
)

write_table_pair(
  collection_holdout_training_counts,
  "collection_holdout_training_counts",
  opt$out_dir
)

# =============================================================================
# P1 consensus call completeness and duplicate profile audit
# =============================================================================

duplicate_panel_dir <- file.path(
  opt$step2_dir,
  "panels",
  opt$duplicate_panel
)

duplicate_matrix_path <- file.path(duplicate_panel_dir, "snp_matrix.rds")

if (!file.exists(duplicate_matrix_path)) {
  stop("Missing SNP matrix for duplicate audit: ", duplicate_matrix_path)
}

X <- readRDS(duplicate_matrix_path)

if (!is.matrix(X)) X <- as.matrix(X)
storage.mode(X) <- "integer"

if (is.null(rownames(X)) || is.null(colnames(X))) {
  stop("Duplicate audit SNP matrix requires SNP row names and sample columns.")
}

missing_reference_columns <- setdiff(refs$sample_id, colnames(X))

if (length(missing_reference_columns) > 0L) {
  stop(
    "Duplicate audit matrix is missing reference sample(s): ",
    paste(head(missing_reference_columns, 20L), collapse = ", ")
  )
}

X_ref <- X[, refs$sample_id, drop = FALSE]

observed_values <- sort(unique(as.vector(X_ref)))
bad_values <- observed_values[
  !is.na(observed_values) & !(observed_values %in% c(0L, 1L))
]

if (length(bad_values) > 0L) {
  stop(
    "Duplicate audit expects binary 0, 1, or missing calls. Unexpected value(s): ",
    paste(bad_values, collapse = ", ")
  )
}

genotype_completeness <- data.table(
  sample_id = colnames(X_ref),
  n_panel_snps = nrow(X_ref),
  n_genotypes_called = colSums(!is.na(X_ref)),
  n_genotypes_missing = colSums(is.na(X_ref))
)

genotype_completeness[, fraction_genotypes_called := safe_fraction(
  n_genotypes_called,
  n_panel_snps
)]

genotype_completeness <- merge(
  genotype_completeness,
  refs[, .(
    sample_id,
    macroregion_3,
    subregion,
    country,
    site,
    collection,
    sample_accession
  )],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

setorder(genotype_completeness, fraction_genotypes_called, sample_id)

write_table_pair(
  genotype_completeness,
  "reference_P1_genotype_completeness",
  opt$out_dir
)

message(
  "Calculating pairwise SNP profile concordance for ",
  ncol(X_ref),
  " references across ",
  nrow(X_ref),
  " P1 SNPs..."
)

valid <- 1.0 * (!is.na(X_ref))
G <- X_ref
G[is.na(G)] <- 0L
storage.mode(G) <- "double"

n_compared_matrix <- crossprod(valid)
ones_valid_matrix <- crossprod(G, valid)
ones_both_matrix <- crossprod(G)

mismatch_matrix <- ones_valid_matrix +
  t(ones_valid_matrix) -
  2 * ones_both_matrix

mismatch_matrix[abs(mismatch_matrix) < 1e-8] <- 0
mismatch_matrix <- round(pmax(mismatch_matrix, 0))

n_valid_per_sample <- diag(n_compared_matrix)
upper_idx <- which(upper.tri(n_compared_matrix), arr.ind = TRUE)

pairwise_similarity <- data.table(
  sample_a = colnames(X_ref)[upper_idx[, 1]],
  sample_b = colnames(X_ref)[upper_idx[, 2]],
  n_compared = as.integer(n_compared_matrix[upper_idx]),
  n_mismatches = as.integer(mismatch_matrix[upper_idx]),
  n_valid_a = as.integer(n_valid_per_sample[upper_idx[, 1]]),
  n_valid_b = as.integer(n_valid_per_sample[upper_idx[, 2]])
)

pairwise_similarity[, concordance := fifelse(
  n_compared > 0L,
  1 - n_mismatches / n_compared,
  NA_real_
)]

pairwise_similarity[, genotype_identical_on_overlap :=
  n_compared > 0L & n_mismatches == 0L]

pairwise_similarity[, full_profile_identical :=
  genotype_identical_on_overlap &
    n_compared == n_valid_a &
    n_compared == n_valid_b]

meta_a <- refs[, .(
  sample_a = sample_id,
  accession_a = sample_accession,
  macroregion_a = macroregion_3,
  subregion_a = subregion,
  country_a = country,
  site_a = site,
  collection_a = collection
)]

meta_b <- refs[, .(
  sample_b = sample_id,
  accession_b = sample_accession,
  macroregion_b = macroregion_3,
  subregion_b = subregion,
  country_b = country,
  site_b = site,
  collection_b = collection
)]

pairwise_similarity <- merge(
  pairwise_similarity,
  meta_a,
  by = "sample_a",
  all.x = TRUE,
  sort = FALSE
)

pairwise_similarity <- merge(
  pairwise_similarity,
  meta_b,
  by = "sample_b",
  all.x = TRUE,
  sort = FALSE
)

pairwise_similarity[, `:=`(
  same_macroregion = macroregion_a == macroregion_b,
  same_subregion = subregion_a == subregion_b,
  same_country = country_a == country_b,
  same_site = country_a == country_b & site_a == site_b,
  same_collection = collection_a == collection_b
)]

pairwise_similarity_eligible <- pairwise_similarity[
  n_compared >= opt$min_pairwise_snps & !is.na(concordance)
]

setorder(
  pairwise_similarity_eligible,
  -concordance,
  -n_compared,
  n_mismatches,
  sample_a,
  sample_b
)

top_pairwise_similarity <- head(
  pairwise_similarity_eligible,
  opt$top_pair_rows
)

near_duplicate_pairs <- pairwise_similarity_eligible[
  concordance >= opt$near_duplicate_concordance
]

exact_overlap_pairs <- pairwise_similarity_eligible[
  genotype_identical_on_overlap %in% TRUE
]

full_profile_identical_pairs <- pairwise_similarity_eligible[
  full_profile_identical %in% TRUE
]

write_table_pair(
  top_pairwise_similarity,
  "reference_pairwise_similarity_top",
  opt$out_dir
)

write_table_pair(
  near_duplicate_pairs,
  "reference_near_duplicate_pairs",
  opt$out_dir
)

write_table_pair(
  exact_overlap_pairs,
  "reference_identical_on_overlap_pairs",
  opt$out_dir
)

write_table_pair(
  full_profile_identical_pairs,
  "reference_full_profile_identical_pairs",
  opt$out_dir
)

# =============================================================================
# Consolidated sample audit table
# =============================================================================

flag_aggregate <- retained_published_flags[, .(
  n_published_flags = uniqueN(flag_category),
  published_flags = collapse_values(flag_category)
), by = sample_id]

issue_aggregate <- metadata_issues[, .(
  n_metadata_issues = .N,
  metadata_issue_types = collapse_values(issue_type),
  highest_metadata_severity = fifelse(
    any(severity == "critical"),
    "critical",
    fifelse(
      any(severity == "high"),
      "high",
      fifelse(any(severity == "documentation"), "documentation", NA_character_)
    )
  )
), by = sample_id]

reference_sample_audit <- merge(
  refs[, .(
    sample_id,
    sample_accession,
    organism,
    macroregion_3,
    subregion,
    country,
    site,
    collection,
    collection_date,
    population,
    geo_loc_name
  )],
  reference_ortholog_coverage[, .(
    sample_id,
    n_ogs_total,
    n_ogs_missing,
    n_ogs_seen,
    frac_ogs_seen,
    below_completeness_threshold
  )],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

reference_sample_audit <- merge(
  reference_sample_audit,
  genotype_completeness[, .(
    sample_id,
    n_panel_snps,
    n_genotypes_called,
    n_genotypes_missing,
    fraction_genotypes_called
  )],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

reference_sample_audit <- merge(
  reference_sample_audit,
  flag_aggregate,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

reference_sample_audit <- merge(
  reference_sample_audit,
  issue_aggregate,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

reference_sample_audit[is.na(n_published_flags), n_published_flags := 0L]
reference_sample_audit[is.na(n_metadata_issues), n_metadata_issues := 0L]

reference_sample_audit[, in_hybrid_exclusion_tier :=
  sample_id %in% scenario_sets$exclude_published_hybrid_candidates]

reference_sample_audit[, in_nuclear_outlier_exclusion_tier :=
  sample_id %in% scenario_sets$exclude_all_nuclear_PCA_outliers]

reference_sample_audit[, in_extended_caution_exclusion_tier :=
  sample_id %in% scenario_sets$exclude_extended_taxonomic_caution]

setorder(
  reference_sample_audit,
  -n_published_flags,
  -n_metadata_issues,
  macroregion_3,
  subregion,
  country,
  sample_id
)

write_table_pair(
  reference_sample_audit,
  "reference_sample_audit",
  opt$out_dir
)

# =============================================================================
# Audit summary
# =============================================================================

n_pca_outlier_refs <- uniqueN(
  retained_published_flags[
    flag_category == "nuclear_PCA_outlier",
    sample_id
  ]
)

n_hybrid_candidate_refs <- uniqueN(
  retained_published_flags[
    flag_category == "published_hybrid_candidate",
    sample_id
  ]
)

n_coi_caution_refs <- uniqueN(
  retained_published_flags[
    flag_category == "COI_dorsalis_kandiensis_ambiguity",
    sample_id
  ]
)

n_flagged_reference_union <- uniqueN(retained_published_flags$sample_id)

n_critical_metadata_issues <- nrow(metadata_issues[severity == "critical"])
n_high_metadata_issues <- nrow(metadata_issues[severity == "high"])
n_documentation_issues <- nrow(metadata_issues[severity == "documentation"])

n_single_collection_classes <- nrow(
  study_class_summary[single_collection_class %in% TRUE]
)

n_collection_holdout_evaluable <- nrow(
  collection_holdout_folds[fold_evaluable %in% TRUE]
)

n_collection_holdout_total <- nrow(collection_holdout_folds)

summary_lines <- c(
  "PESTFLY REFERENCE SET AUDIT",
  "",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Repository root: ", repo_root),
  paste0("Metadata source: ", metadata_path),
  paste0("Duplicate profile panel: ", opt$duplicate_panel),
  "",
  "FROZEN SNAPSHOT",
  paste0("Total samples: ", nrow(meta)),
  paste0("References: ", nrow(refs)),
  paste0("Queries: ", nrow(queries)),
  paste0(
    "Reference macroregions: ",
    paste(
      reference_class_counts[, paste0(macroregion_3, "=", sum(N)), by = macroregion_3]$V1,
      collapse = "; "
    )
  ),
  "",
  "IDENTIFIER AND METADATA CHECKS",
  paste0("Duplicated sample identifiers: ", nrow(duplicate_sample_ids)),
  paste0("Duplicated reference accessions: ", nrow(duplicate_accessions)),
  paste0("Critical metadata issues: ", n_critical_metadata_issues),
  paste0("High metadata issues: ", n_high_metadata_issues),
  paste0("Documentation completeness issues: ", n_documentation_issues),
  paste0(
    "References below ortholog completeness threshold ",
    format(opt$min_fraction_ogs_seen, scientific = FALSE),
    ": ",
    sum(reference_ortholog_coverage$below_completeness_threshold, na.rm = TRUE)
  ),
  "",
  "PUBLISHED FLAGS RETAINED IN THE REFERENCE SET",
  paste0("Nuclear PCA outliers: ", n_pca_outlier_refs),
  paste0("Published hybrid candidates: ", n_hybrid_candidate_refs),
  paste0("Additional COI taxonomic caution specimens: ", n_coi_caution_refs),
  paste0("Unique retained references with at least one published flag: ", n_flagged_reference_union),
  "",
  "GENOTYPE DUPLICATE CHECK",
  paste0("P1 SNPs examined: ", nrow(X_ref)),
  paste0("Reference profiles examined: ", ncol(X_ref)),
  paste0("Minimum compared SNPs per reported pair: ", opt$min_pairwise_snps),
  paste0(
    "Near duplicate threshold: ",
    format(opt$near_duplicate_concordance, digits = 6)
  ),
  paste0("Near duplicate pairs: ", nrow(near_duplicate_pairs)),
  paste0("Identical on overlap pairs: ", nrow(exact_overlap_pairs)),
  paste0("Fully identical profiles including missingness: ", nrow(full_profile_identical_pairs)),
  "",
  "STUDY SUPPORT",
  paste0("Panel classes represented by only one source collection: ", n_single_collection_classes),
  paste0(
    "Evaluable leave one collection out folds by class support: ",
    n_collection_holdout_evaluable,
    "/",
    n_collection_holdout_total
  ),
  "",
  "PREDECLARED SENSITIVITY TIERS",
  paste0(
    "Published hybrid candidate exclusion: ",
    length(scenario_sets$exclude_published_hybrid_candidates),
    " references"
  ),
  paste0(
    "All nuclear PCA outlier exclusion: ",
    length(scenario_sets$exclude_all_nuclear_PCA_outliers),
    " references"
  ),
  paste0(
    "Extended taxonomic caution exclusion: ",
    length(scenario_sets$exclude_extended_taxonomic_caution),
    " references"
  ),
  "",
  "INTERPRETATION RULE",
  paste(
    "No specimen was removed by this audit.",
    "Definite metadata errors require correction before later analyses.",
    "Published biological cautions are retained in the primary analysis and",
    "tested only through the predeclared sensitivity tiers."
  )
)

writeLines(
  summary_lines,
  con = file.path(opt$out_dir, "REFERENCE_SET_AUDIT_SUMMARY.txt"),
  useBytes = TRUE
)

run_info <- list(
  analysis = "validation_03_reference_set_audit",
  repo_root = repo_root,
  metadata_path = metadata_path,
  qc_dir = opt$qc_dir,
  step2_dir = opt$step2_dir,
  out_dir = opt$out_dir,
  duplicate_panel = opt$duplicate_panel,
  parameters = list(
    min_pairwise_snps = opt$min_pairwise_snps,
    near_duplicate_concordance = opt$near_duplicate_concordance,
    top_pair_rows = opt$top_pair_rows,
    min_fraction_ogs_seen = opt$min_fraction_ogs_seen,
    min_ref_per_group = opt$min_ref_per_group
  ),
  counts = list(
    n_samples = nrow(meta),
    n_references = nrow(refs),
    n_queries = nrow(queries),
    n_pca_outlier_references = n_pca_outlier_refs,
    n_hybrid_candidate_references = n_hybrid_candidate_refs,
    n_coi_caution_references = n_coi_caution_refs,
    n_near_duplicate_pairs = nrow(near_duplicate_pairs),
    n_full_profile_identical_pairs = nrow(full_profile_identical_pairs)
  ),
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(run_info, file.path(opt$out_dir, "reference_set_audit_run_info.rds"))

message("")
message("Reference set audit complete.")
message("Summary: ", file.path(opt$out_dir, "REFERENCE_SET_AUDIT_SUMMARY.txt"))
message("Sample audit: ", file.path(opt$out_dir, "reference_sample_audit.tsv"))
message("No reference samples were removed or relabelled.")
