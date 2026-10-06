# Reports for plant health authorities

Translate saved Step 04 assignments into explicit six criterion reports. PASS requires all six criteria at that resolution. A macroregion PASS with subregion FAIL reports macroregion affinity only. Macroregion FAIL reports an unresolved assignment and the subregion is NOT_EVALUATED.

## Inputs

The required files are `results/04_origin_assignment/final_assignment.tsv` and its associated `step4_params.tsv`. The formatter reads thresholds from that parameter record rather than assuming the manuscript defaults apply to every run. It verifies that reconstructed criterion eligibility and posterior categories agree with Step 04. A contradictory input stops generation.

Optional inputs are RF tables in `results/06_rf_corroboration/`, validation summaries in `results/05_loo_validation/`, and a declared commodity country map. Missing RF panels remain unavailable and do not affect the six criterion decision.

## Outputs

| File in `results/07_authority_report/` | Contents |
| :--- | :--- |
| `FINAL_origin_tracing_authority_REPORT.xlsx` | Simplified and extended reports, Criteria and Reporting_Thresholds sheets, interpretation and retained reference and validation sheets |
| `FINAL_origin_tracing_authority_simplified.tsv` | One row per query, finest supported class, status and criterion counts at each resolution, every unmet criterion, separate RF and commodity comparisons |
| `FINAL_origin_tracing_authority_extended.tsv` | Saved Step 04 metrics, joined corroboration and metadata, criterion summaries |
| `FINAL_origin_tracing_authority_criteria.tsv` | Six rows per sample and resolution, with observed value, threshold, result and reason |
| `reporting_thresholds.tsv` | Thresholds and source field names for both resolutions |

The R formatter also writes current simplified and extended RDS objects and `step7_run_info.rds`. Its run record contains the threshold source, reporting parameters, criterion counts and session information.

## Run and check

Run from the repository root. Reporting requires no sequence alignments and does not refit any model.

```bash
Rscript steps/07_authority_facing_report/check_reporting.R
Rscript steps/07_authority_facing_report/run.R --out_dir results/07_authority_report_rerun
```

The first command uses base R to check the supplied benchmark, inclusive threshold boundaries, twelve single criterion failures, conditional branches, missing metrics and invalid input. The second writes to a new directory for comparison. Once checked, the default command can regenerate the standard report location:

```bash
Rscript steps/07_authority_facing_report/run.R
```

## Checkpoint provenance

The current workbook, four TSV tables, simplified and extended RDS objects and `step7_run_info.rds` are the byte unchanged native R outputs from 5 October 2026. R 4.5.1 used data.table 1.18.0, optparse 1.7.5 and openxlsx 4.2.8.1. The isolated implementation and reporting checks passed, as did the formatter itself. All 22 reported classes, resolutions and confidence categories agree with Step 04 and the independent checkpoint diagnostic.

The earlier reporting reconstruction is preserved in `reconstruction_checkpoint/`. Original colour workbook, TSV tables, RDS objects and historical run record remain unchanged in `legacy_colour_report/`.

The original release checker stopped on 25 descriptive label differences in the extended table: “Step 4” versus “Step 04”. The same labels appear in the simplified table. `verify_release.R` now normalises that exact spelling only in `rf_corroboration_detail` and `commodity_origin_mismatch_detail`. The original FAIL log is preserved. Independent review completed the workbook and RDS comparisons and checked all 1,214 archived result file checksums. See [native R verification](../../docs/NATIVE_R_VERIFICATION.md).

The corrected release checker can be run from the repository root with `Rscript verify_release.R`. It writes a new `release_check/` directory and leaves archived results unchanged. The corrected wrapper has been syntax screened; its original component checks and formatter were executed in the supplied R session.

## Parameters

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--step4_final` | `results/04_origin_assignment/final_assignment.tsv` | Saved Step 04 assignment table |
| `--step4_params` | Empty | Uses `step4_params.tsv` beside the assignment file. Supply explicitly for another parameter record |
| `--step5_dir` | `results/06_rf_corroboration` | Optional RF output directory |
| `--step4b_dir` | `results/05_loo_validation` | Optional individual LOO output directory |
| `--step6b_dir` | Empty | Optional internal reduced panel diagnostic directory |
| `--intercept_origin_file` | Empty | Uses the supplied built in intercept map. An explicitly supplied path must exist and have unique populated sample identifiers |
| `--out_dir` | `results/07_authority_report` | Output directory |
| `--out_xlsx` | `FINAL_origin_tracing_authority_REPORT.xlsx` | Workbook filename |
| `--include_validation_sheets` | `TRUE` | Include available validation sheets |

Required packages for `run.R` are `data.table`, `openxlsx` and `optparse`. `reporting_criteria.R` and `check_reporting.R` use base R. RF remains a separate algorithm comparison. Declared commodity country is a consistency check and does not establish the biological source or transport route.

See [conceptual guide](../../docs/CONCEPTS.md), [installation](../../docs/INSTALL.md), and [reporting criteria](../../docs/REPORTING_CRITERIA.md).
