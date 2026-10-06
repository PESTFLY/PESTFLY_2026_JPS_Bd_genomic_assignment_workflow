# Reports for plant health authorities

Format saved Step 04 assignments into reports showing all six criterion outcomes and the finest supported resolution. RF and commodity comparisons remain separate.

## Inputs and outputs

Required inputs are `final_assignment.tsv` and its matching `step4_params.tsv`. The formatter reads their thresholds and checks consistency with Step 04. RF tables, individual validation summaries and a declared commodity country map are optional.

Files in `results/07_authority_report/`:

| File | Contents |
| :--- | :--- |
| `FINAL_origin_tracing_authority_REPORT.xlsx` | Reports, criterion values, thresholds and supporting sheets |
| `FINAL_origin_tracing_authority_simplified.tsv` | One row per query, supported class, criterion counts and unmet criteria |
| `FINAL_origin_tracing_authority_extended.tsv` | Assignment metrics and joined supporting information |
| `FINAL_origin_tracing_authority_criteria.tsv` | Six criterion rows per specimen and resolution |
| `reporting_thresholds.tsv` | Thresholds and source fields |

The script also writes simplified and extended RDS tables and `step7_run_info.rds`.

## Run and check

From the repository root, check the rules and write reports to a separate directory:

```bash
Rscript steps/07_authority_facing_report/check_reporting.R
Rscript steps/07_authority_facing_report/run.R --out_dir results/07_authority_report_rerun
```

No sequence inputs or model fitting are needed. The default output directory contains the archived native R exports from 5 October 2026; use a different directory for comparison. Earlier reports remain in `legacy_colour_report/` and `reconstruction_checkpoint/`. See [native R verification](../../docs/NATIVE_R_VERIFICATION.md) for the original checker mismatch and completed export review.

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

Packages: `data.table`, `openxlsx` and `optparse`. The reporting rules and their checks use base R.

See [reporting criteria](../../docs/REPORTING_CRITERIA.md), [concepts](../../docs/CONCEPTS.md) and [installation](../../docs/INSTALL.md). The manuscript provides the full interpretation.
