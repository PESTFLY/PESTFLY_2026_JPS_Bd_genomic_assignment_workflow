# Random Forest corroboration

Train a second classification algorithm on the same fixed SNP resources and compare its predictions with the primary hierarchical assignments.

## Conceptual interpretation

Independence here refers to the classifier, not to the underlying data or feature discovery. Cross validation uses fixed Step 02 SNP resources. RF agrees with, qualifies or contradicts the primary call but does not override it or contribute a seventh reporting criterion.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels/<panel_id>/snp_matrix.rds
results/04_origin_assignment/final_assignment.tsv, optional comparison
```

## Main outputs

```text
results/06_rf_corroboration/RF_all_panels_cv_summary.tsv
results/06_rf_corroboration/RF_all_panels_predictions_queries.tsv
results/06_rf_corroboration/step5_run_info.rds
```

## Run

Run from the repository root:

```bash
Rscript steps/06_random_forest_corroboration/run.R
Rscript steps/06_random_forest_corroboration/run.R --help
```

## Parameters

These are script defaults, transcribed from the actual optparse definitions. Archived parameter and run records describe the benchmark execution.

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 output directory containing metadata_clean.tsv [default %default] |
| `--step3_dir` | `results/02_snp_panels` | Step 02 output directory containing panels/ [default %default] |
| `--step4_dir` | `results/04_origin_assignment` | Step 04 output directory containing final_assignment.tsv, if available [default %default] |
| `--out_dir` | `results/06_rf_corroboration` | Step 06 output directory [default %default] |
| `--panels` | `P1_macroregion_africa_vs_asia,P2_subregion_within_africa,P3_subregion_within_asia` | Comma-separated Step 02 panel IDs to run [default %default] |
| `--query_is_reference_value` | `FALSE` | Metadata is_reference value identifying query/intercept samples [default %default] |
| `--keep_snps_frac` | `0.70` | Keep SNPs with at least this non-missing fraction across panel samples [default %default] |
| `--min_snps_ref` | `500` | Drop reference samples with fewer than this many non-missing SNPs [default %default] |
| `--min_snps_query` | `200` | Flag query samples with fewer than this many non-missing SNPs as uncertain [default %default] |
| `--min_class_n` | `5` | Drop reference classes with fewer than this many samples [default %default] |
| `--trees` | `1000` | Number of trees [default %default] |
| `--mtry` | `0` | mtry; 0 = sqrt(number of retained SNPs) [default %default] |
| `--min_node_size` | `1` | Minimum node size for ranger [default %default] |
| `--seed` | `1` | Random seed [default %default] |
| `--cores` | `0` | Threads for ranger. 0 = all detected cores minus one [default %default] |
| `--folds` | `5` | Number of stratified CV folds [default %default] |
| `--reps` | `5` | Number of repeated CV rounds [default %default] |
| `--write_cv_rows` | `TRUE` | Write per-sample CV rows: TRUE/FALSE [default %default] |
| `--high_prob` | `0.95` | High-confidence RF threshold [default %default] |
| `--moderate_prob` | `0.85` | Moderate-confidence RF threshold [default %default] |
| `--min_gap` | `0.20` | Flag uncertain if top probability minus second probability is below this [default %default] |

## Technical notes

The archived run used ranger, 1000 trees, fivefold cross validation repeated five times and base seed 1. The script prefers ranger when available and otherwise uses randomForest. Record the backend when reproducing results. Training fold means impute missing values for training and test matrices; held out samples do not supply those imputation means.

Required packages: `data.table`, `openxlsx`, `optparse`, `randomForest`, `ranger`.

The script records parameters and session information with its outputs. Rerunning with defaults may replace archived files. Preserve a separate archive copy for comparison.

See [conceptual guide](../../docs/CONCEPTS.md), [installation](../../docs/INSTALL.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md).
