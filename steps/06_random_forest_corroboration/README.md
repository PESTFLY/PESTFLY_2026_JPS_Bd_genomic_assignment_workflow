# Random Forest corroboration

Fit a second classifier to the fixed SNP resources and compare query predictions with Step 04 assignments.

The benchmark used ranger, 1000 trees, fivefold validation repeated five times and seed 1. Training fold means supply missing value imputation. The script falls back to randomForest if ranger is unavailable; record the backend used. RF remains separate from the six reporting criteria.

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

From the repository root:

```bash
Rscript steps/06_random_forest_corroboration/run.R
Rscript steps/06_random_forest_corroboration/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 QC directory |
| `--step3_dir` | `results/02_snp_panels` | Step 02 SNP directory |
| `--step4_dir` | `results/04_origin_assignment` | Step 04 assignment directory |
| `--out_dir` | `results/06_rf_corroboration` | Output directory |
| `--panels` | `P1_macroregion_africa_vs_asia,P2_subregion_within_africa,P3_subregion_within_asia` | Comma separated panel identifiers |
| `--query_is_reference_value` | `FALSE` | Metadata value identifying queries |
| `--keep_snps_frac` | `0.70` | Keep SNPs with at least this non-missing fraction across panel samples |
| `--min_snps_ref` | `500` | Drop reference samples with fewer than this many non-missing SNPs |
| `--min_snps_query` | `200` | Flag query samples with fewer than this many non-missing SNPs as uncertain |
| `--min_class_n` | `5` | Minimum reference specimens per class |
| `--trees` | `1000` | Number of trees |
| `--mtry` | `0` | mtry; 0 = sqrt(number of retained SNPs) |
| `--min_node_size` | `1` | Minimum node size for ranger |
| `--seed` | `1` | Base random seed |
| `--cores` | `0` | RF threads; 0 uses detected cores minus one |
| `--folds` | `5` | Number of stratified CV folds |
| `--reps` | `5` | Number of repeated CV rounds |
| `--write_cv_rows` | `TRUE` | Write per-sample CV rows: TRUE/FALSE |
| `--high_prob` | `0.95` | High-confidence RF threshold |
| `--moderate_prob` | `0.85` | Moderate-confidence RF threshold |
| `--min_gap` | `0.20` | Flag uncertain if top probability minus second probability is below this |

Required packages: `data.table`, `openxlsx`, `optparse`, `randomForest`, `ranger`.

See [installation](../../docs/INSTALL.md), [concepts](../../docs/CONCEPTS.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md). The manuscript provides the full methods and interpretation.
