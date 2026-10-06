# Individual fixed panel leave one out validation

Validate reference prediction with each focal specimen excluded from class allele frequency estimation, using the fixed SNP panels selected in Step 02.

## Conceptual interpretation

The focal specimen contributed to the initial SNP discovery and ranking. This validates prediction with fixed panels, not a fully nested marker discovery procedure. Use the geographically grouped analysis for a stricter assessment of reference geography.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels/<panel_id>/snp_map.tsv
results/02_snp_panels/panels/<panel_id>/snp_matrix.rds
```

## Main outputs

```text
results/05_loo_validation/step4b_loo_summary_all_panels.tsv
results/05_loo_validation/step4b_loo_predictions_all_panels.tsv
results/05_loo_validation/step4b_loo_raw_all_panels.tsv.gz
results/05_loo_validation/step4b_run_info.rds
```

## Run

Run from the repository root:

```bash
Rscript steps/05_full_panel_loo_validation/run.R
Rscript steps/05_full_panel_loo_validation/run.R --help
```

## Parameters

These are script defaults, transcribed from the actual optparse definitions. Archived parameter and run records describe the benchmark execution.

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 output directory containing metadata_clean.tsv [default %default] |
| `--step3_dir` | `results/02_snp_panels` | Step 02 output directory containing panels/ [default %default] |
| `--out_dir` | `results/05_loo_validation` | Step 05 output directory [default %default] |
| `--panels` | `P1_macroregion_africa_vs_asia,P2_subregion_within_africa,P3_subregion_within_asia` | Comma-separated panel IDs to validate [default %default] |
| `--min_snps_query_macroregion` | `200` | Minimum non-missing SNPs required for P1 macroregion LOO assignment [default %default] |
| `--min_snps_query_subregion` | `300` | Minimum non-missing SNPs required for P2/P3 subregion LOO assignment [default %default] |
| `--min_ref_per_group` | `3` | Minimum reference samples required per class before LOO [default %default] |
| `--keep_snps_frac` | `0.70` | Keep SNPs with at least this non-missing fraction across panel references [default %default] |
| `--pseudocount` | `0.5` | Laplace smoothing pseudocount for reference allele frequencies [default %default] |
| `--epsilon` | `0.02` | Per-site genotype error/flip probability [default %default] |
| `--high_posterior` | `0.95` | High posterior threshold [default %default] |
| `--moderate_posterior` | `0.85` | Moderate posterior lower bound [default %default] |
| `--min_gap_macroregion` | `0.20` | Minimum top-second posterior gap for P1 [default %default] |
| `--min_gap_subregion` | `0.25` | Minimum top-second posterior gap for P2/P3 [default %default] |
| `--base_K` | `1,2,3,4,5,10,20,50,100,200,500` | Comma-separated base K values [default %default] |
| `--extra_K` | `1000,2000,5000,10000` | Comma-separated extra K values tested when available [default %default] |
| `--auto_extend_K` | `TRUE` | Auto-extend K grid beyond base_K when enough SNPs are available: TRUE/FALSE [default %default] |
| `--include_all_snps_K` | `TRUE` | Include all available SNPs as final K value: TRUE/FALSE [default %default] |
| `--max_K_points` | `30` | Maximum number of K points after auto-extension [default %default] |
| `--min_agreement` | `0.90` | Minimum fraction of usable K calls matching final largest-K call [default %default] |
| `--min_K_available` | `6` | Minimum number of usable K results needed to judge stability [default %default] |
| `--tail_fraction` | `0.50` | Fraction of largest usable K values used for tail-stability check [default %default] |
| `--min_tail_agreement` | `1.00` | Required agreement with final call among largest usable-K tail values [default %default] |

## Technical notes

The likelihood, query SNP requirements, posterior, gap and multi K stability settings mirror Step 04. Raw top class accuracy and the fraction of reportable calls are different quantities; report both when assessing operational performance.

Required packages: `data.table`, `openxlsx`, `optparse`.

The script records parameters and session information with its outputs. Rerunning with defaults may replace archived files. Preserve a separate archive copy for comparison.

See [conceptual guide](../../docs/CONCEPTS.md), [installation](../../docs/INSTALL.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md).
