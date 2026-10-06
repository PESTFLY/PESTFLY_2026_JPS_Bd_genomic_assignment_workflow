# Hierarchical genomic assignment

Compare each query with represented macroregion classes, then evaluate a conditional subregion branch when the macroregion is reportable.

## Conceptual interpretation

K is the number of ranked SNPs included, not a population count or an ADMIXTURE parameter. Posterior support is conditional on represented classes and the likelihood model. High support does not show that the true source was sampled, and does not identify a transport route.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels_index.tsv
results/02_snp_panels/panels/<panel_id>/snp_map.tsv
results/02_snp_panels/panels/<panel_id>/snp_matrix.rds
```

## Main outputs

```text
results/04_origin_assignment/final_assignment.tsv
results/04_origin_assignment/P1_macroregion_raw.tsv
results/04_origin_assignment/P1_macroregion_stability.tsv
results/04_origin_assignment/step4_params.tsv
```

## Run

Run from the repository root:

```bash
Rscript steps/04_hierarchical_origin_assignment/run.R
Rscript steps/04_hierarchical_origin_assignment/run.R --help
```

## Parameters

These are script defaults, transcribed from the actual optparse definitions. Archived parameter and run records describe the benchmark execution.

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 output directory containing metadata_clean.tsv [default %default] |
| `--step3_dir` | `results/02_snp_panels` | Step 02 output directory containing panels_index.tsv and panels/ [default %default] |
| `--out_dir` | `results/04_origin_assignment` | Output directory [default %default] |
| `--out_xlsx` | `step4_assignment_multilevel_macroregion_first.xlsx` | Output Excel workbook filename [default %default] |
| `--query_is_reference_value` | `FALSE` | Metadata is_reference value identifying query/intercept samples [default %default] |
| `--min_snps_query_macroregion` | `200` | Minimum non-missing SNPs required for P1 macroregion assignment [default %default] |
| `--min_snps_query_subregion` | `300` | Minimum non-missing SNPs required for P2/P3 subregion assignment [default %default] |
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
| `--min_agreement` | `0.90` | Minimum fraction of usable K calls matching the final largest-K call [default %default] |
| `--min_K_available` | `6` | Minimum number of usable K results needed to judge stability [default %default] |
| `--tail_fraction` | `0.50` | Fraction of largest usable K values used for tail-stability check [default %default] |
| `--min_tail_agreement` | `1.00` | Required agreement with final call among largest usable-K tail values [default %default] |

## Technical notes

Smoothed consensus allele frequencies use pseudocount 0.5 and allele flip allowance epsilon 0.02, with equal class priors. The six reporting criteria and the cumulative K grid are documented in docs/REPORTING_CRITERIA.md. Thresholds are specified operational settings, not independently calibrated error rates.

Required packages: `data.table`, `openxlsx`, `optparse`.

The script records parameters and session information with its outputs. Rerunning with defaults may replace archived files. Preserve a separate archive copy for comparison.

See [conceptual guide](../../docs/CONCEPTS.md), [installation](../../docs/INSTALL.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md).
