# Individual fixed panel leave one out validation

Withhold each reference from allele frequency estimation and predict its class using fixed Step 02 SNP panels.

The focal specimen contributed to the original marker discovery and ranking. Report accuracy together with the fraction of reportable calls. Use grouped geographic validation to assess country and site holdouts.

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

From the repository root:

```bash
Rscript steps/05_full_panel_loo_validation/run.R
Rscript steps/05_full_panel_loo_validation/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 QC directory |
| `--step3_dir` | `results/02_snp_panels` | Step 02 SNP directory |
| `--out_dir` | `results/05_loo_validation` | Output directory |
| `--panels` | `P1_macroregion_africa_vs_asia,P2_subregion_within_africa,P3_subregion_within_asia` | Comma separated panel identifiers |
| `--min_snps_query_macroregion` | `200` | Minimum usable SNPs for P1 |
| `--min_snps_query_subregion` | `300` | Minimum usable SNPs for P2 and P3 |
| `--min_ref_per_group` | `3` | Minimum references per class |
| `--keep_snps_frac` | `0.70` | Keep SNPs with at least this non-missing fraction across panel references |
| `--pseudocount` | `0.5` | Reference allele frequency pseudocount |
| `--epsilon` | `0.02` | Per site allele flip allowance |
| `--high_posterior` | `0.95` | High posterior threshold |
| `--moderate_posterior` | `0.85` | Moderate posterior threshold |
| `--min_gap_macroregion` | `0.20` | Minimum P1 posterior gap |
| `--min_gap_subregion` | `0.25` | Minimum P2 and P3 posterior gap |
| `--base_K` | `1,2,3,4,5,10,20,50,100,200,500` | Base K values |
| `--extra_K` | `1000,2000,5000,10000` | Additional K values when available |
| `--auto_extend_K` | `TRUE` | Include additional K values when available |
| `--include_all_snps_K` | `TRUE` | Add all eligible SNPs as final K |
| `--max_K_points` | `30` | Maximum number of K values |
| `--min_agreement` | `0.90` | Minimum agreement with final class across usable K |
| `--min_K_available` | `6` | Minimum usable K values |
| `--tail_fraction` | `0.50` | Largest K fraction used for tail agreement |
| `--min_tail_agreement` | `1.00` | Required tail agreement |

Required packages: `data.table`, `openxlsx`, `optparse`.

See [installation](../../docs/INSTALL.md), [concepts](../../docs/CONCEPTS.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md). The manuscript provides the full methods and interpretation.
