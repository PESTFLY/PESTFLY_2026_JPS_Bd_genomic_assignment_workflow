# Grouped geographic validation

Withhold a country or a collection site within a country and predict the held out reference specimens.

Each fold filters and reranks retained candidate SNPs and estimates frequencies using training references. It retains the original one SNP per ortholog candidate set. Folds without adequate true class support are non evaluable; report their counts alongside accuracy.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels/
```

## Outputs

Archived output directory: `results/validation/02_grouped_geographic_cv/`.

## Run

From the repository root:

```bash
Rscript steps/validation/02_grouped_geographic_cv/run.R
Rscript steps/validation/02_grouped_geographic_cv/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 QC directory |
| `--step2_dir` | `results/02_snp_panels` | Step 02 SNP directory |
| `--loo_dir` | `results/05_loo_validation` | Step 05 validation directory |
| `--out_dir` | `results/validation/02_grouped_geographic_cv` | Output directory |
| `--panels` | `P1_macroregion_africa_vs_asia,P2_subregion_within_africa,P3_subregion_within_asia` | Comma separated panel identifiers |
| `--fold_units` | `country,site` | Country or site holdouts, or both |
| `--n_cores` | `0` | Workers; 0 uses detected cores minus one, capped at eight |
| `--min_ref_per_group` | `3` | Minimum references per class |
| `--max_site_missing` | `0.20` | Maximum reference missing fraction per SNP |
| `--min_mac` | `2` | Minimum reference minor allele count |
| `--foldwise_reranking` | `TRUE` | Recalculate Hudson scores and ranking in each fold |
| `--pseudocount` | `0.5` | Reference allele frequency pseudocount |
| `--epsilon` | `0.02` | Per site allele flip allowance |
| `--min_snps_query_macroregion` | `200` | Minimum usable SNPs for P1 |
| `--min_snps_query_subregion` | `300` | Minimum usable SNPs for P2 and P3 |
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
| `--top_markers_per_fold` | `20` | Fold ranked markers retained for inspection |

Required packages: `data.table`, `openxlsx`, `optparse`, `parallel`.

See [installation](../../../docs/INSTALL.md), [concepts](../../../docs/CONCEPTS.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md). The manuscript provides the full methods and interpretation.
