# Alternative marker set resampling

Partition the top 2K Hudson ranked markers into rank matched A and B sets of K markers, giving thirty panels from fifteen partitions.

A and B are disjoint within a partition; different partitions can overlap. Each fixed size panel uses SNP, posterior and gap requirements. Stability is assessed across alternative sets rather than through the cumulative K criteria.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels/
results/04_origin_assignment/final_assignment.tsv
results/05_loo_validation/
```

## Outputs

Archived output directory: `results/validation/06_marker_resampling/`.

## Run

From the repository root:

```bash
Rscript steps/validation/06_marker_resampling/run.R
Rscript steps/validation/06_marker_resampling/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 QC directory |
| `--step2_dir` | `results/02_snp_panels` | Step 02 SNP directory |
| `--assignment_dir` | `results/04_origin_assignment` | Step 04 assignment directory |
| `--loo_dir` | `results/05_loo_validation` | Step 05 validation directory |
| `--out_dir` | `results/validation/06_marker_resampling` | Output directory |
| `--panels` | `P1_macroregion_africa_vs_asia,P2_subregion_within_africa,P3_subregion_within_asia` | Comma separated panel identifiers |
| `--panel_sizes` | `500,1000,2000,5000` | Marker panel sizes |
| `--n_partitions` | `15` | Complementary A/B partitions per panel and K |
| `--n_cores` | `0` | Workers; 0 uses detected cores minus one, capped at eight |
| `--seed` | `20261001` | Base random seed |
| `--pseudocount` | `0.5` | Reference allele frequency pseudocount |
| `--epsilon` | `0.02` | Per site allele flip allowance |
| `--min_snps_query_macroregion` | `200` | Minimum usable SNPs for P1 |
| `--min_snps_query_subregion` | `300` | Minimum usable SNPs for P2 and P3 |
| `--high_posterior` | `0.95` | High posterior threshold |
| `--moderate_posterior` | `0.85` | Moderate posterior threshold |
| `--min_gap_macroregion` | `0.20` | Minimum P1 posterior gap |
| `--min_gap_subregion` | `0.25` | Minimum P2 and P3 posterior gap |

Required packages: `data.table`, `optparse`, `parallel`.

See [installation](../../../docs/INSTALL.md), [concepts](../../../docs/CONCEPTS.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md). The manuscript provides the full methods and interpretation.
