# Hierarchical genomic assignment

Assign each query to a macroregion, then evaluate its subregion when macroregion reporting criteria pass.

K is the number of highest ranked SNPs. The likelihood uses equal class priors, pseudocount 0.5 and allele flip allowance 0.02. See the shared reporting criteria for convergence and reportability.

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

From the repository root:

```bash
Rscript steps/04_hierarchical_origin_assignment/run.R
Rscript steps/04_hierarchical_origin_assignment/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 QC directory |
| `--step3_dir` | `results/02_snp_panels` | Step 02 SNP directory |
| `--out_dir` | `results/04_origin_assignment` | Output directory |
| `--out_xlsx` | `step4_assignment_multilevel_macroregion_first.xlsx` | Workbook filename |
| `--query_is_reference_value` | `FALSE` | Metadata value identifying queries |
| `--min_snps_query_macroregion` | `200` | Minimum usable SNPs for P1 |
| `--min_snps_query_subregion` | `300` | Minimum usable SNPs for P2 and P3 |
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
