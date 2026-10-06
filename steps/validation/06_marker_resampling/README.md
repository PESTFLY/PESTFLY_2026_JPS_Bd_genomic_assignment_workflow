# Alternative marker set resampling

For each K, partition the top 2K corrected Hudson ranked markers into complementary rank matched panels of K markers, with thirty alternative panels from fifteen partitions.

## Scope and interpretation

Paired A and B sets have no marker overlap and share a rank distribution. Different partitions can overlap, so the thirty sets are not thirty independent discovery datasets. This tests sensitivity to nonnested marker choice conditional on the existing Hudson ranking.

Base seed is 20261001. Default K values are 500, 1000, 2000 and 5000. Fixed size panels retain the likelihood, posterior, gap and minimum SNP requirements, but do not apply within panel multi K convergence criteria. Their stability is assessed across alternative marker sets.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels/
results/04_origin_assignment/final_assignment.tsv
results/05_loo_validation/
```

## Outputs

Archived output directory: `results/validation/06_marker_resampling/`.

The analysis retains parameter or design records, summary tables, detailed specimen or replicate results, figures where produced and a run record. Scenario analyses also retain outputs in `scenarios/`. Existing `return_bundle/` directories contain selected exports, not every full output.

## Run

```bash
Rscript steps/validation/06_marker_resampling/run.R
Rscript steps/validation/06_marker_resampling/run.R --help
```

Use the repository root as working directory. Scenario analyses 04 and 07 require regenerated FASTA alignments; the other analyses use archived QC and panel resources. A rerun may replace outputs or reuse already validated scenario tables as recorded by the script.

## Parameters

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Directory containing metadata_clean.tsv [default %default] |
| `--step2_dir` | `results/02_snp_panels` | Corrected Step 02 panel directory [default %default] |
| `--assignment_dir` | `results/04_origin_assignment` | Corrected Step 04 assignment directory [default %default] |
| `--loo_dir` | `results/05_loo_validation` | Corrected individual LOO directory [default %default] |
| `--out_dir` | `results/validation/06_marker_resampling` | Supplementary validation output directory [default %default] |
| `--panels` | `paste( c( "P1_macroregion_africa_vs_asia", "P2_subregion_within_africa", "P3_subregion_within_asia" ), collapse = "," )` | Comma separated panel identifiers [default %default] |
| `--panel_sizes` | `500,1000,2000,5000` | Comma separated fixed marker panel sizes [default %default] |
| `--n_partitions` | `15L` | Complementary A/B partitions per panel and K [default %default] |
| `--n_cores` | `0L` | Parallel workers; zero uses up to eight detected cores minus one [default %default] |
| `--seed` | `20261001L` | Base random seed [default %default] |
| `--pseudocount` | `0.5` | Allele frequency pseudocount [default %default] |
| `--epsilon` | `0.02` | Per site allele flip probability [default %default] |
| `--min_snps_query_macroregion` | `200L` | Minimum usable SNPs for P1 [default %default] |
| `--min_snps_query_subregion` | `300L` | Minimum usable SNPs for P2 and P3 [default %default] |
| `--high_posterior` | `0.95` | High posterior threshold [default %default] |
| `--moderate_posterior` | `0.85` | Moderate posterior threshold [default %default] |
| `--min_gap_macroregion` | `0.20` | Minimum P1 posterior gap [default %default] |
| `--min_gap_subregion` | `0.25` | Minimum P2 and P3 posterior gap [default %default] |

Required packages: `data.table`, `optparse`, `parallel`.

These script defaults are distinct from archived execution records. For reproducibility retain those records and fixed seeds. See [conceptual guide](../../../docs/CONCEPTS.md), [installation](../../../docs/INSTALL.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md).
