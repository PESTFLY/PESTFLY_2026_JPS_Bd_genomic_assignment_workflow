# Grouped geographic validation

Withhold all reference specimens from each country or country nested collection site, then evaluate predictions for the held out geographic unit.

## Scope and interpretation

Each evaluable fold filters the retained SNP set, recalculates corrected Hudson scores, reranks markers and estimates allele frequencies using training references only. The retained one SNP per ortholog candidates were discovered beforehand. Alternative sites are not rediscovered from alignments, so this is not fully nested de novo feature discovery.

A fold is not evaluable if the true class is absent after withholding, or training representation fails required support. The benchmark Central Africa country and site folds and East Asia country fold illustrate missing class coverage. Report evaluable sample counts alongside accuracy.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels/
```

## Outputs

Archived output directory: `results/validation/02_grouped_geographic_cv/`.

The analysis retains parameter or design records, summary tables, detailed specimen or replicate results, figures where produced and a run record. Scenario analyses also retain outputs in `scenarios/`. Existing `return_bundle/` directories contain selected exports, not every full output.

## Run

```bash
Rscript steps/validation/02_grouped_geographic_cv/run.R
Rscript steps/validation/02_grouped_geographic_cv/run.R --help
```

Use the repository root as working directory. Scenario analyses 04 and 07 require regenerated FASTA alignments; the other analyses use archived QC and panel resources. A rerun may replace outputs or reuse already validated scenario tables as recorded by the script.

## Parameters

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Directory containing metadata_clean.tsv [default %default] |
| `--step2_dir` | `results/02_snp_panels` | Corrected Step 02 panel directory [default %default] |
| `--loo_dir` | `results/05_loo_validation` | Existing individual LOO directory used for comparison [default %default] |
| `--out_dir` | `results/validation/02_grouped_geographic_cv` | Separate grouped validation output directory [default %default] |
| `--panels` | `paste( c( "P1_macroregion_africa_vs_asia", "P2_subregion_within_africa", "P3_subregion_within_asia" ), collapse = "," )` | Comma separated panel identifiers [default %default] |
| `--fold_units` | `country,site` | Grouped validation units: country, site, or both [default %default] |
| `--n_cores` | `0L` | Parallel workers. Zero uses up to eight detected cores minus one [default %default] |
| `--min_ref_per_group` | `3L` | Minimum training references required for a candidate class [default %default] |
| `--max_site_missing` | `0.20` | Maximum training missing fraction for a retained SNP [default %default] |
| `--min_mac` | `2L` | Minimum training minor allele count for a retained SNP [default %default] |
| `--foldwise_reranking` | `TRUE` | Recalculate Hudson FST and rerank retained SNPs per fold [default %default] |
| `--pseudocount` | `0.5` | Allele frequency pseudocount [default %default] |
| `--epsilon` | `0.02` | Per site allele flip probability [default %default] |
| `--min_snps_query_macroregion` | `200L` | Minimum usable SNPs for P1 [default %default] |
| `--min_snps_query_subregion` | `300L` | Minimum usable SNPs for P2 and P3 [default %default] |
| `--high_posterior` | `0.95` | High posterior threshold [default %default] |
| `--moderate_posterior` | `0.85` | Moderate posterior threshold [default %default] |
| `--min_gap_macroregion` | `0.20` | Minimum P1 posterior gap [default %default] |
| `--min_gap_subregion` | `0.25` | Minimum P2 and P3 posterior gap [default %default] |
| `--base_K` | `1,2,3,4,5,10,20,50,100,200,500` | Comma separated base K values [default %default] |
| `--extra_K` | `1000,2000,5000,10000` | Additional K values when available [default %default] |
| `--auto_extend_K` | `TRUE` | Use extra K values when available [default %default] |
| `--include_all_snps_K` | `TRUE` | Include all eligible SNPs as the final K [default %default] |
| `--max_K_points` | `30L` | Maximum number of K values [default %default] |
| `--min_agreement` | `0.90` | Minimum agreement across usable K values [default %default] |
| `--min_K_available` | `6L` | Minimum number of usable K values [default %default] |
| `--tail_fraction` | `0.50` | Largest K fraction used for tail stability [default %default] |
| `--min_tail_agreement` | `1.00` | Required agreement in the largest K tail [default %default] |
| `--top_markers_per_fold` | `20L` | Number of fold ranked markers retained for audit [default %default] |

Required packages: `data.table`, `openxlsx`, `optparse`, `parallel`.

These script defaults are distinct from archived execution records. For reproducibility retain those records and fixed seeds. See [conceptual guide](../../../docs/CONCEPTS.md), [installation](../../../docs/INSTALL.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md).
