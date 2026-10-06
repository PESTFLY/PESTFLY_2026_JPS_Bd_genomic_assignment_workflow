# Balanced reference downsampling

Repeat balanced sampling of class references at multiple training sizes, using independently held out reference test sets and thirty replicates per design.

## Scope and interpretation

Training only filtering, Hudson ranking and allele frequency estimation are repeated within the retained candidate SNP set. The analysis does not rediscover alternative ortholog sites. Raw accuracy and uncertainty rate should be read together; high accuracy among calls does not imply that all specimens were reportable.

Base seed is 20260930. The seed rule adds panel index times 1000000, training size index times 10000 and replicate. Default requested sizes are 3, 5, 10, 20, 40, 80 and 100 references per class, plus each panel maximum evaluable balanced size. Class size constraints limit the actual designs, recorded in downsampling_design.tsv.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels/
results/04_origin_assignment/final_assignment.tsv
results/05_loo_validation/
```

## Outputs

Archived output directory: `results/validation/05_reference_downsampling/`.

The analysis retains parameter or design records, summary tables, detailed specimen or replicate results, figures where produced and a run record. Scenario analyses also retain outputs in `scenarios/`. Existing `return_bundle/` directories contain selected exports, not every full output.

## Run

```bash
Rscript steps/validation/05_reference_downsampling/run.R
Rscript steps/validation/05_reference_downsampling/run.R --help
```

Use the repository root as working directory. Scenario analyses 04 and 07 require regenerated FASTA alignments; the other analyses use archived QC and panel resources. A rerun may replace outputs or reuse already validated scenario tables as recorded by the script.

## Parameters

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Directory containing metadata_clean.tsv [default %default] |
| `--step2_dir` | `results/02_snp_panels` | Corrected Step 02 panel directory [default %default] |
| `--assignment_dir` | `results/04_origin_assignment` | Corrected Step 04 assignment directory [default %default] |
| `--loo_dir` | `results/05_loo_validation` | Corrected individual LOO directory [default %default] |
| `--out_dir` | `results/validation/05_reference_downsampling` | Supplementary validation output directory [default %default] |
| `--panels` | `paste( c( "P1_macroregion_africa_vs_asia", "P2_subregion_within_africa", "P3_subregion_within_asia" ), collapse = "," )` | Comma separated panel identifiers [default %default] |
| `--n_replicates` | `30L` | Random replicates per panel and training size [default %default] |
| `--n_cores` | `0L` | Parallel workers; zero uses up to eight detected cores minus one [default %default] |
| `--seed` | `20260930L` | Base random seed [default %default] |
| `--training_sizes` | `3,5,10,20,40,80,100` | Candidate references per class; each panel also includes its maximum evaluable balanced size [default %default] |
| `--test_per_group` | `5L` | Maximum balanced held out references per class and replicate [default %default] |
| `--min_ref_per_group` | `3L` | Minimum training references per class [default %default] |
| `--max_site_missing` | `0.20` | Maximum training missing fraction for a retained SNP [default %default] |
| `--min_mac` | `2L` | Minimum training minor allele count [default %default] |
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
| `--min_agreement` | `0.90` | Minimum global agreement across usable K values [default %default] |
| `--min_K_available` | `6L` | Minimum usable K values for a reliable call [default %default] |
| `--tail_fraction` | `0.50` | Largest K fraction used for tail stability [default %default] |
| `--min_tail_agreement` | `1.00` | Required agreement in the largest K tail [default %default] |
| `--marker_overlap_K` | `20,100,500,2000,5000` | K values for overlap with the complete Step 02 ranking [default %default] |

Required packages: `data.table`, `optparse`, `parallel`.

These script defaults are distinct from archived execution records. For reproducibility retain those records and fixed seeds. See [conceptual guide](../../../docs/CONCEPTS.md), [installation](../../../docs/INSTALL.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md).
