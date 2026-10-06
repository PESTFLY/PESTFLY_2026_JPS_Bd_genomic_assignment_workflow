# SNP panel discovery and Hudson ranking

Discover biallelic consensus SNPs among eligible references, select one SNP per retained ortholog for each analytical panel and rank those SNPs using corrected Hudson FST scores.

## Conceptual interpretation

Reference and alternate consensus nucleotides are encoded 0 and 1. Ambiguous and IUPAC states are missing; these are not diploid genotype dosage calls. P1 compares Africa with Asia; P2 and P3 rank by the maximum pairwise subregion Hudson score. Query samples do not select or order markers.

## Inputs

```text
results/00_fasta/OG*.fasta
results/01_qc/metadata_clean.tsv
results/01_qc/ogs_pass_qc.txt
```

## Main outputs

```text
results/02_snp_panels/panels_index.tsv
results/02_snp_panels/panels/<panel_id>/snp_map.tsv
results/02_snp_panels/panels/<panel_id>/snp_matrix.rds
results/02_snp_panels/step3_run_info.rds
```

## Run

Run from the repository root:

```bash
Rscript steps/02_snp_panel_discovery_and_ranking/run.R
Rscript steps/02_snp_panel_discovery_and_ranking/run.R --help
```

## Parameters

These are script defaults, transcribed from the actual optparse definitions. Archived parameter and run records describe the benchmark execution.

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--align_dir` | `results/00_fasta` | Directory with Step 00 FASTA alignments [default %default] |
| `--file_glob` | `OG*.fasta` | FASTA file pattern [default %default] |
| `--qc_dir` | `results/01_qc` | Step 01 output directory [default %default] |
| `--out_dir` | `results/02_snp_panels` | Step 02 output directory [default %default] |
| `--reference_value` | `auto` | Which is_reference value means reference DB: auto, TRUE, or FALSE [default %default] |
| `--macroregions` | `Africa,Asia` | Comma-separated macroregion_3 values used for P1/P2/P3 [default %default] |
| `--max_site_missing` | `0.20` | Maximum missing fraction among reference samples at a site [default %default] |
| `--min_mac` | `2` | Minimum minor allele count among reference samples [default %default] |
| `--min_ref_per_group` | `3` | Minimum reference samples per class/group for a panel [default %default] |
| `--max_snps_per_og` | `1` | Maximum SNPs retained per OG per panel [default %default] |
| `--exclude_regex` | `^$` | Regex of taxa to exclude before SNP extraction, if needed [default %default] |
| `--allowed_extra_labels` | `Bdors,Blati` | Comma-separated labels allowed in FASTA but absent from metadata [default %default] |
| `--cores` | `0` | Parallel workers. 0 = all detected cores minus one [default %default] |
| `--chunk_size` | `300` | Number of FASTA files per processing chunk [default %default] |
| `--debug_n` | `0` | Process first N passing OGs only; 0 = all [default %default] |
| `--panel_sizes` | `20,50,100,200,500,1000,2000,5000` | Comma-separated Top-K values for og_panel_topK files [default %default] |
| `--write_tsv_matrix` | `FALSE` | Also write snp_matrix.tsv.gz for each panel [default %default] |
| `--run_p4` | `FALSE` | Also run optional global population-level panel P4 [default %default] |

## Technical notes

Defaults include reference missingness at most 0.20, minor allele count at least 2 and at least three references per class. The highest scoring SNP per ortholog is retained. The estimator uses nonmissing consensus call counts, returns missing for insufficient counts or invalid denominators, and retains finite negative estimates. Scores rank assignment markers; they are not genome wide demographic differentiation estimates.

Required packages: `Biostrings`, `data.table`, `optparse`, `parallel`.

The script records parameters and session information with its outputs. Rerunning with defaults may replace archived files. Preserve a separate archive copy for comparison.

See [conceptual guide](../../docs/CONCEPTS.md), [installation](../../docs/INSTALL.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md).
