# SNP discovery and Hudson FST ranking

Discover biallelic reference SNPs, retain the highest scoring SNP per ortholog and rank each panel by Hudson FST.

Consensus alleles are encoded as 0 and 1; ambiguous states are missing. P1 uses Africa versus Asia; P2 and P3 use the maximum pairwise subregion score. Queries do not select markers. Finite negative estimates are retained.

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

From the repository root:

```bash
Rscript steps/02_snp_panel_discovery_and_ranking/run.R
Rscript steps/02_snp_panel_discovery_and_ranking/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--align_dir` | `results/00_fasta` | Step 00 FASTA directory |
| `--file_glob` | `OG*.fasta` | Input file pattern |
| `--qc_dir` | `results/01_qc` | Step 01 QC directory |
| `--out_dir` | `results/02_snp_panels` | Output directory |
| `--reference_value` | `auto` | Which is_reference value means reference DB: auto, TRUE, or FALSE |
| `--macroregions` | `Africa,Asia` | Comma-separated macroregion_3 values used for P1/P2/P3 |
| `--max_site_missing` | `0.20` | Maximum reference missing fraction per SNP |
| `--min_mac` | `2` | Minimum reference minor allele count |
| `--min_ref_per_group` | `3` | Minimum references per class |
| `--max_snps_per_og` | `1` | Maximum SNPs per ortholog and panel |
| `--exclude_regex` | `^$` | Taxon exclusion pattern |
| `--allowed_extra_labels` | `Bdors,Blati` | Alignment labels allowed outside metadata |
| `--cores` | `0` | Workers; 0 uses detected cores minus one |
| `--chunk_size` | `300` | FASTA files per processing chunk |
| `--debug_n` | `0` | Process first N orthologs; 0 means all |
| `--panel_sizes` | `20,50,100,200,500,1000,2000,5000` | Marker panel sizes |
| `--write_tsv_matrix` | `FALSE` | Also write compressed TSV SNP matrices |
| `--run_p4` | `FALSE` | Generate optional population panel P4 |

Required packages: `Biostrings`, `data.table`, `optparse`, `parallel`.

See [installation](../../docs/INSTALL.md), [concepts](../../docs/CONCEPTS.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md). The manuscript provides the full methods and interpretation.
