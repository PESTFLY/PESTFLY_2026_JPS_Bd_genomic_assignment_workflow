# Metadata snapshot and ortholog QC

Harmonise specimen identifiers and analytical classes, check alignment labels and retain orthologs meeting QC thresholds.

The benchmark has 330 references and 22 queries. `Bdors` and `Blati` are additional reference genome labels. Retain the supplied class coding when reproducing the benchmark.

## Inputs

```text
data/000_input_data/metadata.xlsx, sheet Selection
data/000_input_data/phy/OG*.fa
results/00_fasta/OG*.fasta
```

## Main outputs

```text
results/01_qc/metadata_clean.tsv
results/01_qc/og_qc.tsv
results/01_qc/ogs_pass_qc.txt
results/01_qc/step1_2_run_info.rds
```

## Run

From the repository root:

```bash
Rscript steps/01_metadata_snapshot_and_ortholog_qc/run.R
Rscript steps/01_metadata_snapshot_and_ortholog_qc/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--align_dir` | `results/00_fasta` | Step 00 FASTA directory |
| `--file_glob` | `OG*.fasta` | Input file pattern |
| `--metadata_xlsx` | `data/000_input_data/metadata.xlsx` | Metadata Excel file |
| `--metadata_sheet` | `Selection` | Metadata Excel sheet name |
| `--phy_dir` | `data/000_input_data/phy` | Raw PHYLIP directory used for label consistency checks |
| `--phy_glob` | `OG*.fa` | Raw PHYLIP file pattern for label consistency checks |
| `--out_dir` | `results/01_qc` | Output directory |
| `--min_len` | `300` | Minimum ortholog length |
| `--min_nonmissing_frac` | `0.50` | Per-sample minimum non-missing fraction used in occupancy checks |
| `--min_pop_occupancy` | `0.70` | Minimum represented fraction per population |
| `--min_pop_occupancy_frac_pops` | `0.80` | Fraction of populations meeting occupancy threshold |
| `--max_mean_ambig` | `0.05` | Maximum mean ambiguous fraction |
| `--exclude_regex` | `^$` | Taxon exclusion pattern |
| `--allowed_extra_labels` | `Bdors,Blati` | Alignment labels allowed outside metadata |
| `--cores` | `0` | Workers; 0 uses detected cores minus one |
| `--debug_n` | `0` | Process first N orthologs; 0 means all |
| `--no_sample_missingness` | `FALSE` | Skip writing sample_missingness.tsv/rds |
| `--label_check` | `all` | Label consistency check against raw PHYLIP files: all, first, none |

Required packages: `Biostrings`, `data.table`, `optparse`, `parallel`, `readxl`.

See [installation](../../docs/INSTALL.md), [concepts](../../docs/CONCEPTS.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md). The manuscript provides the full methods and interpretation.
