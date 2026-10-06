# Metadata snapshot and ortholog QC

Harmonise specimen identifiers and analytical class metadata, verify alignment labels and retain orthologs that satisfy the documented completeness and ambiguity filters.

## Conceptual interpretation

Analytical class labels can represent source lineage affinity rather than physical collection geography. The baseline retains 330 references and 22 queries; the ten Other references are not trained as Africa or Asia. Bdors and Blati are allowed extra alignment labels, not additional metadata specimens.

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

Run from the repository root:

```bash
Rscript steps/01_metadata_snapshot_and_ortholog_qc/run.R
Rscript steps/01_metadata_snapshot_and_ortholog_qc/run.R --help
```

## Parameters

These are script defaults, transcribed from the actual optparse definitions. Archived parameter and run records describe the benchmark execution.

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--align_dir` | `results/00_fasta` | Directory containing Step 00 FASTA alignments [default %default] |
| `--file_glob` | `OG*.fasta` | FASTA file pattern [default %default] |
| `--metadata_xlsx` | `data/000_input_data/metadata.xlsx` | Metadata Excel file [default %default] |
| `--metadata_sheet` | `Selection` | Metadata Excel sheet name [default %default] |
| `--phy_dir` | `data/000_input_data/phy` | Raw PHYLIP directory used for label consistency checks [default %default] |
| `--phy_glob` | `OG*.fa` | Raw PHYLIP file pattern for label consistency checks [default %default] |
| `--out_dir` | `results/01_qc` | Output directory [default %default] |
| `--min_len` | `300` | Minimum alignment length to keep OG [default %default] |
| `--min_nonmissing_frac` | `0.50` | Per-sample minimum non-missing fraction used in occupancy checks [default %default] |
| `--min_pop_occupancy` | `0.70` | Population occupancy threshold: fraction of samples in population represented in OG [default %default] |
| `--min_pop_occupancy_frac_pops` | `0.80` | Required fraction of populations passing min_pop_occupancy [default %default] |
| `--max_mean_ambig` | `0.05` | Maximum mean ambiguous fraction in OG [default %default] |
| `--exclude_regex` | `^$` | Regex of taxa to exclude before QC, if needed [default %default] |
| `--allowed_extra_labels` | `Bdors,Blati` | Comma-separated labels allowed in PHYLIP/FASTA but absent from metadata [default %default] |
| `--cores` | `0` | Parallel workers. 0 = all detected cores minus one [default %default] |
| `--debug_n` | `0` | Process only first N FASTA files; 0 = all [default %default] |
| `--no_sample_missingness` | `FALSE` | Skip writing sample_missingness.tsv/rds [default %default] |
| `--label_check` | `all` | Label consistency check against raw PHYLIP files: all, first, none [default %default] |

## Technical notes

QC defaults are minimum length 300 bp; sample nonmissing fraction 0.50 for population occupancy; population occupancy 0.70 in at least 0.80 of populations; mean ambiguity at most 0.05. Congo and Reunion harmonisation is explicit in the script. Preserve the supplied benchmark classifications when reproducing the paper.

Required packages: `Biostrings`, `data.table`, `optparse`, `parallel`, `readxl`.

The script records parameters and session information with its outputs. Rerunning with defaults may replace archived files. Preserve a separate archive copy for comparison.

See [conceptual guide](../../docs/CONCEPTS.md), [installation](../../docs/INSTALL.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md).
