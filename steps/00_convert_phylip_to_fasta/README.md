# PHYLIP to FASTA conversion

Convert the supplied OMA/read2tree interleaved PHYLIP alignments to FASTA while checking taxon counts and alignment lengths.

## Conceptual interpretation

This is format conversion of existing alignments. It does not process raw reads, estimate orthology or rediscover loci. The .fa suffix of the input files does not identify their format: these inputs are PHYLIP.

## Inputs

```text
data/000_input_data/phy/OG*.fa
```

## Main outputs

```text
results/00_fasta/OG*.fasta
results/00_fasta/convert_summary.tsv
results/00_fasta/step0_run_info.rds
```

## Run

Run from the repository root:

```bash
Rscript steps/00_convert_phylip_to_fasta/run.R
Rscript steps/00_convert_phylip_to_fasta/run.R --help
```

## Parameters

These are script defaults, transcribed from the actual optparse definitions. Archived parameter and run records describe the benchmark execution.

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--in_dir` | `data/000_input_data/phy` | Input directory with OMA/read2tree PHYLIP MSA files [default %default] |
| `--out_dir` | `results/00_fasta` | Output directory for converted FASTA files [default %default] |
| `--file_glob` | `OG*.fa` | Input PHYLIP file pattern [default %default] |
| `--max_files` | `0` | For testing: process first N files only; 0 = all [default %default] |
| `--cores` | `0` | Parallel workers. 0 = all detected cores minus one [default %default] |
| `--wrap_width` | `80` | FASTA line width [default %default] |
| `--strict_length` | `FALSE` | If TRUE, stop when reconstructed sequence length != header n_sites [default %default] |

## Technical notes

Restoring the separately deposited PHYLIP inputs is required. The generated FASTA alignments are omitted from the public package.

Required packages: `data.table`, `optparse`, `parallel`.

The script records parameters and session information with its outputs. Rerunning with defaults may replace archived files. Preserve a separate archive copy for comparison.

See [conceptual guide](../../docs/CONCEPTS.md), [installation](../../docs/INSTALL.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md).
