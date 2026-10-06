# PHYLIP to FASTA conversion

Convert existing OMA/read2tree interleaved PHYLIP alignments to FASTA, checking taxon counts and alignment lengths. The input files use the `.fa` suffix but contain PHYLIP data.

Restore the separately supplied PHYLIP inputs before running.

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

From the repository root:

```bash
Rscript steps/00_convert_phylip_to_fasta/run.R
Rscript steps/00_convert_phylip_to_fasta/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--in_dir` | `data/000_input_data/phy` | Input directory with OMA/read2tree PHYLIP MSA files |
| `--out_dir` | `results/00_fasta` | Output directory |
| `--file_glob` | `OG*.fa` | Input file pattern |
| `--max_files` | `0` | For testing: process first N files only; 0 = all |
| `--cores` | `0` | Workers; 0 uses detected cores minus one |
| `--wrap_width` | `80` | FASTA line width |
| `--strict_length` | `FALSE` | If TRUE, stop when reconstructed sequence length != header n_sites |

Required packages: `data.table`, `optparse`, `parallel`.

See [installation](../../docs/INSTALL.md), [concepts](../../docs/CONCEPTS.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md). The manuscript provides the full methods and interpretation.
