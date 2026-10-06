# Mascarene coding sensitivity

Compare baseline Asia and Southeast Asia coding of 40 Reunion and five Mauritius references with Africa and East Africa coding and with Other coding.

The baseline follows the lineage interpretation in the manuscript. Each scenario repeats Steps 02, 04 and 05 in its own directory. FASTA alignments are required; completed scenario outputs can be reused as recorded by the script.

## Inputs

```text
results/01_qc/
results/00_fasta/OG*.fasta
results/02_snp_panels/
results/04_origin_assignment/
results/05_loo_validation/
```

## Outputs

Archived output directory: `results/validation/04_mascarene_coding_sensitivity/`.

## Run

From the repository root:

```bash
Rscript steps/validation/04_mascarene_coding_sensitivity/run.R
Rscript steps/validation/04_mascarene_coding_sensitivity/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 QC directory |
| `--align_dir` | `results/00_fasta` | Step 00 FASTA directory |
| `--baseline_step2_dir` | `results/02_snp_panels` | Baseline Step 02 directory |
| `--baseline_step4_dir` | `results/04_origin_assignment` | Baseline Step 04 directory |
| `--baseline_step5_dir` | `results/05_loo_validation` | Baseline Step 05 directory |
| `--out_dir` | `results/validation/04_mascarene_coding_sensitivity` | Output directory |
| `--cores` | `0` | Workers; 0 uses detected cores minus one |
| `--chunk_size` | `300` | FASTA files per processing chunk |
| `--resume` | `TRUE` | Reuse complete scenario steps |

Required packages: `data.table`, `optparse`.

See [installation](../../../docs/INSTALL.md), [concepts](../../../docs/CONCEPTS.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md). The manuscript provides the full methods and interpretation.
