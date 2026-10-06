# Reference exclusion sensitivity

Compare the primary snapshot with exclusions of twelve published caution references, two lower completeness references and their combined set of fourteen.

Each predefined scenario repeats Steps 02, 04 and 05 and requires FASTA. The run record identifies reused outputs. The shorter `07_exclusion` result directory accommodates Windows paths.

## Inputs

```text
results/01_qc/
results/00_fasta/OG*.fasta
results/02_snp_panels/
results/04_origin_assignment/
results/05_loo_validation/
```

## Outputs

Archived output directory: `results/validation/07_exclusion/`.

## Run

From the repository root:

```bash
Rscript steps/validation/07_reference_exclusion_sensitivity/run.R
Rscript steps/validation/07_reference_exclusion_sensitivity/run.R --help
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
| `--out_dir` | `results/validation/07_exclusion` | Output directory |
| `--cores` | `0` | Workers; 0 uses detected cores minus one |
| `--chunk_size` | `300` | FASTA files per processing chunk |
| `--resume` | `TRUE` | Reuse complete scenario steps |

Required packages: `data.table`, `optparse`.

See [installation](../../../docs/INSTALL.md), [concepts](../../../docs/CONCEPTS.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md). The manuscript provides the full methods and interpretation.
