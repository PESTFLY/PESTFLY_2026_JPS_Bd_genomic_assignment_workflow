# Reference exclusion sensitivity

Compare the retained primary benchmark with exclusions of twelve published caution references, two lower completeness references, and their combined set of fourteen.

## Scope and interpretation

These are predefined sensitivity scenarios, not post hoc removals chosen to improve query assignments. Changes quantify dependence on the retained references; unchanged calls do not prove that every retained specimen is taxonomically or biologically unproblematic.

Scenario runs repeat Steps 02, 04 and 05 and require FASTA alignments. Results use the short 07_exclusion directory to limit Windows path length. Preexisting validated scenario outputs may be reused and are identified by run records.

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

The analysis retains parameter or design records, summary tables, detailed specimen or replicate results, figures where produced and a run record. Scenario analyses also retain outputs in `scenarios/`. Existing `return_bundle/` directories contain selected exports, not every full output.

## Run

```bash
Rscript steps/validation/07_reference_exclusion_sensitivity/run.R
Rscript steps/validation/07_reference_exclusion_sensitivity/run.R --help
```

Use the repository root as working directory. Scenario analyses 04 and 07 require regenerated FASTA alignments; the other analyses use archived QC and panel resources. A rerun may replace outputs or reuse already validated scenario tables as recorded by the script.

## Parameters

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Primary Step 01 QC directory [default %default] |
| `--align_dir` | `results/00_fasta` | Step 00 FASTA directory [default %default] |
| `--baseline_step2_dir` | `results/02_snp_panels` | Corrected baseline Step 02 directory [default %default] |
| `--baseline_step4_dir` | `results/04_origin_assignment` | Corrected baseline Step 04 directory [default %default] |
| `--baseline_step5_dir` | `results/05_loo_validation` | Corrected baseline Step 05 directory [default %default] |
| `--out_dir` | `results/validation/07_exclusion` | Supplementary validation output directory [default %default] |
| `--cores` | `0L` | Step 02 workers; 0 uses all detected cores minus one [default %default] |
| `--chunk_size` | `300L` | Step 02 FASTA files per processing chunk [default %default] |
| `--resume` | `TRUE` | Reuse complete scenario steps after an interrupted run [default %default] |

Required packages: `data.table`, `optparse`.

These script defaults are distinct from archived execution records. For reproducibility retain those records and fixed seeds. See [conceptual guide](../../../docs/CONCEPTS.md), [installation](../../../docs/INSTALL.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md).
