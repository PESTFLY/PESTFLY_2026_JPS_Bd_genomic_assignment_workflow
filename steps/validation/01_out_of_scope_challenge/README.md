# Out of scope geographic challenge

Treat the five Hawaiian and five Papua New Guinean reference specimens labelled Other as unknown queries against the fixed Africa and Asia candidate classes.

## Scope and interpretation

This is a closed set challenge. High Asian affinity for an unrepresented collection population is evidence that support criteria do not reject an unrepresented source, rather than successful validation of true geographic origin.

Uses existing fixed SNP panels and the Step 04 model. The challenge script can reuse existing validated assignment tables; its run record identifies that execution mode. It records the original Excel writer limitation when a conditional branch is empty.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/
steps/04_hierarchical_origin_assignment/run.R
```

## Outputs

Archived output directory: `results/validation/01_out_of_scope_challenge/`.

The analysis retains parameter or design records, summary tables, detailed specimen or replicate results, figures where produced and a run record. Scenario analyses also retain outputs in `scenarios/`. Existing `return_bundle/` directories contain selected exports, not every full output.

## Run

```bash
Rscript steps/validation/01_out_of_scope_challenge/run.R
Rscript steps/validation/01_out_of_scope_challenge/run.R --help
```

Use the repository root as working directory. Scenario analyses 04 and 07 require regenerated FASTA alignments; the other analyses use archived QC and panel resources. A rerun may replace outputs or reuse already validated scenario tables as recorded by the script.

## Parameters

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--metadata` | `results/01_qc/metadata_clean.tsv` | Original cleaned metadata [default %default] |
| `--step3_dir` | `results/02_snp_panels` | Corrected Step 02 panel directory [default %default] |
| `--step4_script` | `steps/04_hierarchical_origin_assignment/run.R` | Existing Step 04 assignment script [default %default] |
| `--out_dir` | `results/validation/01_out_of_scope_challenge` | Separate challenge output directory [default %default] |

Required packages: `data.table`, `optparse`.

These script defaults are distinct from archived execution records. For reproducibility retain those records and fixed seeds. See [conceptual guide](../../../docs/CONCEPTS.md), [installation](../../../docs/INSTALL.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md).
