# Geographic out of scope challenge

Treat five Hawaiian and five Papua New Guinean Other references as unknown queries against Africa and Asia classes.

All ten benchmark specimens received high confidence Asian affinities. This demonstrates that support criteria cannot detect an unrepresented source. The run record identifies reused assignment tables or an interrupted workbook stage.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/
steps/04_hierarchical_origin_assignment/run.R
```

## Outputs

Archived output directory: `results/validation/01_out_of_scope_challenge/`.

## Run

From the repository root:

```bash
Rscript steps/validation/01_out_of_scope_challenge/run.R
Rscript steps/validation/01_out_of_scope_challenge/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--metadata` | `results/01_qc/metadata_clean.tsv` | Cleaned benchmark metadata |
| `--step3_dir` | `results/02_snp_panels` | Step 02 SNP directory |
| `--step4_script` | `steps/04_hierarchical_origin_assignment/run.R` | Step 04 assignment script |
| `--out_dir` | `results/validation/01_out_of_scope_challenge` | Output directory |

Required packages: `data.table`, `optparse`.

See [installation](../../../docs/INSTALL.md), [concepts](../../../docs/CONCEPTS.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md). The manuscript provides the full methods and interpretation.
