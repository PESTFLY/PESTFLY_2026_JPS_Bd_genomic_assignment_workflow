# Mutual information diagnostics

Calculate reference only mutual information between the consensus allele and target class, in bits, as a diagnostic comparison with Hudson ranking.

## Conceptual interpretation

Mutual information does not select, reorder or replace the Step 02 marker panels and does not enter the Step 04 likelihood. It is a descriptive comparison between two informativeness statistics.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/02_snp_panels/panels_index.tsv
results/02_snp_panels/panels/<panel_id>/snp_matrix.rds
```

## Main outputs

```text
results/03_mi_diagnostics/panels_mi_index.tsv
results/02_snp_panels/panels/<panel_id>/snp_mi.tsv
results/03_mi_diagnostics/step03_run_info.rds
```

## Run

Run from the repository root:

```bash
Rscript steps/03_mutual_information_snp_diagnostics/run.R
Rscript steps/03_mutual_information_snp_diagnostics/run.R --help
```

## Parameters

These are script defaults, transcribed from the actual optparse definitions. Archived parameter and run records describe the benchmark execution.

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--meta_path` | `results/01_qc/metadata_clean.tsv` | Path to metadata_clean.tsv [default %default] |
| `--step3_dir` | `results/02_snp_panels` | paste( "Step 02 output directory containing panels_index.tsv", "[default %default]" ) |
| `--out_dir` | `results/03_mi_diagnostics` | paste( "Step 03 output directory for MI index and run information", "[default %default]" ) |
| `--min_class_n` | `3` | paste( "Drop target classes with fewer than this many reference samples", "[default %default]" ) |
| `--min_snps_keep` | `50` | Skip panels with fewer than this many SNPs [default %default] |
| `--include_p4_if_present` | `FALSE` | Also compute MI for optional P4 if present [default %default] |

## Technical notes

The script removes unused class levels and includes MI self checks. MI output tables are written next to their Step 02 panels; the global index is written under results/03_mi_diagnostics.

Required packages: `data.table`, `optparse`.

The script records parameters and session information with its outputs. Rerunning with defaults may replace archived files. Preserve a separate archive copy for comparison.

See [conceptual guide](../../docs/CONCEPTS.md), [installation](../../docs/INSTALL.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md).
