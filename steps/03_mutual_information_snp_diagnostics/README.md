# Mutual information diagnostics

Calculate reference based mutual information in bits as a diagnostic comparison with Hudson FST ranking.

MI does not enter marker selection or assignment. Per panel MI tables are written beside the Step 02 panels; the index is written in `results/03_mi_diagnostics/`.

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

From the repository root:

```bash
Rscript steps/03_mutual_information_snp_diagnostics/run.R
Rscript steps/03_mutual_information_snp_diagnostics/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--meta_path` | `results/01_qc/metadata_clean.tsv` | Cleaned metadata file |
| `--step3_dir` | `results/02_snp_panels` | Step 02 SNP directory |
| `--out_dir` | `results/03_mi_diagnostics` | Output directory |
| `--min_class_n` | `3` | Minimum reference specimens per class |
| `--min_snps_keep` | `50` | Minimum SNPs per panel |
| `--include_p4_if_present` | `FALSE` | Include P4 when available |

Required packages: `data.table`, `optparse`.

See [installation](../../docs/INSTALL.md), [concepts](../../docs/CONCEPTS.md) and [reporting criteria](../../docs/REPORTING_CRITERIA.md). The manuscript provides the full methods and interpretation.
