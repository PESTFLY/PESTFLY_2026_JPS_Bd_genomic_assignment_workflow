# Reference set audit

Check identifiers, metadata, class representation, ortholog completeness, consensus profile similarity and published specimen flags.

The audit annotates the primary snapshot without removing references. It identifies twelve published caution references and two with lower completeness for predefined sensitivity scenarios.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/01_qc/
results/02_snp_panels/panels/P1_macroregion_africa_vs_asia/
```

## Outputs

Archived output directory: `results/validation/03_reference_set_audit/`.

## Run

From the repository root:

```bash
Rscript steps/validation/03_reference_set_audit/run.R
Rscript steps/validation/03_reference_set_audit/run.R --help
```

## Parameters

Defaults are listed below; saved parameter and run records describe the benchmark execution.

| Option | Default | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 QC directory |
| `--step2_dir` | `results/02_snp_panels` | Step 02 SNP directory |
| `--out_dir` | `results/validation/03_reference_set_audit` | Output directory |
| `--duplicate_panel` | `P1_macroregion_africa_vs_asia` | Panel for consensus profile comparison |
| `--min_pairwise_snps` | `1000` | Minimum shared SNPs for pairwise comparison |
| `--near_duplicate_concordance` | `0.995` | Near duplicate concordance threshold |
| `--top_pair_rows` | `200` | Most similar pairs retained |
| `--min_fraction_ogs_seen` | `0.95` | Ortholog completeness review threshold |
| `--min_ref_per_group` | `3` | Minimum references per class |

Required packages: `data.table`, `optparse`.

See [installation](../../../docs/INSTALL.md), [concepts](../../../docs/CONCEPTS.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md). The manuscript provides the full methods and interpretation.
