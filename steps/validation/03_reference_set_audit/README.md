# Reference set audit

Audit identifiers, metadata, reference composition, ortholog completeness, SNP profile concordance, source collections and published caution flags.

## Scope and interpretation

Caution flags are annotations from the source literature, not proof that a specimen is invalid. The primary benchmark keeps these specimens; sensitivity analysis tests predefined exclusions. Geography and source collection can be confounded.

Near duplicate screening uses P1 consensus calls, at least 1000 shared calls and concordance at least 0.995. Two references have lower ortholog completeness and twelve carry published caution flags. The audit itself does not change the primary dataset.

## Inputs

```text
results/01_qc/metadata_clean.tsv
results/01_qc/
results/02_snp_panels/panels/P1_macroregion_africa_vs_asia/
```

## Outputs

Archived output directory: `results/validation/03_reference_set_audit/`.

The analysis retains parameter or design records, summary tables, detailed specimen or replicate results, figures where produced and a run record. Scenario analyses also retain outputs in `scenarios/`. Existing `return_bundle/` directories contain selected exports, not every full output.

## Run

```bash
Rscript steps/validation/03_reference_set_audit/run.R
Rscript steps/validation/03_reference_set_audit/run.R --help
```

Use the repository root as working directory. Scenario analyses 04 and 07 require regenerated FASTA alignments; the other analyses use archived QC and panel resources. A rerun may replace outputs or reuse already validated scenario tables as recorded by the script.

## Parameters

| Option | Default R expression | Meaning |
| :--- | :--- | :--- |
| `--qc_dir` | `results/01_qc` | Step 01 QC directory [default %default] |
| `--step2_dir` | `results/02_snp_panels` | Corrected Step 02 panel directory [default %default] |
| `--out_dir` | `results/validation/03_reference_set_audit` | Audit output directory [default %default] |
| `--duplicate_panel` | `P1_macroregion_africa_vs_asia` | Panel matrix used for genotype duplicate checks [default %default] |
| `--min_pairwise_snps` | `1000L` | Minimum pairwise nonmissing SNPs for similarity checks [default %default] |
| `--near_duplicate_concordance` | `0.995` | Concordance threshold for near duplicate flags [default %default] |
| `--top_pair_rows` | `200L` | Most similar reference pairs retained for audit [default %default] |
| `--min_fraction_ogs_seen` | `0.95` | Ortholog completeness review threshold [default %default] |
| `--min_ref_per_group` | `3L` | Minimum training references per class for holdout feasibility [default %default] |

Required packages: `data.table`, `optparse`.

These script defaults are distinct from archived execution records. For reproducibility retain those records and fixed seeds. See [conceptual guide](../../../docs/CONCEPTS.md), [installation](../../../docs/INSTALL.md) and [S6 source map](../../../docs/S6_SOURCE_MAP.md).
