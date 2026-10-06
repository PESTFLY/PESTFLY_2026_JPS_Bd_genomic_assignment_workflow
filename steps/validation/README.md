# Supplementary validation and sensitivity analyses

The seven analyses form Supplementary File S6. Each module guide lists inputs, outputs, commands and defaults; the manuscript describes the full methods.

| Analysis | Script directory | Result directory | Main prerequisites |
| :--- | :--- | :--- | :--- |
| 01 Out of scope challenge | `01_out_of_scope_challenge/` | `results/validation/01_out_of_scope_challenge/` | QC metadata, Step 02 panels and the Step 04 script |
| 02 Grouped geographic validation | `02_grouped_geographic_cv/` | `results/validation/02_grouped_geographic_cv/` | QC metadata and Step 02 panels |
| 03 Reference set audit | `03_reference_set_audit/` | `results/validation/03_reference_set_audit/` | QC records and Step 02 panels |
| 04 Mascarene coding sensitivity | `04_mascarene_coding_sensitivity/` | `results/validation/04_mascarene_coding_sensitivity/` | QC, regenerated FASTA and baseline Steps 02, 04 and 05 |
| 05 Balanced reference downsampling | `05_reference_downsampling/` | `results/validation/05_reference_downsampling/` | QC, Step 02 panels and Step 04 assignments |
| 06 Alternative marker set resampling | `06_marker_resampling/` | `results/validation/06_marker_resampling/` | QC, Step 02 panels and Step 04 assignments |
| 07 Reference exclusion sensitivity | `07_reference_exclusion_sensitivity/` | `results/validation/07_exclusion/` | QC, regenerated FASTA and baseline Steps 02, 04 and 05 |

Run one module from the repository root, for example:

```bash
Rscript steps/validation/02_grouped_geographic_cv/run.R
```

`Rscript run_validation_pipeline.R` runs all seven with defaults. It requires their listed resources, including FASTA for 04 and 07, and does not run the core pipeline first. Use a separate copy for recomputation.

Full results are in each analysis directory and its `scenarios/` folders; `return_bundle/` contains selected exports. S6 omits only identical export copies mapped in `duplicate_exports.tsv`; the repository retains them.

Use short Windows paths for nested scenarios. Analysis 07 accepts a shorter `--out_dir`, such as `results/v07`. See [validation scope](../../docs/CONCEPTS.md), [installation](../../docs/INSTALL.md) and [S6 source map](../../docs/S6_SOURCE_MAP.md).
