# Supplementary validation and sensitivity analyses

These analyses supply the material described in Supplementary File S6. Public names describe analytical purpose; original numerical outputs and binary records are retained.

| Analysis | Script directory | Result directory | Main prerequisites |
| :--- | :--- | :--- | :--- |
| 01 Out of scope challenge | `01_out_of_scope_challenge/` | `results/validation/01_out_of_scope_challenge/` | QC metadata, Step 02 panels and the Step 04 script |
| 02 Grouped geographic validation | `02_grouped_geographic_cv/` | `results/validation/02_grouped_geographic_cv/` | QC metadata and Step 02 panels |
| 03 Reference set audit | `03_reference_set_audit/` | `results/validation/03_reference_set_audit/` | QC records and Step 02 panels |
| 04 Mascarene coding sensitivity | `04_mascarene_coding_sensitivity/` | `results/validation/04_mascarene_coding_sensitivity/` | QC, regenerated FASTA and baseline Steps 02, 04 and 05 |
| 05 Balanced reference downsampling | `05_reference_downsampling/` | `results/validation/05_reference_downsampling/` | QC, Step 02 panels and Step 04 assignments |
| 06 Alternative marker set resampling | `06_marker_resampling/` | `results/validation/06_marker_resampling/` | QC, Step 02 panels and Step 04 assignments |
| 07 Reference exclusion sensitivity | `07_reference_exclusion_sensitivity/` | `results/validation/07_exclusion/` | QC, regenerated FASTA and baseline Steps 02, 04 and 05 |

Run each `run.R` from the repository root, for example:

```bash
Rscript steps/validation/02_grouped_geographic_cv/run.R
```

The separate `run_validation_pipeline.R` runner invokes all seven with their defaults. It requires every prerequisite, including FASTA for analyses 04 and 07, and does not run the core workflow first. Use the archived outputs for inspection when alignments are unavailable.

Some analyses contain `return_bundle/` exports of selected summaries. Complete outputs are also retained in the parent directory or `scenarios/`. S6 omits only byte identical export copies with retained counterparts in the same analysis; its duplicate_exports.tsv maps every omission. The full repository retains those copies.

Analysis 07 checks long paths on Windows. Use a short repository location. Its `--out_dir` option also accepts a shorter directory such as `results/v07` when the default path is too long.

Each analysis directory now contains a README with its purpose, scope, complete CLI defaults, inputs and interpretation. See docs/CONCEPTS.md for the distinction between fixed panel validation, training only reranking and complete marker discovery, and docs/S6_SOURCE_MAP.md for the assembled supplementary structure.
