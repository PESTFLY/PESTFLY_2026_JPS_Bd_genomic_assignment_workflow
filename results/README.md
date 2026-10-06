# Archived results

Core folders correspond to Steps 00 to 07. `validation/` contains the [seven supplementary analyses](../steps/validation/README.md).

The current authority report is in `07_authority_report/`:

```text
FINAL_origin_tracing_authority_REPORT.xlsx
FINAL_origin_tracing_authority_simplified.tsv
FINAL_origin_tracing_authority_extended.tsv
FINAL_origin_tracing_authority_criteria.tsv
reporting_thresholds.tsv
```

These are native R exports from 5 October 2026. The earlier reports are retained in `legacy_colour_report/` and `reconstruction_checkpoint/`. See [native R verification](../docs/NATIVE_R_VERIFICATION.md).

Step 04 retains complete query trajectories and stability tables. Step 05 contains individual fixed panel leave one out validation; Step 06 contains Random Forest corroboration. Scenario outputs and selected `return_bundle/` exports are retained under the supplementary analyses.

FASTA alignments and the large `01_qc/sample_missingness.tsv` export are excluded; conversion records and the complete sample missingness RDS remain available.

In Step 02 `panels_index.tsv` files, `panel_dir` is relative to its index directory. Other adapted TSV paths are relative to the repository root. Historical binary records, parameter files and some filenames retain original paths and labels. Changes are recorded in `docs/location_index_changes.tsv`, `docs/textual_output_changes.tsv` and `docs/source_file_manifest.tsv`.

Use a separate copy for recomputation because default output paths can replace archived files.
