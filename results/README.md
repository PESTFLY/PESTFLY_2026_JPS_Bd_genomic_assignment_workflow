# Archived results

The current Step 07 workbook, TSV tables, RDS objects and run record are verified native R exports with explicit six criterion PASS or FAIL reporting. The checks and formatter passed in R 4.5.1 on 5 October 2026. Original colour reports remain under `07_authority_report/legacy_colour_report/`, and the earlier independent reconstruction remains under `07_authority_report/reconstruction_checkpoint/`. See [native R verification](../docs/NATIVE_R_VERIFICATION.md) and [reporting criteria](../docs/REPORTING_CRITERIA.md).

Core folders correspond to Steps 00 to 07. The separate `validation/` folder contains the seven supplementary analyses listed in `steps/validation/README.md`.

Main report files:

```text
results/07_authority_report/FINAL_origin_tracing_authority_REPORT.xlsx
results/07_authority_report/FINAL_origin_tracing_authority_simplified.tsv
results/07_authority_report/FINAL_origin_tracing_authority_extended.tsv
```

Complete query trajectories are retained in Step 04 raw and stability tables. Step 05 supplies individual fixed panel leave one out validation; Step 06 supplies Random Forest corroboration.

Generated FASTA alignments and `results/01_qc/sample_missingness.tsv` are excluded. The conversion summaries, QC RDS counterpart and other QC outputs are retained. The original uploaded revision archive preserves the excluded TSV.

## Portable locations

In each Step 02 `panels_index.tsv`, `panel_dir` is relative to the directory containing that index. Other adapted TSV index locations are relative to the repository root. Changed cells are recorded in `docs/location_index_changes.tsv`.

Original binary objects, workbooks, figures and parameter records are retained unchanged, including the relocated historical Step 07 report files. The current Step 07 exports retain the supplied native R bytes. Stored locations and earlier analysis labels describe the actual runs. They have not been rewritten as if the analyses were generated in this packaged directory.

Six summary text files have presentation edits recorded in docs/textual_output_changes.tsv. Their numerical lines are unchanged. The seven supplementary analyses and all original statistical TSV, RDS, workbook, figure and compressed outputs are unchanged from the preceding checkpoint. Native Step 07 exports replace the current reporting reconstruction, which is retained separately. S6 omits only byte identical summary export copies, mapped in its duplicate_exports.tsv; the full repository retains those copies.

Rerunning with default output directories can replace archived files. Use a separate copy for recomputation.
