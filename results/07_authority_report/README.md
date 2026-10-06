# Six criterion authority report

The current workbook, TSV tables, RDS objects and `step7_run_info.rds` are the native R exports from 5 October 2026. They report 19 macroregion PASS and three FAIL outcomes, plus 12 subregion PASS, seven FAIL and three NOT_EVALUATED outcomes. All 22 reported classes, resolutions and confidence categories are unchanged.

The updated formatter and its base R checks passed in the supplied R 4.5.1 session. Independent review verified all criterion and threshold values, every native main workbook sheet and both RDS tables. The original release checker stopped on descriptive step label spelling; its logs and the completed comparison are recorded in [native R verification](../../docs/NATIVE_R_VERIFICATION.md). See the [Step 07 guide](../../steps/07_authority_facing_report/README.md). No genomic or statistical model was recomputed.

`legacy_colour_report/` preserves all six original colour report files, including the workbook, TSV tables, RDS objects and historical R run record. `reconstruction_checkpoint/` preserves the independently assembled Step 4 and Step 5 report and its reconstruction information. Use the parent directory files for the current native criterion reports.
