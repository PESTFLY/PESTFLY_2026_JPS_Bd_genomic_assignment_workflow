# Native R reporting verification

On 5 October 2026, the isolated Hudson FST and MI checks, six criterion reporting checks and Step 07 formatter passed in Windows R 4.5.1, using optparse 1.7.5, data.table 1.18.0 and openxlsx 4.2.8.1. Package build version warnings did not prevent completion.

## Original checker result

The original release checker returned FAIL on 25 descriptive label differences: 21 RF detail cells and four commodity mismatch cells used “Step 4” in the earlier checkpoint and “Step 04” in the native export. The same differences appear in both extended and simplified tables.

`verify_release.R` now normalises that exact spelling only in `rf_corroboration_detail` and `commodity_origin_mismatch_detail`. Other text and scientific comparisons retain their checks. The corrected wrapper has received syntax screening but has not been executed in R; the original FAIL logs remain available.

## Independent export comparison

| Comparison | Result |
| :--- | :--- |
| Criteria, 264 rows and 12 fields | PASS |
| Thresholds, 12 rows and 5 fields | PASS |
| Extended report, 22 rows and 105 compared fields | PASS with exact label normalisation |
| Simplified report, 22 rows and 37 compared fields | PASS with exact label normalisation |
| Four main workbook sheets against native TSV files | PASS across every field |
| Seven retained reference and validation sheets against checkpoint | PASS |
| Both native RDS tables against TSV files | PASS across every field |
| Native run record against Step 04 | PASS |
| All 1,214 archived result MD5 values against the supplied record | PASS |

The review used TSV, OOXML and native R object readers without changing exports or rerunning statistical models. All 22 classes, resolutions and confidence labels were preserved.

## Files and command

`native_R_verification/original_release_check.zip` preserves both supplied run folders. The accompanying logs, accepted label differences and independent comparison records document the review of run `20261005_175802_26116`. Current exports are in `results/07_authority_report/`; earlier reports remain in its two provenance subdirectories.

To run the corrected checker from the repository root:

```bash
Rscript verify_release.R
```

It writes a new `release_check/` folder without replacing archived results. Workbook previews used a disposable copy for unused drawing relationships and empty shared strings; the archived workbook retains the supplied bytes. Excel desktop rendering was not tested.
