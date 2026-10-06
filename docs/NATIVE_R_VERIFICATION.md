# Native R reporting verification

On 5 October 2026, the isolated Hudson FST and MI implementation checks, six criterion reporting checks and updated Step 07 formatter passed in a supplied Windows R 4.5.1 session. The native run used optparse 1.7.5, data.table 1.18.0 and openxlsx 4.2.8.1, matching the recorded main reporting package versions. The session printed warnings that those packages had been built under R 4.5.2; all three tasks completed successfully.

## Checker correction

The original checker returned FAIL after its extended table comparison found 25 descriptive cells spelling the step label differently. Twenty one RF detail cells and four commodity mismatch detail cells used “Step 4” in the independently assembled checkpoint and “Step 04” in the native export. Those same 25 label differences also appear in the simplified export.

The correction in `verify_release.R` normalises exactly `Step 4 ` to `Step 04 ` only in `rf_corroboration_detail` and `commodity_origin_mismatch_detail`. Other wording in these fields is still checked. Scientific values, categories, missing values and thresholds retain their previous comparison rules. The original four reporting prose exclusions apply only to checkpoint comparisons. All fields, including prose, are compared between the native workbook, TSV and RDS outputs.

## Completed export review

An independent review of the supplied files completed the comparisons that the original wrapper did not reach. The review used Python standard library TSV and OOXML readers and an XDR reader for the native R vectors and records. It made no changes to the supplied exports.

| Comparison | Result |
| :--- | :--- |
| Criterion table against checkpoint, 264 rows and 12 fields | PASS |
| Threshold table against checkpoint, 12 rows and 5 fields | PASS |
| Extended report against checkpoint, 22 rows and 105 compared fields | PASS with the exact step label normalisation |
| Simplified report against checkpoint, 22 rows and 37 compared fields | PASS with the exact step label normalisation |
| All four main workbook sheets against native TSV exports | PASS across every field |
| Seven retained reference and validation workbook sheets against checkpoint | PASS |
| Simplified and extended native RDS tables against native TSV exports | PASS across every field |
| Native run record, 22 samples and consistency with Step 04 | PASS |
| All 1,214 archived result file MD5 values against supplied pre run record | PASS |

Macroregion outcomes remain 19 PASS and three FAIL. Subregion outcomes remain 12 PASS, seven FAIL and three NOT_EVALUATED. All 22 reported classes, resolutions and confidence labels are unchanged. No genomic or statistical model was rerun.

## Records and current files

`docs/native_R_verification/original_release_check.zip` preserves both uploaded run directories without changing their FAIL logs. The completed review uses `20261005_175802_26116`, corresponding to the supplied console output. Its original log, session, package versions and comparison files are also available beside the independent comparison summary, accepted label differences and JSON verification report. The original overall FAIL is preserved and is separate from the subsequent export review PASS.

Current native exports are in `results/07_authority_report/`. The previous independent reconstruction is in `reconstruction_checkpoint/`; original colour reports remain in `legacy_colour_report/`. Native run records retain the machine paths of the actual execution as provenance.

Workbook previews used a disposable copy to accommodate 26 unused drawing relationships and empty shared strings in the openxlsx export that the strict preview reader could not import directly. The worksheet values and styles were retained for preview. The archived native workbook remains byte identical to the supplied export. Excel desktop rendering was not tested.

The corrected wrapper has been syntax screened in the assembly workspace, where R is unavailable. It has not been executed there. The supplied native component checks and formatter execution plus the completed independent file comparisons establish this checkpoint without another user run. Future reporting checks can use `Rscript verify_release.R` from the extracted repository root.

The selected code and content licences and alignment access policy are recorded in [LICENSING.md](../LICENSING.md). External release record checks remain outstanding.
