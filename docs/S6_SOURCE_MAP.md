# Source map for Supplementary File S6

S6 is assembled as `PESTFLY_Supplementary_File_S6.zip`. Its extracted root is `PESTFLY_S6`. The archive contains all seven supplementary analyses, complete baseline query trajectories, and the QC, SNP and script resources needed for inspection and panel based reproduction.

The reader guide explains each design, its findings and its limitations. `source_map.tsv` supplies the exact script and key output paths. In the repository these resources are under `supplementary/S6/`; in the S6 archive they are at the extracted root.

| Section | Analysis | Result directory under `results/validation/` |
| :--- | :--- | :--- |
| S6.1 | Geographic out of scope challenge | `01_out_of_scope_challenge` |
| S6.2 | Grouped geographic cross validation | `02_grouped_geographic_cv` |
| S6.3 | Reference set audit | `03_reference_set_audit` |
| S6.4 | Mascarene coding sensitivity | `04_mascarene_coding_sensitivity` |
| S6.5 | Balanced reference downsampling | `05_reference_downsampling` |
| S6.6 | Alternative marker set resampling | `06_marker_resampling` |
| S6.7 | Reference exclusion sensitivity | `07_exclusion` |

Scenario analyses 04 and 07 retain complete `scenarios/` outputs and unique comparison tables in `return_bundle/`. In analyses 05 and 06, use the main result directory for summary tables and figures. Only byte identical export copies with a retained counterpart inside the same analysis are omitted from S6. `duplicate_exports.tsv` lists every alias. The full repository checkpoint retains those copies.

Query trajectories are the Step 04 P1, P2 and P3 raw tables, their stability tables, final assignments, K grid and parameter records. They contain 656 rows: 352 P1, 192 P2 and 112 P3. Subregion branches exist only for accepted macroregion assignments. Rows at low K remain present, including unavailable posterior values and explicit status fields.

Hudson FST and MI implementation evidence is recorded in the Step 02 and Step 03 run information RDS files. `tools/check_implementation.R` reruns their isolated function checks in base R. The archived self tests passed. The isolated wrapper, six criterion reporting checks and updated Step 07 formatter also passed in native R 4.5.1 on 5 October 2026. Current Step 07 reports are the verified native exports; [native R verification](NATIVE_R_VERIFICATION.md) records the comparison and original checker label mismatch.

Build S6 from a complete repository checkout with Python 3:

```bash
python3 tools/build_s6.py --output PESTFLY_Supplementary_File_S6.zip
```

The builder uses the standard library, selects the resources, removes verified export duplicates, writes an integrity manifest, and checks every ZIP member against its source hash. It does not run R or recompute statistical outputs. Keep the resulting archive outside a public Git history when using the default builder path, or upload it as a release asset.

After extraction, `python3 verify_archive.py` validates all recorded file sizes and hashes. The manifest does not hash itself. Original binary records and parameter logs retain historical paths and internal labels. Presentation changes to six summary text files are documented in `docs/textual_output_changes.tsv`; numerical lines remain unchanged. Alignment inputs are distributed separately upon request. The public code and content licences are included at the S6 archive root.
