# Supplementary File S6 source map

`PESTFLY_Supplementary_File_S6.zip` extracts to `PESTFLY_S6`. It includes seven supplementary analyses, baseline query trajectories, QC records, SNP resources and scripts. The manuscript provides the full analytical methods.


| Section | Analysis | Directory under `results/validation/` |
| :--- | :--- | :--- |
| S6.1 | Geographic out of scope challenge | `01_out_of_scope_challenge` |
| S6.2 | Grouped geographic validation | `02_grouped_geographic_cv` |
| S6.3 | Reference set audit | `03_reference_set_audit` |
| S6.4 | Mascarene coding sensitivity | `04_mascarene_coding_sensitivity` |
| S6.5 | Balanced reference downsampling | `05_reference_downsampling` |
| S6.6 | Alternative marker set resampling | `06_marker_resampling` |
| S6.7 | Reference exclusion sensitivity | `07_exclusion` |


[source_map.tsv](../supplementary/S6/source_map.tsv) lists exact script and key output paths. Analyses 04 and 07 retain scenario outputs in `scenarios/` and comparison summaries in `return_bundle/`. For 05 and 06, use the main result directory for summaries and figures.

Step 04 trajectories contain 656 rows: 352 P1, 192 P2 and 112 P3. The same directory holds stability tables, final assignments, the K grid and parameters. Low K rows remain included; subregion branches depend on a reportable macroregion call.

Build a new archive from a complete checkout:

```bash
python3 tools/build_s6.py --output PESTFLY_Supplementary_File_S6.zip
```

The builder records size and SHA256 for each archive file except the manifest itself. It maps omitted identical exports in `duplicate_exports.tsv` and checks ZIP contents. It does not run R. Keep generated archives outside Git history or attach them to a new release.

After extraction, `python3 verify_archive.py` checks the archive manifest. The manifest and PDF retained in this repository describe the published v1.04 S6 snapshot; rebuilding generates a manifest for the current files. Historical run records preserve original paths and labels.

See the [S6 guide](../supplementary/S6/README.md), [installation](INSTALL.md), [native R checks](NATIVE_R_VERIFICATION.md) and [alignment access](../ALIGNMENT_ACCESS.md). Sequence alignments are supplied separately upon request. Public licences are included at the archive root.
