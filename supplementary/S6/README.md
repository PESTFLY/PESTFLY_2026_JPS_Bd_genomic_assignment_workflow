# Supplementary File S6

Validation and sensitivity analyses for the PESTFLY benchmark. The published v1.04 archive is attached to the [GitHub release](https://github.com/molecular-lab-RMCA/PESTFLY_2026_JPS_Bd_genomic_assignment_workflow/releases/tag/v1.04); its workflow DOI is [10.5281/zenodo.23192148](https://doi.org/10.5281/zenodo.23192148).

S6 contains scripts, parameters, tables, figures and complete outputs for seven supplementary analyses, together with baseline query trajectories. The manuscript provides the full methods and interpretation. `source_map.tsv` lists exact paths; the retained PDF reader guide describes the v1.04 snapshot.

## Analyses


| Section | Analysis | Directory under `results/validation/` |
| :--- | :--- | :--- |
| S6.1 | Geographic out of scope challenge | `01_out_of_scope_challenge` |
| S6.2 | Grouped geographic validation | `02_grouped_geographic_cv` |
| S6.3 | Reference set audit | `03_reference_set_audit` |
| S6.4 | Mascarene coding sensitivity | `04_mascarene_coding_sensitivity` |
| S6.5 | Balanced reference downsampling | `05_reference_downsampling` |
| S6.6 | Alternative marker set resampling | `06_marker_resampling` |
| S6.7 | Reference exclusion sensitivity | `07_exclusion` |


Matching scripts are in `steps/validation/`. Analyses 04 and 07 include complete `scenarios/` outputs and unique comparison exports in `return_bundle/`. For 05 and 06, read summaries and figures in the main result directory. Analysis 07 uses a short directory name for Windows paths.

## Query trajectories and reporting

`results/04_origin_assignment/` contains all evaluated Top K rows:

| Table | Queries | Rows |
| :--- | ---: | ---: |
| `P1_macroregion_raw.tsv` | 22 | 352 |
| `P2_africa_subregion_raw.tsv` | 12 | 192 |
| `P3_asia_subregion_raw.tsv` | 7 | 112 |

Rows retain class support, posterior gap, usable SNP count and status. Low K rows are included. Three unresolved macroregion queries have no evaluated subregion branch. Stability tables, final assignments, the K grid and parameters are in the same directory.

Assignments express affinity to represented classes. The reporting criteria and validation scope are summarised in `docs/REPORTING_CRITERIA.md` and `docs/CONCEPTS.md`; RF remains separate corroboration.

## Inspect or reproduce

Paths below refer to the extracted `PESTFLY_S6` root. TSV files, PDF figures and Excel workbooks can be opened directly; RDS files contain R objects. See `docs/INSTALL.md` and module READMEs for dependencies and commands.

Supplementary analyses 01, 02, 03, 05 and 06 use the included QC and SNP resources. Analyses 04 and 07 require regenerated FASTA. PHYLIP and FASTA inputs are supplied separately upon request from Massimiliano Virgilio at <massimiliano.virgilio@africamuseum.be>. Raw data are deposited at [10.5281/zenodo.20340447](https://doi.org/10.5281/zenodo.20340447) with restricted file access; see `ALIGNMENT_ACCESS.md`. The omitted `sample_missingness.tsv` has a complete RDS counterpart.

Run from a separate extracted copy to preserve archived outputs:

```bash
Rscript run_validation_pipeline.R
```

The runner requires every analysis prerequisite. Implementation and reporting checks need no sequence inputs:

```bash
Rscript tools/check_implementation.R
Rscript steps/07_authority_facing_report/check_reporting.R
```

These checks and the formatter passed in R 4.5.1. `docs/NATIVE_R_VERIFICATION.md` records the original checker mismatch and independent export review.

## Integrity and licensing

After extracting the published archive:

```bash
python3 verify_archive.py
```

`archive_manifest.tsv` records file sizes and SHA256 hashes; `duplicate_exports.tsv` maps omitted identical export copies. The retained manifest, PDF and assembly records describe the published snapshot. Documentation on GitHub's main branch can include later editorial updates; a new build generates its own manifest.

Historical records retain original paths and labels. Code uses MIT; covered documentation, public metadata and results use CC BY 4.0. Copyright (c) 2026 Massimiliano Virgilio and Wannes Dermauw. Separately supplied alignments are outside these public licences. Licence and access notices are included at the archive root.
