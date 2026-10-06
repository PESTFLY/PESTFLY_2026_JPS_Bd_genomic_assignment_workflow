# PESTFLY genomic assignment workflow

This repository supports the manuscript “Genomic assignment of the invasive fruit fly pest *Bactrocera dorsalis*: a repeatable, uncertainty-aware SNP workflow for biosecurity reporting”, submitted to Journal of Pest Science.

PESTFLY compares previously identified *Bactrocera dorsalis* specimens with represented reference classes using existing population genomic data. Assignments describe genomic affinity to those classes. The workflow does not establish actual collection origin or a transport pathway, or determine whether the true source is represented.

## Package status

Release [v1.03](https://github.com/PESTFLY/PESTFLY_origin_tracing_publication_pipeline/releases/tag/v1.03) was published on 6 October 2026 from the Step 8 checkpoint prepared on 5 October 2026. All eight core modules and seven supplementary analyses have reader guides, conceptual notes and CLI tables. Supplementary File S6 contains complete outputs, query trajectories and integrity checks. The isolated Hudson FST and MI checks, six criterion reporting checks and Step 07 formatter passed in R 4.5.1 on Windows. The current Step 07 workbook, TSV tables, RDS objects and run record are the verified native R exports. The recorded outcomes remain 19 macroregion passes and 12 subregion passes, with all 22 original reported classes and confidence labels preserved. Original colour reports and the independently assembled checkpoint remain available as provenance. The code uses MIT and the project documentation and public results use CC BY 4.0, with copyright notices naming Massimiliano Virgilio and Wannes Dermauw. Input alignments remain available upon request. The workflow DOI for v1.03 is confirmed; the alignment deposit metadata remains to be confirmed. See [native R verification](docs/NATIVE_R_VERIFICATION.md) for the original checker label mismatch and the completed export review. Statistical models were not recomputed.

## Structure

| Location | Contents |
| :--- | :--- |
| `run_public_pipeline.R` | Runner for core Steps 00 to 07 |
| `steps/00_*` to `steps/07_*` | Eight core scripts with existing README files |
| `steps/validation/` | Seven supplementary validation and sensitivity scripts |
| `run_validation_pipeline.R` | Separate runner for the seven supplementary analyses |
| `data/000_input_data/metadata.xlsx` | Public sample metadata |
| `data/000_input_data/phy/` | Expected location for separately supplied PHYLIP alignments |
| `results/00_fasta/` to `results/07_authority_report/` | Archived core outputs |
| `results/validation/` | Archived supplementary outputs |
| `docs/` | Packaging notes, source hashes and location change records |
| `supplementary/S6/` | Reader guide, source map and archive records |
| `tools/build_s6.py` | Rebuild the S6 archive using Python 3 |
| `verify_release.R` | Check implementation and reporting without rerunning models |

## Reader guides

[Concepts](docs/CONCEPTS.md), [installation and recorded versions](docs/INSTALL.md), [metadata](docs/METADATA.md), [six reporting criteria](docs/REPORTING_CRITERIA.md), and [S6 source map](docs/S6_SOURCE_MAP.md).

## Core workflow

| Step | Purpose | Result directory |
| :--- | :--- | :--- |
| 00 | PHYLIP to FASTA conversion | `results/00_fasta/` |
| 01 | Metadata and ortholog QC | `results/01_qc/` |
| 02 | SNP panel discovery and ranking | `results/02_snp_panels/` |
| 03 | Mutual information diagnostics | `results/03_mi_diagnostics/` |
| 04 | Hierarchical genomic assignment | `results/04_origin_assignment/` |
| 05 | Individual fixed panel leave one out validation | `results/05_loo_validation/` |
| 06 | Random Forest corroboration | `results/06_rf_corroboration/` |
| 07 | Reports for plant health authorities | `results/07_authority_report/` |

Run from the repository root after installing the R dependencies and restoring the required alignment inputs:

```bash
Rscript run_public_pipeline.R
```

Run supplementary analyses separately:

```bash
Rscript run_validation_pipeline.R
```

See [the supplementary analysis guide](steps/validation/README.md) for prerequisites. Both runners use the scripts' defaults. Rerunning can replace archived outputs; use a separate copy for recomputation.

## Data access

Public metadata, derived SNP panels and archived reports can be inspected with the supplied files. Reproducing Steps 00 to 02 and supplementary analyses 04 and 07 requires the separately distributed alignments or regenerated FASTA files.

The workflow archive for release v1.03 is identified by [10.5281/zenodo.23187057](https://doi.org/10.5281/zenodo.23187057). The revised manuscript identifies the separately distributed alignment deposit as <https://doi.org/10.5281/zenodo.20283931>, with files available upon request. Its publication metadata and access terms remain to be confirmed. The public code repository is <https://github.com/PESTFLY/PESTFLY_origin_tracing_publication_pipeline>. The earlier manuscript workflow identifier and its unconfirmed relationship to the current record are retained in [publication records](docs/GITHUB_PUBLICATION.md).

## Archived names and provenance

The three Step 00 conversion records from the original results upload are retained unchanged and included in the source file manifest. They record original machine locations as run provenance.

Directories use the public numbering. Some output basenames and compatible CLI options retain earlier identifiers, including `step3`, `step4b` and `step5`, because downstream scripts use them. See `docs/packaging_changes.txt`.

Adapted TSV indexes have portable locations. Original RDS objects and historical parameter records are preserved and can contain original machine locations and analysis labels. Location changes are recorded in `docs/location_index_changes.tsv`. Source and packaged hashes are recorded in `docs/source_file_manifest.tsv`.

## Supplementary File S6

S6 contains the seven supplementary analyses with scripts, parameters, tables, figures and complete outputs, plus query trajectories across all evaluated SNP subsets. See [the archive reader guide](supplementary/S6/README.md), [the PDF guide](supplementary/S6/S6_reader_guide.pdf), and [the source map](docs/S6_SOURCE_MAP.md).

Build the archive with `python3 tools/build_s6.py`. The archive retains the repository layout and baseline resources. It omits only verified byte identical summary export copies, with each omission mapped to the retained file. The full repository retains those export copies. No statistical analysis was recomputed for this assembly.

## GitHub publication

The current release is [v1.03](https://github.com/PESTFLY/PESTFLY_origin_tracing_publication_pipeline/releases/tag/v1.03), published on 6 October 2026 from commit `f397f8e0a269fa089140eff22a00618d3583dd3c`. Both attached S5 and S6 files match the verified publication assets. The package preserves file bytes through Git checkin and checkout so its recorded hashes remain valid on Windows. See [publication records](docs/GITHUB_PUBLICATION.md). Citation and publication documentation on `main` record the completed release; the release tag identifies its original snapshot.

## Licensing

Copyright (c) 2026 Massimiliano Virgilio and Wannes Dermauw. Code is licensed under [MIT](LICENSE). Project documentation and public metadata and derived results use [CC BY 4.0](LICENSE_CONTENT.txt), except where a separate notice applies. See [licensing scope](LICENSING.md) and [alignment access upon request](ALIGNMENT_ACCESS.md). Full reproduction of the alignment stages requires obtaining the separately distributed inputs.
