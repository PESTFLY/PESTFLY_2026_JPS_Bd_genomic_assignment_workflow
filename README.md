# PESTFLY genomic assignment workflow

PESTFLY assigns previously identified *Bactrocera dorsalis* specimens to represented reference classes using existing population genomic data. Assignments describe genomic affinity to those classes.

This repository accompanies the Journal of Pest Science manuscript “Genomic assignment of the invasive fruit fly pest *Bactrocera dorsalis*: a repeatable, uncertainty-aware SNP workflow for biosecurity reporting”. The manuscript provides the full methods and biological interpretation; these guides explain how to use the scripts and archived files.

## Release and citation

Release [v1.04](https://github.com/PESTFLY/PESTFLY_2026_JPS_Bd_genomic_assignment_workflow/releases/tag/v1.04), published on 6 October 2026, includes Supplementary Files S5 and S6. Cite the workflow using [10.5281/zenodo.23192148](https://doi.org/10.5281/zenodo.23192148) and [CITATION.cff](CITATION.cff). Release details are in [publication records](docs/GITHUB_PUBLICATION.md).

## Start here

[Installation and versions](docs/INSTALL.md), [conceptual guide](docs/CONCEPTS.md), [metadata](docs/METADATA.md), [reporting criteria](docs/REPORTING_CRITERIA.md) and [S6 source map](docs/S6_SOURCE_MAP.md).

| Location | Contents |
| :--- | :--- |
| `steps/00_*` to `steps/07_*` | Core R scripts and module guides |
| `steps/validation/` | Seven supplementary analyses |
| `data/000_input_data/metadata.xlsx` | Public specimen metadata |
| `results/` | Archived tables, figures, R objects and run records |
| `docs/` | Shared guides, parameter tables and provenance |
| `supplementary/S6/` | S6 reader guide, source map and archive records |

## Core workflow

| Step | Purpose | Result directory |
| :--- | :--- | :--- |
| 00 | PHYLIP to FASTA conversion | `results/00_fasta/` |
| 01 | Metadata and ortholog QC | `results/01_qc/` |
| 02 | SNP discovery and Hudson FST ranking | `results/02_snp_panels/` |
| 03 | Mutual information diagnostics | `results/03_mi_diagnostics/` |
| 04 | Hierarchical genomic assignment | `results/04_origin_assignment/` |
| 05 | Individual fixed panel leave one out validation | `results/05_loo_validation/` |
| 06 | Random Forest corroboration | `results/06_rf_corroboration/` |
| 07 | Reports for plant health authorities | `results/07_authority_report/` |

Run from the repository root after installing the dependencies and restoring the required inputs:

```bash
Rscript run_public_pipeline.R
Rscript run_validation_pipeline.R
```

The runners use module defaults. Use a separate copy for recomputation because default output paths can replace archived files. See [supplementary prerequisites](steps/validation/README.md) before running the second command.

## Data and supplementary files

Public metadata, SNP resources and results can be inspected directly. Raw data are deposited at [10.5281/zenodo.20340447](https://doi.org/10.5281/zenodo.20340447) with restricted file access. Request processed PHYLIP inputs from Massimiliano Virgilio at <massimiliano.virgilio@africamuseum.be>; see [alignment access](ALIGNMENT_ACCESS.md). PHYLIP and generated FASTA files are excluded from the public repository and S6.

S6 contains the seven supplementary analyses and complete query trajectories. See [its reader guide](supplementary/S6/README.md). Build a new archive with `python3 tools/build_s6.py`; the builder records file hashes and maps omitted identical export copies. Published release assets retain their release contents.

The archived reporting checks and formatter passed in R 4.5.1. The original release checker label mismatch and subsequent export review are documented in [native R verification](docs/NATIVE_R_VERIFICATION.md). Historical run records and some filenames retain earlier paths and step labels; [results documentation](results/README.md) explains them.

## Licensing

Copyright (c) 2026 Massimiliano Virgilio and Wannes Dermauw. Code uses [MIT](LICENSE); project documentation, public metadata and derived results use [CC BY 4.0](LICENSE_CONTENT.txt), except where a separate notice applies. See [licensing scope](LICENSING.md). Separately distributed alignments have their own access terms.
