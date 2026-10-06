# Runtime and installation

## Recorded benchmark environment

These versions match the retained R session records and the revised manuscript. The detailed record locations are in `runtime_version_evidence.tsv`.

| Component | Recorded version |
| :--- | :--- |
| R | 4.5.1 |
| data.table | 1.18.0 |
| optparse | 1.7.5 |
| readxl | 1.4.5 |
| Biostrings | 2.76.0 |
| openxlsx | 4.2.8.1 |
| ranger | 0.17.0 |
| randomForest | 4.7-1.2 |

`parallel`, `stats`, `utils`, `tools` and `grDevices` are supplied with R. Archived runs used R 4.5.1 on Windows. The RF record identifies ranger as the classifier backend; the script otherwise falls back to randomForest when ranger is unavailable. Reproducing its fitted results requires retaining the backend and recorded seeds.

## Dependency setup

Use an R installation appropriate for the recorded package versions. The following R commands illustrate dependency setup using configured CRAN repositories and Bioconductor:

```r
install.packages(c("data.table", "optparse", "readxl", "openxlsx",
                   "ranger", "randomForest", "BiocManager"))
BiocManager::install("Biostrings", version = "3.21",
                     ask = FALSE, update = FALSE)
```

Bioconductor 3.21 is compatible with R 4.5. The example installs versions available from the configured repositories; it is not a complete restoration of the original dependency snapshot. Compare installed versions with `environment/benchmark_versions.tsv` before claiming exact benchmark reproduction. A complete executable environment or transitive dependency lock has not been reconstructed from these session records.

Check installed components from the repository root:

```bash
Rscript environment/check_environment.R
Rscript environment/check_environment.R --strict
```

The default check reports version differences and fails when required dependencies are missing. Strict mode also fails for any mismatch with the recorded main versions. Installation commands and the new checker have not been executed in this packaging environment; they have received syntax screening.

## Reproduction routes

| Route | Required resources | What can be done |
| :--- | :--- | :--- |
| Inspect saved outputs | Supplied metadata and result archive | Read tables, figures, parameters and reports |
| Rerun from archived SNP panels | Required R packages and archived QC and SNP resources | Run Steps 03 to 07 and supplementary analyses 01, 02, 03, 05 and 06 |
| Complete core reproduction | PHYLIP inputs, metadata and R environment | Regenerate Step 00 FASTA and Steps 01 to 07 |
| Scenario panel reconstruction | Regenerated FASTA, QC and baseline outputs | Run Mascarene coding and reference exclusion scenarios |

The associated raw data deposit is [restricted Zenodo record 10.5281/zenodo.20340447](https://doi.org/10.5281/zenodo.20340447). Its metadata are public; its files require authorised access. For complete core reproduction, obtain the processed PHYLIP inputs from Massimiliano Virgilio and copy them into `data/000_input_data/phy/`. See [alignment access](../ALIGNMENT_ACCESS.md) for requests and the distinction between the raw data deposit and processed alignment inputs. A full run uses `Rscript run_public_pipeline.R`; supplementary analyses use `Rscript run_validation_pipeline.R`.

Both runners use defaults. A rerun can overwrite archived outputs. Keep a separate copy for inspection and comparison. Step 07 passed native R execution on 5 October 2026 using R 4.5.1, data.table 1.18.0, optparse 1.7.5 and openxlsx 4.2.8.1. Its exported tables, workbook and RDS objects were checked against the assembled checkpoint and each other. See [native R verification](NATIVE_R_VERIFICATION.md). Its formatter runs directly from saved Step 04 assignments and parameters, without alignment inputs.

## Computing and seeds

Core alignment stages accept `--cores`; zero generally requests detected cores minus one. Supplementary grouped validation, downsampling and marker resampling accept `--n_cores`; their automatic worker cap is eight. Actual parallel strategies differ by module and operating system, so specify workers explicitly when recording a reproduction run. Keep repository paths short on Windows for nested scenario outputs.

Core RF uses base seed 1, fold seed offsets and model seed offsets stored in its script and run record. Downsampling uses base seed 20260930; marker resampling uses 20261001. Individual design and replicate seeds are retained with their outputs.

Upstream fastp 0.23.4, OMA 2.7.0 and read2tree 2.0.1 generated the supplied ortholog alignments, as described in the manuscript. That upstream processing is outside this R repository and is not required merely to inspect the archived SNP outputs.

Official installation references: [Bioconductor installation](https://bioconductor.org/install/) and [Bioconductor 3.21 release](https://bioconductor.org/news/bioc_3_21_release/).
