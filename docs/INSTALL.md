# Installation and reproduction

## Recorded benchmark versions

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

These versions match the manuscript and retained session records; sources are listed in `runtime_version_evidence.tsv`. Base R supplies `parallel`, `stats`, `utils`, `tools` and `grDevices`. The archived RF run used ranger and recorded seeds; retain that backend for reproduction.

## Install and check

In R, using configured CRAN repositories and Bioconductor:

```r
install.packages(c("data.table", "optparse", "readxl", "openxlsx",
                   "ranger", "randomForest", "BiocManager"))
BiocManager::install("Biostrings", version = "3.21",
                     ask = FALSE, update = FALSE)
```

This example installs available packages, rather than restoring a locked environment. Compare versions with `environment/benchmark_versions.tsv`; a complete dependency lock is not supplied. These setup commands have not been executed in the packaging environment. Official references: [Bioconductor installation](https://bioconductor.org/install/) and [release 3.21](https://bioconductor.org/news/bioc_3_21_release/).

From the repository root:

```bash
Rscript environment/check_environment.R
Rscript environment/check_environment.R --strict
```

Both checks fail for missing required packages; strict mode also fails for version differences. The checker has received syntax screening.

## Choose a reproduction route

| Route | Required resources | Modules |
| :--- | :--- | :--- |
| Inspect outputs | Supplied metadata and archived files | Tables, figures, reports and run records |
| Rerun from SNP panels | R packages, archived QC and SNP resources | Core 03 to 07; supplementary 01, 02, 03, 05 and 06 |
| Complete core workflow | Processed PHYLIP inputs, metadata and R packages | Core 00 to 07 |
| Rebuild scenario panels | Regenerated FASTA and baseline resources | Supplementary 04 and 07 |

Request alignment inputs as described in [alignment access](../ALIGNMENT_ACCESS.md) and place PHYLIP files in `data/000_input_data/phy/`. The raw data deposit has [restricted file access](https://doi.org/10.5281/zenodo.20340447).

```bash
Rscript run_public_pipeline.R
Rscript run_validation_pipeline.R
```

Use a separate copy for recomputation because runners use defaults and can replace archived outputs. Individual module guides list inputs, outputs and options. [cli_parameters.tsv](cli_parameters.tsv) gives all 234 documented options; `Rscript <module>/run.R --help` displays the script's full help. Saved parameter and run records describe the actual benchmark settings.

Step 07 can regenerate reports from saved Step 04 assignments and parameters without sequence inputs. See [native R verification](NATIVE_R_VERIFICATION.md) for the completed checks and original checker mismatch.

## Workers and seeds

Core alignment modules use `--cores`; zero requests detected cores minus one. Grouped validation, downsampling and marker resampling use `--n_cores`, with an automatic cap of eight. Specify workers explicitly when recording a run and use short paths on Windows for nested scenarios.

Base seeds are 1 for RF, 20260930 for downsampling and 20261001 for marker resampling. Derived seeds are retained in scripts and run records.

Upstream fastp 0.23.4, OMA 2.7.0 and read2tree 2.0.1 processing is described in the manuscript and lies outside this R workflow.
