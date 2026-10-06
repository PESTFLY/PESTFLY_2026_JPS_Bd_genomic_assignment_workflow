# Supplementary File S6

## Validation and sensitivity analyses for the PESTFLY benchmark

This archive supports the population genomic analysis of *Bactrocera dorsalis* using previously generated sequence data and completed statistical outputs. It contains seven supplementary analyses with their R scripts, parameters, summary tables, figures, specimen and replicate outputs, plus complete baseline query trajectories across the evaluated SNP subsets.

Open [S6_reader_guide.pdf](S6_reader_guide.pdf) for the design, findings and interpretation of each analysis. [source_map.tsv](source_map.tsv) maps sections to scripts and key outputs. All paths below are relative to the extracted `PESTFLY_S6` directory. Script defaults require running from that directory.

## Analysis index

Result directories are under `results/validation/`.

| Section | Analysis | Result directory |
| :--- | :--- | :--- |
| S6.1 | Geographic out of scope challenge | `01_out_of_scope_challenge` |
| S6.2 | Grouped geographic cross validation | `02_grouped_geographic_cv` |
| S6.3 | Reference set audit | `03_reference_set_audit` |
| S6.4 | Mascarene coding sensitivity | `04_mascarene_coding_sensitivity` |
| S6.5 | Balanced reference downsampling | `05_reference_downsampling` |
| S6.6 | Alternative marker set resampling | `06_marker_resampling` |
| S6.7 | Reference exclusion sensitivity | `07_exclusion` |

The matching scripts are in `steps/validation/`. Analysis 07 uses the shorter output directory `07_exclusion` for compatibility with nested Windows paths. Complete scenario outputs for analyses 04 and 07 are retained under `scenarios/`. Their unique cross scenario summaries remain in `return_bundle/`. For analyses 05 and 06, use the main result directory for summaries and figures.

## Benchmark and interpretation

The snapshot contains 352 specimens: 330 references and 22 queries. P1 compares 105 African with 215 Asian references; ten references from Hawaii and Papua New Guinea are outside these represented macroregion classes. P2 compares four African subregions and P3 three Asian subregions.

Assignments indicate genomic affinity to represented classes. They do not establish actual collection origin or a transport pathway. In the closed set challenge, all ten samples from geographically unrepresented locations received confident Asian affinities. Confidence within the represented classes therefore cannot establish that the actual source is represented.

Raw accuracy is the fraction of evaluable reference specimens whose top call matches their recorded class. Reportability is the fraction meeting the applicable reporting criteria. Accuracy among reportable calls uses only accepted calls in its denominator. Query agreement measures agreement with a saved baseline call; it is not accuracy against known query origin. Grouped folds without the held out true class in training are recorded as non evaluable, rather than counted as errors or correct calls.

Fixed panel leave one out validation refits allele frequencies after withholding each focal specimen, conditional on the full reference marker panel. Grouped geographic validation and downsampling filter and rerank within the retained Step 02 candidate SNPs using training references. They do not rediscover an alternative SNP from each original ortholog. Marker resampling uses complementary sets within each partition of the existing ranking; sets from different partitions can overlap.

## Complete query trajectories

The three baseline raw tables in `results/04_origin_assignment/` contain 656 rows over 16 K values per evaluated branch:

| File | Queries | Rows |
| :--- | :--- | :--- |
| `P1_macroregion_raw.tsv` | 22 | 352 |
| `P2_africa_subregion_raw.tsv` | 12 | 192 |
| `P3_asia_subregion_raw.tsv` | 7 | 112 |

Each row retains the top and second class, their posterior support, posterior gap, usable SNP count, confidence and status. Low K rows and unavailable values are retained. P2 and P3 are conditional on an accepted macroregion call; three failed macroregion queries have no evaluated subregion branch. The stability tables, final assignments, K grid and parameter records are in the same directory.

The six baseline reporting requirements are the minimum usable SNP count, posterior support, posterior gap, global agreement across usable K values, tail agreement and the number of usable K values. The archived thresholds are documented in `docs/REPORTING_CRITERIA.md`. Random Forest corroboration remains a separate comparison.

## Inspection and reproduction

TSV tables, gzipped TSV tables, PDF figures and Excel summaries can be inspected directly. RDS files retain the original R objects. Core metadata, QC records, SNP panels and scripts are included so the archive can be inspected and the panel based analyses can be rerun without a separate repository checkout.

See `docs/INSTALL.md`, `docs/CONCEPTS.md`, and the module READMEs for dependencies and scope. Recorded benchmark versions include R 4.5.1. A complete dependency lock has not been reconstructed.

The PHYLIP and FASTA alignments are distributed separately upon request from Massimiliano Virgilio at <massimiliano.virgilio@africamuseum.be>. Restoring alignments is required for full core reproduction and scenario panel rebuilding in analyses 04 and 07. The large text export `results/01_qc/sample_missingness.tsv` is omitted; its complete RDS counterpart is included. The manuscript identifies the alignment deposit as <https://doi.org/10.5281/zenodo.20283931>. External release records remain to be verified.

After restoring required inputs, the separate supplementary runner uses all seven modules' defaults:

```bash
Rscript run_validation_pipeline.R
```

Use a separate extracted copy for recomputation because default outputs can replace archived files. Run one module directly to restrict computation to that analysis. Modules 01, 02, 03, 05 and 06 use supplied QC and SNP resources. Modules 04 and 07 also require regenerated FASTA.

Implementation checks can be run without sequence inputs or third party R packages:

```bash
Rscript tools/check_implementation.R
Rscript steps/07_authority_facing_report/check_reporting.R
```

The isolated implementation wrapper parses the original scripts and evaluates only their named Hudson FST and MI functions. It passed in native R 4.5.1 on 5 October 2026, alongside the six criterion reporting checks and updated Step 07 formatter. Current Step 07 workbook, TSV tables, RDS objects and run record are the verified native R outputs. The earlier reconstruction and original colour reports remain available as provenance. `docs/NATIVE_R_VERIFICATION.md` records the original checker label mismatch and the completed independent export review. `Rscript verify_release.R` runs the corrected checker without refitting models. No statistical model was rerun for this S6 assembly.

## Integrity and provenance

`archive_manifest.tsv` records the size and SHA256 of every archive file except itself. Verify an extracted archive with Python 3 using only its standard library:

```bash
python3 verify_archive.py
```

Only byte identical `return_bundle/` exports with a retained counterpart in the same analysis have been omitted. [duplicate_exports.tsv](duplicate_exports.tsv) maps every omitted copy to its retained file and hash. All scientifically distinct outputs, including per specimen and per replicate tables, scenario records, RDS objects and figures, are included.

Original machine paths, dates and early internal analysis labels remain in historical records. Adapted TSV indexes use the extracted archive layout; `docs/location_index_changes.tsv` records those earlier packaging changes. Six summary text files received presentation changes to describe their analytical purpose and the complete archive. Their numerical lines are unchanged; `docs/textual_output_changes.tsv` records the exact text replacements and hashes. Original records remain available in the source revision archive and earlier repository checkpoints.

`assembly_info.json` records packaging status. Code uses MIT; project documentation and public metadata and derived results use CC BY 4.0, except where a separate notice applies. Copyright (c) 2026 Massimiliano Virgilio and Wannes Dermauw. See `LICENSE`, `LICENSE_CONTENT.txt`, `LICENSING.md` and `ALIGNMENT_ACCESS.md` at the archive root. The input alignments are outside these public licence grants and remain available upon request. External release records remain to be verified.
