# GitHub publication records

Release [v1.04](https://github.com/PESTFLY/PESTFLY_2026_JPS_Bd_genomic_assignment_workflow/releases/tag/v1.04) was published on 6 October 2026 from commit `3700ac653d08ae7a4e016c41777d231294205b5e`. The repository is public at <https://github.com/PESTFLY/PESTFLY_2026_JPS_Bd_genomic_assignment_workflow>. The original v1.03 snapshot and assets are retained.

## Published release v1.03

The Step 8 checkpoint was prepared on 5 October 2026 and uploaded through branch `publication_v1_03`. [Pull request #1](https://github.com/PESTFLY/PESTFLY_2026_JPS_Bd_genomic_assignment_workflow/pull/1) merged the preparation commit `7cb8661f91c767d3a06cb9db8bf56fa4b8e3eaf3` into `main` on 6 October 2026.

The release tag `v1.03` identifies the merged snapshot, with tree `4e0873284333879f98c6d88b7381d9ebf5d46f2b`. All 1,319 checkpoint files matched their GitHub blob contents and sizes exactly. The three additional Step 00 conversion records matched the original supplied `results.zip` files exactly. The source file manifest on `main` also records those three conversion files. Ten empty directory placeholders have no scientific content.

The earlier release [v1.02](https://github.com/PESTFLY/PESTFLY_2026_JPS_Bd_genomic_assignment_workflow/releases/tag/v1.02), published on 19 May 2026, identifies commit `62e9396205a65b33d3c8118485b4122c5c1bccbe` and had no attached assets when inspected.

## Verified v1.03 supplementary assets

| Release asset | Bytes | SHA256 |
| :--- | ---: | :--- |
| `PESTFLY_Supplementary_File_S5.xlsx` | 67636 | `7acba0de94d06097afe723381185b9e9e92f7ddc1dc2be052ee77aaa81fb1d96` |
| `PESTFLY_Supplementary_File_S6.zip` | 81181907 | `2024ec3275c9b8eed8441b3a694a8817a6eae20e69bf2416c20430499b8eb9c1` |

Both assets were uploaded successfully, and GitHub's recorded SHA256 digests match the verified publication files. They retain the checkpoint bytes. Citation and publication documentation on `main` record the completed release; the tag and supplementary assets identify the original published snapshot.

Code uses MIT. Project documentation, public metadata and derived results use CC BY 4.0, except where separate notices apply. Copyright holders are Massimiliano Virgilio and Wannes Dermauw. Input alignments are distributed separately upon request and are excluded from the public licence grants.

## File integrity

The root `.gitattributes` preserves recorded file bytes through Git checkin and checkout, including on Windows. Text differences remain enabled for source files and documentation.

Before upload, the checkpoint was checked using a separate local Git index with `core.autocrlf=true`. Every staged blob was compared with its source bytes, and every exported checkout file was compared with the package manifest. After upload, the complete GitHub tree was compared with the verified checkpoint and the original Step 00 conversion records. The citation file and documentation have their current hashes in `docs/source_file_manifest.tsv`.

The native R implementation checks, six criterion reporting checks and Step 07 formatter passed in the supplied R session. The original wrapper's descriptive label mismatch and the independent report comparison are documented in [native R verification](NATIVE_R_VERIFICATION.md). Publication packaging did not recompute statistical models.

## Citation metadata

`CITATION.cff` on `main` records version `v1.04`, published on `2026-10-06`, with its own release URL. Its version specific Zenodo DOI remains unverified, so the earlier v1.03 DOI is not assigned to v1.04. The original v1.03 archive remains identified by [10.5281/zenodo.23187057](https://doi.org/10.5281/zenodo.23187057). The current software authors are Massimiliano Virgilio and Lore Esselens. Manuscript authorship is maintained separately. The final published manuscript citation can be recorded when available.

## Zenodo records

The author supplied a [published Zenodo record](https://zenodo.org/records/23187057) showing version `v1.03`, publication date `6 October 2026`, the matching workflow title and release description, followed by the DOI link [10.5281/zenodo.23187057](https://doi.org/10.5281/zenodo.23187057). This author supplied evidence confirms the workflow citation. Automated public page and DOI registry retrievals remained unavailable.

| Record | Identifier | Verification status |
| :--- | :--- | :--- |
| Workflow archive for v1.03 | `10.5281/zenodo.23187057` | Confirmed from the author supplied published record and DOI link |
| Earlier workflow identifier from the manuscript | `10.5281/zenodo.20289875` | Relationship to the current record remains to be confirmed |
| Restricted raw data deposit | `10.5281/zenodo.20340447` | Public Zenodo metadata confirm restricted file access, title, creator, date and version |
| Earlier alignment identifier from the manuscript | `10.5281/zenodo.20283931` | Relationship to the confirmed raw data deposit remains unconfirmed; processed alignment inputs are available upon request |

For the released workflow, cite the confirmed v1.03 DOI. For the raw data, cite [restricted record 10.5281/zenodo.20340447](https://zenodo.org/records/20340447). Zenodo identifies it as `PESTFLY_diagnostic_snp_raw_data_v1.0_2026-05-22`, published on 22 May 2026, version `v1`, with Massimiliano Virgilio as creator. Its metadata are public, while files require authorised access. Restricted files were not downloaded or inspected during this documentation update. No relationship to the earlier alignment identifier, all versions DOI, or earlier workflow identifier has been inferred.

## Creator metadata correction

On 6 October 2026, the author requested a correction to the software creator list. The citation on `main` now lists Massimiliano Virgilio and Lore Esselens. The manuscript author list is unaffected. The correction to the older v1.03 Zenodo creator metadata remains unverified; the new v1.04 Zenodo metadata also await verification.

The corrected Supplementary File S6 was published with v1.04 with the current software citation, restricted raw data access notices, current public repository address and a dated assembly metadata revision. It includes the author's current conceptual guide edits. All R scripts, statistical outputs and figures retain their original bytes. Page 9 of the reader guide records the prepared v1.04 release and distinguishes its pending DOI from the earlier v1.03 archive; its first eight pages are unchanged in text and rendered appearance. All 1,208 archive manifest entries pass their size and SHA256 checks. The corrected ZIP has 81183572 bytes and SHA256 `14171ee89e2376991d4f9ee708bd2f2f607ab3887ec3404d98931b63e0876e9b`. It has not replaced the published v1.03 asset.

The v1.03 tag and its original supplementary assets remain the published snapshot, including their original citation metadata. The v1.04 GitHub release distributes the corrected citation files. Its Zenodo record still requires verification. Historical Git records have not been rewritten. Zenodo permits creator metadata edits on the existing DOI; see [edit published records](https://help.zenodo.org/docs/deposit/manage-records/#edit).

## Documentation used

* [Git attributes](https://git-scm.com/docs/gitattributes)
* [Managing GitHub releases](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)
* [Citation files on GitHub](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files)

## Repository rename

On 6 October 2026, the repository was renamed to `PESTFLY_2026_JPS_Bd_genomic_assignment_workflow`. The GitHub repository identifier remains `1241595394`. Current repository, release and pull request links use the new address. The published v1.03 tag and its attached files retain their original contents. The S6 builder preserves maintained software creator, data access and repository metadata while recalculating its inventory fields. The prepared S6 archive was rebuilt with this builder and compared with the preceding prepared archive; every analysis script and result file was identical.

## Published release v1.04

Status: **PUBLISHED_ON_GITHUB**. GitHub release identifier: `404885223`. Publication timestamp: `2026-10-06T15:30:34Z`. Tag `v1.04` identifies commit `3700ac653d08ae7a4e016c41777d231294205b5e` and tree `d315857fa04de73b5cfd7327ff5dc9170ba76c04`. The release is neither a draft nor a prerelease. Both uploaded supplementary files match their prepared byte counts and SHA256 digests.

| Published v1.04 asset | Bytes | SHA256 |
| :--- | ---: | :--- |
| `PESTFLY_Supplementary_File_S5.xlsx` | 67636 | `7acba0de94d06097afe723381185b9e9e92f7ddc1dc2be052ee77aaa81fb1d96` |
| `PESTFLY_Supplementary_File_S6.zip` | 81183572 | `14171ee89e2376991d4f9ee708bd2f2f607ab3887ec3404d98931b63e0876e9b` |

The v1.04 source snapshot and S6 contain preparation records created before publication. The tag and uploaded files retain those verified bytes. Current publication status is recorded here and in `external_publication_records.json`. The v1.04 Zenodo record, DOI and creator metadata await verification. The DOI `10.5281/zenodo.23187057` continues to identify the earlier v1.03 archive.
