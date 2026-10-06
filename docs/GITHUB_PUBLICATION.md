# GitHub publication records

Release [v1.03](https://github.com/PESTFLY/PESTFLY_2026_JPS_Bd_genomic_assignment_workflow/releases/tag/v1.03) was published on 6 October 2026 from commit `f397f8e0a269fa089140eff22a00618d3583dd3c`. The repository is public at <https://github.com/PESTFLY/PESTFLY_2026_JPS_Bd_genomic_assignment_workflow>.

## Published release

The Step 8 checkpoint was prepared on 5 October 2026 and uploaded through branch `publication_v1_03`. [Pull request #1](https://github.com/PESTFLY/PESTFLY_2026_JPS_Bd_genomic_assignment_workflow/pull/1) merged the preparation commit `7cb8661f91c767d3a06cb9db8bf56fa4b8e3eaf3` into `main` on 6 October 2026.

The release tag `v1.03` identifies the merged snapshot, with tree `4e0873284333879f98c6d88b7381d9ebf5d46f2b`. All 1,319 checkpoint files matched their GitHub blob contents and sizes exactly. The three additional Step 00 conversion records matched the original supplied `results.zip` files exactly. The source file manifest on `main` also records those three conversion files. Ten empty directory placeholders have no scientific content.

The earlier release [v1.02](https://github.com/PESTFLY/PESTFLY_2026_JPS_Bd_genomic_assignment_workflow/releases/tag/v1.02), published on 19 May 2026, identifies commit `62e9396205a65b33d3c8118485b4122c5c1bccbe` and had no attached assets when inspected.

## Verified supplementary assets

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

`CITATION.cff` on `main` records version `v1.03`, publication date `2026-10-06`, the release URL and workflow DOI [10.5281/zenodo.23187057](https://doi.org/10.5281/zenodo.23187057). The current software authors are Massimiliano Virgilio and Lore Esselens. Manuscript authorship is maintained separately. The final published manuscript citation can be recorded when available.

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

On 6 October 2026, the author requested a correction to the software creator list. The citation on `main` now lists Massimiliano Virgilio and Lore Esselens. The manuscript author list is unaffected. The corresponding creator metadata correction on Zenodo remains pending.

A corrected Supplementary File S6 has been prepared with the current software citation, restricted raw data access notices, current public repository address and a dated assembly metadata revision. It also incorporates the author's current conceptual guide edits. All R scripts, statistical outputs and figures retain their original bytes. The reader guide adds the raw data, workflow and repository links on page 9; its first eight pages are unchanged in text and rendered appearance. All 1,208 archive manifest entries pass their size and SHA256 checks. The corrected ZIP has 81183349 bytes and SHA256 `dcc9cf7efadb4fd5f5ab347abb8d6222794d0705287d7407d8ffdd4bee35f87b`. It has not replaced the published v1.03 asset.

The v1.03 tag and its original supplementary assets remain the published snapshot, including their original citation metadata. An updated release is required to distribute corrected citation files through the GitHub and Zenodo source archives. Historical Git records have not been rewritten. Zenodo permits creator metadata edits on the existing DOI; see [edit published records](https://help.zenodo.org/docs/deposit/manage-records/#edit).

## Documentation used

* [Git attributes](https://git-scm.com/docs/gitattributes)
* [Managing GitHub releases](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)
* [Citation files on GitHub](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files)

## Repository rename

On 6 October 2026, the repository was renamed to `PESTFLY_2026_JPS_Bd_genomic_assignment_workflow`. The GitHub repository identifier remains `1241595394`. Current repository, release and pull request links use the new address. The published v1.03 tag and its attached files retain their original contents. The S6 builder preserves maintained software creator, data access and repository metadata while recalculating its inventory fields. The prepared S6 archive was rebuilt with this builder and compared with the preceding prepared archive; every analysis script and result file was identical.
