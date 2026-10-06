# GitHub publication records

Release [v1.03](https://github.com/PESTFLY/PESTFLY_origin_tracing_publication_pipeline/releases/tag/v1.03) was published on 6 October 2026 from commit `f397f8e0a269fa089140eff22a00618d3583dd3c`. The repository is public at <https://github.com/PESTFLY/PESTFLY_origin_tracing_publication_pipeline>.

## Published release

The Step 8 checkpoint was prepared on 5 October 2026 and uploaded through branch `publication_v1_03`. [Pull request #1](https://github.com/PESTFLY/PESTFLY_origin_tracing_publication_pipeline/pull/1) merged the preparation commit `7cb8661f91c767d3a06cb9db8bf56fa4b8e3eaf3` into `main` on 6 October 2026.

The release tag `v1.03` identifies the merged snapshot, with tree `4e0873284333879f98c6d88b7381d9ebf5d46f2b`. All 1,319 checkpoint files matched their GitHub blob contents and sizes exactly. The three additional Step 00 conversion records matched the original supplied `results.zip` files exactly. The source file manifest on `main` also records those three conversion files. Ten empty directory placeholders have no scientific content.

The earlier release [v1.02](https://github.com/PESTFLY/PESTFLY_origin_tracing_publication_pipeline/releases/tag/v1.02), published on 19 May 2026, identifies commit `62e9396205a65b33d3c8118485b4122c5c1bccbe` and had no attached assets when inspected.

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

`CITATION.cff` on `main` records version `v1.03`, publication date `2026-10-06` and the release URL. The existing author list is retained. A confirmed workflow DOI and the final published manuscript citation can be recorded when available.

## Zenodo records still to confirm

| Record | Identifier from manuscript | Current verification status |
| :--- | :--- | :--- |
| Workflow archive | `10.5281/zenodo.20289875` | Public page and metadata retrieval unavailable |
| Alignment deposit | `10.5281/zenodo.20283931` | Public page and metadata retrieval unavailable; alignments remain available upon request |

The earlier lookup results do not establish that the identifiers are invalid or that the deposits are unpublished. Confirm each record's publication status, title, creators and access terms in Zenodo. Check whether the workflow identifier represents a specific version or all versions, and confirm which identifier covers this revised release. The manuscript identifiers remain documented as unverified; no unverified DOI has been added to `CITATION.cff`.

## Documentation used

* [Git attributes](https://git-scm.com/docs/gitattributes)
* [Managing GitHub releases](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)
* [Citation files on GitHub](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files)
