# GitHub publication preparation

This is the Step 8 publication checkpoint prepared on 5 October 2026. The existing public repository is <https://github.com/PESTFLY/PESTFLY_origin_tracing_publication_pipeline>. Its inspected main commit is `1b482b87f5d5611f87b413abdc555272422ffe24`.

## Release identity

The current published release is [v1.02](https://github.com/PESTFLY/PESTFLY_origin_tracing_publication_pipeline/releases/tag/v1.02), dated 19 May 2026. It points to commit `62e9396205a65b33d3c8118485b4122c5c1bccbe` and has no attached assets. The next proposed release tag is `v1.03`, following the existing naming pattern. The tag has not been created. The proposed preparation branch is `publication_v1_03`.

Code uses MIT. Project documentation and public results use CC BY 4.0, except where separate notices apply. Copyright holders are Massimiliano Virgilio and Wannes Dermauw. Input alignments remain distributed separately upon request.

## Preserve archived file bytes

The root `.gitattributes` prevents Git from converting line endings. This preserves the SHA256 records for both code and archived outputs even on a Windows checkout. Text differences remain enabled for source files and documentation. Keep this file when copying the checkpoint into a working repository.

The package was checked using a separate local Git index with `core.autocrlf=true`. Every staged blob was compared with its source bytes, and every exported checkout file was compared with the package manifest. No scientific model was rerun. This local check made no public GitHub changes.

## Prepare a branch on Windows

1. Extract `PESTFLY_repository_step8.zip` outside the existing repository. Open the existing repository in GitHub Desktop.
2. Select `main`, fetch the current repository state, then create the proposed `publication_v1_03` branch from it.
3. Copy the contents of the extracted `PESTFLY_repository_step8` folder into the repository root. Keep the existing Git history and local metadata directory. Copy the root files as well as `data`, `docs`, `environment`, `results`, `steps`, `supplementary` and `tools`.
4. Review the changes. The compared public main has nine empty placeholder files that are unnecessary; removing those placeholders is optional. Keep any local alignments outside the public commit. The prepared checkpoint contains no PHYLIP or FASTA alignments.
5. Commit the reviewed update and publish the preparation branch. The comparison inventory records the prepared files and the inspected main commit. Recheck the main commit if it has changed since this checkpoint.

The recorded comparison is a complete package comparison. Local files not included in this checkpoint should be reviewed individually rather than removed automatically.

## Draft the release after the repository update

Create a release draft targeting the updated commit. Use the proposed `v1.03` tag only after choosing the final release identity. Attach the verified `PESTFLY_Supplementary_File_S5.xlsx` and `PESTFLY_Supplementary_File_S6.zip`. The preparation folder contains the draft release text. The final publication date and release version can then be added to `CITATION.cff`.

## Zenodo records still to confirm

| Record | Identifier from manuscript | Current verification status |
| :--- | :--- | :--- |
| Workflow archive | `10.5281/zenodo.20289875` | Public page and metadata retrieval unavailable |
| Alignment deposit | `10.5281/zenodo.20283931` | Public page and metadata retrieval unavailable; alignments remain available upon request |

These lookup results do not establish that the identifiers are invalid or that the deposits are unpublished. Confirm each record's publication status, title, creators and access terms in Zenodo. Check whether the workflow identifier represents a specific version or all versions. A version identifier for the earlier release should not be presented as the identifier for this revised release. The existing manuscript identifiers remain documented, but no unverified DOI has been added to `CITATION.cff`.

## Documentation used

* [Git attributes](https://git-scm.com/docs/gitattributes)
* [Managing branches in GitHub Desktop](https://docs.github.com/en/desktop/making-changes-in-a-branch/managing-branches-in-github-desktop)
* [Managing GitHub releases](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository)
* [GitHub file size limits](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)
