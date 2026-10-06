# Input data

Public sample metadata are supplied as `data/000_input_data/metadata.xlsx`.

The PHYLIP ortholog alignments required by Step 00 belong in `data/000_input_data/phy/`. They are distributed upon request from Massimiliano Virgilio at <massimiliano.virgilio@africamuseum.be>. See [alignment access](../ALIGNMENT_ACCESS.md). The revised manuscript identifies the deposit as <https://doi.org/10.5281/zenodo.20283931>; its metadata remain to be verified.

Step 00 creates FASTA alignments in `results/00_fasta/`. This packaging draft retains conversion summaries and run records, but excludes generated FASTA alignments.

The metadata schema is documented in docs/METADATA.md and runtime requirements in docs/INSTALL.md.

Public metadata and derived SNP resources follow the content licence in [LICENSING.md](../LICENSING.md), except where a separate notice applies. This public licence does not cover the separately distributed input alignments.
