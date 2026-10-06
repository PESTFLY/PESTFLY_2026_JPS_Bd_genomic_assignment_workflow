# Conceptual guide

This guide summarises the interpretation used in the manuscript. The manuscript provides the full methods, equations and biological discussion.

## What an assignment means

Geographic assignment means identifying the represented reference class with greatest genomic affinity to a query. Application requires prior identification as *Bactrocera dorsalis*; uncertain specimens need independent morphological and molecular identification.

PESTFLY is therefore a closed set classifier: support identifies the strongest affinity among represented classes but cannot establish that the true source was represented.

## Benchmark and markers

The fixed snapshot contains 330 references and 22 Belgian queries. P1 compares Africa with Asia; P2 and P3 compare subregions within Africa and Asia. Ten references labelled Other are reserved for the geographic challenge. The 40 Reunion and five Mauritius references use Asia and Southeast Asia coding based on the manuscript's lineage interpretation.

Each specimen supplies one consensus nucleotide per site. Unambiguous reference and alternate alleles are encoded as 0 and 1; ambiguous states are missing. Hudson FST ranks reference SNPs, retaining the highest scoring SNP per ortholog. Mutual information provides a diagnostic comparison.

## Assignment and reporting

Step 04 uses smoothed reference allele frequencies and equal class priors. K denotes the number of highest ranked SNPs in a cumulative panel. Agreement across K measures convergence; alternative marker sets test sensitivity to marker choice separately.

Six criteria determine the finest supported reporting level. Random Forest provides separate corroboration, and declared commodity origin supplies an external consistency check. See [reporting criteria](REPORTING_CRITERIA.md) for thresholds and branch decisions. The thresholds are operational settings that require recalibration as the reference database expands.

## Validation scope

| Analysis | What it evaluates |
| :--- | :--- |
| Individual leave one out | Prediction after withholding a specimen from allele frequency estimation, using fixed SNP panels |
| Country and site holdouts | Geographic generalisation with filtering, ranking and frequency estimation repeated on training references within the retained candidate SNP set |
| Out of scope challenge | Assignment of ten Other specimens against Africa and Asia candidate classes |
| Reference downsampling | Sensitivity to balanced reference size, with training based filtering and reranking within retained candidates |
| Marker resampling | Sensitivity to alternative rank matched marker sets |
| Coding and exclusion scenarios | Sensitivity to analytical class definitions and retained references |

Grouped validation and downsampling retain the original one SNP per ortholog candidate set. Paired marker sets are disjoint within each partition; different partitions can overlap.

Read accuracy together with evaluable sample counts and reportability. Query agreement compares a call with its saved baseline, rather than a known geographic source. Module guides describe the files needed to inspect each analysis.
