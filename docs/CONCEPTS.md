# Conceptual guide

## What an assignment means

PESTFLY estimates the relative genomic affinity of a previously identified *Bactrocera dorsalis* specimen to represented analytical classes. Reference coverage, class definitions and the observed consensus allele calls determine those comparisons.

This is a closed set classifier. Posterior, gap and stability criteria therefore cannot show that a source absent from the candidate classes has been rejected. Correct morphological species identification is a prerequisite; this workflow has no validated taxonomic rejection rule for other members of the *B. dorsalis* complex or putative hybrids.

## Benchmark and analytical classes

The fixed snapshot has 352 specimens: 330 references and 22 Belgian queries. P1 uses 105 African and 215 Asian references. The remaining ten references labelled Other are reserved for the geographic challenge. P2 uses four African subregion classes and P3 uses three Asian subregion classes.

The 40 Reunion and five Mauritius references are coded as Asia and Southeast Asia in the baseline because of the lineage interpretation adopted in the manuscript. These labels express analytical ancestry groupings rather than collection geography. Alternative coding and exclusion are tested separately.

## Consensus calls and SNP ranking

Each specimen contributes one reconstructed nucleotide state per site. Unambiguous reference and alternate states are encoded as 0 and 1; ambiguous and IUPAC states are missing. These are consensus or pseudohaploid observations, with no model of within specimen diploid heterozygosity.

Reference only discovery requires biallelic A, C, G or T sites, missingness no greater than 0.20 and minor allele count at least two. Each panel retains the highest scoring SNP per ortholog.

For two groups, Hudson FST is calculated from alternate allele frequencies p1 and p2 and numbers n1 and n2 of nonmissing consensus observations:

$$
F_{ST} = \frac{(p_1-p_2)^2 - p_1(1-p_1)/(n_1-1) - p_2(1-p_2)/(n_2-1)}{p_1(1-p_2)+p_2(1-p_1)}.
$$

Insufficient counts, invalid or zero denominators and nonfinite results are treated as missing; finite negative estimates are retained. P1 uses its binary contrast. P2 and P3 use the maximum pairwise subregion score. The score ranks candidate assignment features.

Mutual information is calculated in bits as a reference only diagnostic comparison. It does not select or reorder markers, or enter the assignment likelihood.

## Likelihood and convergence

Step 04 uses smoothed allele frequencies with pseudocount 0.5, an operational allele flip allowance epsilon of 0.02 and equal class priors. Contributions from usable loci are combined under conditional independence. Resulting posterior values describe the model's comparison among represented classes, rather than calibrated probabilities that an actual source country has been identified.

Top K means the K highest ranked SNPs. Cumulative panels test whether a prediction converges as more markers are considered. As markers nested a separate marker resampling analysis examines nonnested alternative panels. Six criteria determine reportability, as described in [reporting criteria](REPORTING_CRITERIA.md). Their thresholds are purely operational, and are not estimated as error rate cutoffs from the 22 queries.

## What each validation tests

| Analysis | What is held out or changed | What remains conditional |
| :--- | :--- | :--- |
| Individual leave one out | Focal reference removed from allele frequency estimation | Its contribution to original SNP discovery and ranking |
| Country and site holdout | Whole geographic unit removed; training only filtering, ranking and frequency estimation repeated | The previously discovered one SNP per ortholog candidate set |
| Random Forest CV | Reference specimens held out during model fitting and imputation | The same fixed SNP resources used by the main workflow |
| Out of scope challenge | Ten Other references treated as unknown | Only Africa and Asia candidate classes are available |
| Reference downsampling | Balanced training references and independent test sets resampled | Previously retained candidate SNPs |
| Marker resampling | Complementary rank matched marker sets substituted | The existing Hudson ranked candidate pool |
| Coding and exclusion scenarios | Reference class definitions or retained specimens changed | The upstream alignments and specified scenario designs |

Alternative marker A and B sets are disjoint within each paired partition. Different partitions can share markers. Thirty alternative sets do not constitute thirty independent discovery datasets, and their fixed size evaluations do not apply the cumulative multi K convergence criteria.

Raw accuracy, accuracy among evaluable specimens and the proportion receiving reportable calls answer different questions. Include denominators and uncertainty rates when comparing validation designs.

## Reporting and commodity comparisons

RF is a second algorithm on related data. Its agreement or disagreement is corroboration, not a reporting criterion and not an override of Step 04. Declared commodity country is external metadata used after fitting for a broad consistency check. It can differ from the true biological source. Neither agreement nor disagreement verifies a transport pathway.

Step 07 explicitly reports six criterion PASS or FAIL decisions. Its verified native R report preserves the original supported classes and distinguishes an unevaluated subregion from a tested failure. RF and commodity comparisons remain separate. The reporting checks and formatter passed in R 4.5.1; see [reporting criteria](REPORTING_CRITERIA.md) and [native R verification](NATIVE_R_VERIFICATION.md).
