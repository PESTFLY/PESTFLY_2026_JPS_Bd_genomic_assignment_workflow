# Six reporting criteria

The revised manuscript specifies six criteria for each evaluated resolution. The largest usable K supplies the final candidate group and posterior support; only K values meeting the query SNP requirement enter convergence calculations.

| Criterion | Observed Step 04 field suffix | P1 threshold | P2 or P3 threshold |
| :--- | :--- | :--- | :--- |
| Usable query SNPs at final K | `maxK_nsnps` | At least 200 | At least 300 |
| Posterior of best represented class | `maxK_post` | At least 0.85 | At least 0.85 |
| Posterior gap between first and second class | `maxK_gap` | At least 0.20 | At least 0.25 |
| Agreement with final group across usable K | `agreement` | At least 0.90 | At least 0.90 |
| Agreement across largest half of usable K | `tail_agreement` | 1.00 | 1.00 |
| Number of usable K values | `nK_usable_for_stability` | At least 6 | At least 6 |

The suffix is prefixed by `macroregion_` or `subregion_` in `final_assignment.tsv`. The largest half uses the ceiling of the number of usable K values times 0.5. The high posterior category starts at 0.95; moderate is 0.85 to below 0.95, subject to the other reporting criteria. High and Moderate are model support labels, not independently calibrated source identification probabilities.

The standard evaluated K grid is 1, 2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000 and 10000, with all available SNPs added. Small K outputs are retained for diagnosis but do not count toward stability unless the usable SNP requirement is met.

PASS means all six criteria are met at that resolution. FAIL means at least one is unmet. When the macroregion passes and its conditional subregion fails, report only the macroregion. When the macroregion fails, the sample remains unresolved and the subregion is not evaluated. An unevaluated branch must be distinguished from a tested FAIL. RF remains separate from the six criterion count.

## Audit against the supplied archive

Using the archived Step 04 values and parameter thresholds, the documentation audit reconstructed 19 macroregion PASS outcomes and three FAIL outcomes. Each failed adult met five criteria and failed only the usable K count. Twelve evaluated subregions passed, seven failed, and the three adults with failed macroregions had no evaluated subregion branch. All tested decisions agree with Step 04 confidence eligibility.

The diagnostic table is `docs/reporting_criteria_diagnostic.tsv`. It contains observed values, thresholds, individual results and criterion counts for each specimen and resolution. This is an independently derived documentation check against saved outputs; it is not a rerun or the regenerated final authority report.

## Reconciled reporting layer

The updated Step 07 formatter applies these six rules to saved metrics and reads thresholds from the matching Step 04 parameter table. `reporting_criteria.R` implements the rules in base R. The formatter stops if criterion eligibility or the High or Moderate category differs from Step 04. Evaluated missing metrics fail their criterion, while malformed fields, invalid proportions, inconsistent branch panels and duplicate specimen identifiers stop generation.

The current workbook includes a `Criteria` sheet with 264 rows, six for every specimen and resolution, and a `Reporting_Thresholds` sheet. The simplified report shows both statuses, criterion counts and all unmet criterion identifiers. Unevaluated subregions have missing counts and observed values. They are not counted as tested failures. Original colour categories are retained only with historical reports in `results/07_authority_report/legacy_colour_report/`.

The current workbook, TSV tables and RDS objects are native exports from the updated Step 07 formatter, executed on 5 October 2026 in R 4.5.1. The base R reporting checks passed. All 22 reported classes, resolutions and confidence labels are unchanged, and every criterion result matches the independent documentation diagnostic. The earlier reconstruction is retained under `results/07_authority_report/reconstruction_checkpoint/`. See [native R verification](NATIVE_R_VERIFICATION.md) and the [Step 07 guide](../steps/07_authority_facing_report/README.md).
