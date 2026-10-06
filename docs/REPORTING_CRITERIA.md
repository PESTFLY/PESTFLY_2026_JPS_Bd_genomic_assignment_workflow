# Six reporting criteria

All six criteria must pass at an evaluated resolution. The final class and support come from the largest usable K; convergence uses only K values meeting the specimen's SNP requirement.

| Criterion | Step 04 field suffix | P1 threshold | P2 or P3 threshold |
| :--- | :--- | :--- | :--- |
| Usable SNPs at final K | `maxK_nsnps` | At least 200 | At least 300 |
| Best class posterior | `maxK_post` | At least 0.85 | At least 0.85 |
| First versus second posterior gap | `maxK_gap` | At least 0.20 | At least 0.25 |
| Agreement across usable K | `agreement` | At least 0.90 | At least 0.90 |
| Agreement in the largest half of usable K | `tail_agreement` | 1.00 | 1.00 |
| Number of usable K values | `nK_usable_for_stability` | At least 6 | At least 6 |

In `final_assignment.tsv`, suffixes have the prefix `macroregion_` or `subregion_`. The tail contains `ceiling(n_usable_K * 0.5)` values. Posterior support is High at 0.95 or above and Moderate from 0.85 to below 0.95, with all other criteria passing.

The standard K grid is 1, 2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000 and 10000, followed by all available SNPs. Small K outputs remain available for inspection.

## Reporting decisions

| Macroregion | Conditional subregion | Report |
| :--- | :--- | :--- |
| PASS | PASS | Subregion affinity |
| PASS | FAIL | Macroregion affinity |
| FAIL | NOT_EVALUATED | Unresolved |

Step 07 reads thresholds from the matching `step4_params.tsv` and checks consistency with saved Step 04 decisions. Missing metrics fail individual criteria; malformed or contradictory inputs stop report generation. RF and commodity comparisons are separate from the criterion count.

The benchmark has 19 macroregion PASS and three FAIL outcomes, plus 12 subregion PASS, seven FAIL and three NOT_EVALUATED outcomes. Each unresolved adult fails only the usable K count. The workbook's `Criteria` sheet records all 264 specimen and resolution criterion rows; `Reporting_Thresholds` records their settings.

See the [Step 07 guide](../steps/07_authority_facing_report/README.md) for commands and files, [native R verification](NATIVE_R_VERIFICATION.md) for checks, and the manuscript for the rationale and interpretation.
