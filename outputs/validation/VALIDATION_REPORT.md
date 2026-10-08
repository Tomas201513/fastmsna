# Validation report: fastmsna vs analysistools

Generated: 2026-10-08 01:23 | R 4.5.1 | analysistools 0.0.0.903 | srvyr 1.3.1 | survey 4.5

Rows are matched on `analysis_key`. Values are equal when both are missing or
`|old - new| <= 1e-8 * max(1, |old|)`. `max_rel_diff` is the largest relative
difference over all compared values (floating-point noise only).

| scenario | LOA rows | result rows old / new | value mismatches | same row order | failed rows old / new (same) | max rel diff | old (s) | new (s) | equivalent |
|---|---|---|---|---|---|---|---|---|---|
| template_loa | 13 | 143 / 143 | 0 | TRUE | 0 / 0 (TRUE) | 2.3e-14 | 5.13 | 0.45 | **TRUE** |
| template_loa_with_ratio | 15 | 147 / 147 | 0 | TRUE | 0 / 0 (TRUE) | 2.3e-14 | 4.73 | 0.77 | **TRUE** |
| template_no_loa | 21 | 399 / 399 | 0 | TRUE | 0 / 0 (TRUE) | 2.3e-14 | 8.34 | 0.2 | **TRUE** |
| A_2scs_strata_weights | 171 | 2830 / 2830 | 0 | TRUE | 1 / 1 (TRUE) | 4.2e-13 | 121.95 | 1.3 | **TRUE** |
| G_ratio_documented_filters | 171 | 2802 / 2802 | 0 | TRUE | 1 / 1 (TRUE) | 4.2e-13 | - | - | **TRUE** |
| B_clustered | 171 | 2830 / 2830 | 0 | TRUE | 1 / 1 (TRUE) | 1.4e-13 | 122 | 1.14 | **TRUE** |
| C_quota_unweighted | 171 | 2830 / 2830 | 0 | TRUE | 1 / 1 (TRUE) | 5.1e-14 | 67.83 | 1.02 | **TRUE** |
| D_lonely_fail | 171 | 93 / 93 | 0 | TRUE | 156 / 156 (TRUE) | 0.0e+00 | 19.77 | 1.36 | **TRUE** |
| E_lonely_remove | 171 | 2830 / 2830 | 0 | TRUE | 1 / 1 (TRUE) | 7.8e-15 | 118.5 | 0.97 | **TRUE** |
| E_lonely_certainty | 171 | 2830 / 2830 | 0 | TRUE | 1 / 1 (TRUE) | 7.8e-15 | 122.66 | 0.91 | **TRUE** |
| E_lonely_average | 171 | 2830 / 2830 | 0 | TRUE | 1 / 1 (TRUE) | 5.6e-13 | 119.48 | 0.83 | **TRUE** |
| E_adjust_no_domain | 171 | 2830 / 2830 | 0 | TRUE | 1 / 1 (TRUE) | 1.4e-13 | 105.68 | 0.79 | **TRUE** |
| F_fpc | 171 | 2830 / 2830 | 0 | TRUE | 1 / 1 (TRUE) | 1.1e-13 | 104.39 | 0.84 | **TRUE** |

Scenario descriptions:

- **template_loa**: analysistools template data and LOA (sm separator '/')
- **template_loa_with_ratio**: analysistools template data and LOA (sm separator '/')
- **template_no_loa**: template data, automatic LOA, group_var = c('admin1', 'admin1, admin2')
- **A_2scs_strata_weights**: weights + strata, ids = ~1 (NGA 2SCS Step 2 design), lonely.psu = adjust + domain lonely [ratio_legacy = TRUE]
- **G_ratio_documented_filters**: scenario A with ratio_legacy = FALSE: all rows except the 8 ratio rows with numerator_NA_to_0 = FALSE must be identical; those differ by design (84 values, see ratio_legacy_differences.csv)
- **B_clustered**: weights + strata + cluster ids, adjust + domain lonely [ratio_legacy = TRUE]
- **C_quota_unweighted**: no weights, no strata (NGA quota design) [ratio_legacy = TRUE]
- **D_lonely_fail**: survey defaults: lonely.psu = fail (errors must match) [ratio_legacy = TRUE]
- **E_lonely_remove**: clustered, lonely.psu = remove + domain lonely [ratio_legacy = TRUE]
- **E_lonely_certainty**: clustered, lonely.psu = certainty + domain lonely [ratio_legacy = TRUE]
- **E_lonely_average**: clustered, lonely.psu = average + domain lonely [ratio_legacy = TRUE]
- **E_adjust_no_domain**: clustered, lonely.psu = adjust, adjust.domain.lonely = FALSE [ratio_legacy = TRUE]
- **F_fpc**: weights + strata + fpc (population size) [ratio_legacy = TRUE]
