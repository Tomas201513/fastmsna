#' fastmsna: fast, survey-correct descriptive analysis for MSNA datasets
#'
#' `fastmsna` re-implements the calculation layer of `analysistools`
#' (`create_analysis()` and the `create_analysis_*()` family) for large
#' datasets. Results (estimates, confidence intervals, `n`, `n_total`, `n_w`,
#' `n_w_total`, analysis keys) follow the same methodology as
#' `analysistools`/`srvyr`/`survey`:
#'
#' * weighted means, proportions and ratios with Taylor-linearised variances
#'   for stratified (and optionally clustered, with fpc) designs;
#' * domain (subgroup) estimation that keeps the full-design number of PSUs
#'   per stratum, exactly like `survey`;
#' * `survey.lonely.psu` / `survey.adjust.domain.lonely` handling;
#' * weighted medians with the "school" quantile rule and Woodruff
#'   confidence intervals.
#'
#' The speed-up comes from preparing the design once, computing grouping
#' indexes once per disaggregation, and computing every group of every
#' indicator in a handful of vectorised `rowsum()` passes, instead of building
#' one survey-design subset (with a full copy of the dataset) per group and per
#' indicator.
#'
#' Main entry points: [create_analysis_fast()], [run_disaggregated_analysis_fast()],
#' [as_fast_design()], [compare_analysis_results()].
#'
#' @keywords internal
#' @importFrom data.table data.table as.data.table setDT setorderv frankv rbindlist
#' @importFrom data.table := .N .I .SD .BY copy setnames setcolorder
#' @importFrom stats qt setNames
"_PACKAGE"

# Lets data.table syntax work inside the package namespace.
.datatable.aware <- TRUE

utils::globalVariables(c(
  ".row", ".ord", ".gs", ".level", ".g", ".l", ".na_row", ".key_group",
  "analysis_type", "analysis_var", "analysis_var_value", "group_var",
  "group_var_value", "stat", "stat_low", "stat_upp", "n", "n_total", "n_w",
  "n_w_total", "analysis_key", "se", "df", "g", "y", "w", "q", "lo", "up",
  ".c", ".unit", ".unit_label", ".occ", ".pos", "wt", "count", "population", "s",
  "analysis_var_numerator", "analysis_var_denominator", "numerator_NA_to_0",
  "filter_denominator_0", "level"
))
