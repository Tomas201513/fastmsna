# Single-indicator wrappers ------------------------------------------------
#
# Same signatures as analysistools::create_analysis_*(); each returns the
# results table for one indicator. They are convenient for ad-hoc checks; for
# many indicators use create_analysis_fast() with a LOA (one pass for all).

single_analysis <- function(design, row, ...) {
  des <- as_fast_design(design)
  loa <- check_loa_quiet(as.data.frame(row, stringsAsFactors = FALSE), des)
  run_loa(des, loa, on_error = "stop", ...)$results_table
}

check_loa_quiet <- function(loa, des) suppressWarnings(check_loa_fast(loa, des))

#' Fast single-indicator analyses
#'
#' Faster equivalents of `analysistools::create_analysis_mean()`,
#' `create_analysis_median()`, `create_analysis_prop_select_one()`,
#' `create_analysis_prop_select_multiple()` and `create_analysis_ratio()`,
#' with the same arguments and the same output table.
#'
#' @param design A srvyr/survey design, a `fastmsna_design` or (with
#'   [as_fast_design()] arguments) a data frame.
#' @param group_var Grouping variable(s) as one string, e.g. `"admin1"` or
#'   `"admin1, pop_group"`; `NA` (default) for no grouping.
#' @param analysis_var Analysis variable (select-multiple: the parent column).
#' @param level Confidence level, default 0.95.
#' @param sm_separator Select-multiple separator, default ".".
#' @param analysis_var_numerator,analysis_var_denominator Ratio numerator and
#'   denominator columns.
#' @param numerator_NA_to_0 Turn missing numerators into 0 (default `TRUE`).
#' @param filter_denominator_0 Remove rows with a denominator of 0 (default `TRUE`).
#' @param ratio_legacy See [create_analysis_fast()].
#' @param qrule Quantile rule for medians: `"school"` (analysistools default,
#'   average of the two middle values) or `"math"` (survey default, lower value).
#'
#' @return A results table (analysistools long format).
#' @name single_analysis_fast
#' @examples
#' df <- data.frame(x = c(1, 5, 3, NA, 8, 2), g = c("a", "a", "a", "b", "b", "b"))
#' create_analysis_mean_fast(df, analysis_var = "x")
#' create_analysis_median_fast(df, group_var = "g", analysis_var = "x")
NULL

#' @rdname single_analysis_fast
#' @export
create_analysis_mean_fast <- function(design, group_var = NA, analysis_var, level = .95) {
  single_analysis(design, list(analysis_type = "mean", analysis_var = analysis_var,
                               group_var = as.character(group_var), level = level))
}

#' @rdname single_analysis_fast
#' @export
create_analysis_median_fast <- function(design, group_var = NA, analysis_var, level = .95,
                                        qrule = c("school", "math")) {
  qrule <- match.arg(qrule)
  single_analysis(design, list(analysis_type = "median", analysis_var = analysis_var,
                               group_var = as.character(group_var), level = level),
                  qrule = qrule)
}

#' @rdname single_analysis_fast
#' @export
create_analysis_prop_select_one_fast <- function(design, group_var = NA, analysis_var, level = .95) {
  single_analysis(design, list(analysis_type = "prop_select_one", analysis_var = analysis_var,
                               group_var = as.character(group_var), level = level))
}

#' @rdname single_analysis_fast
#' @export
create_analysis_prop_select_multiple_fast <- function(design, group_var = NA, analysis_var,
                                                      level = .95, sm_separator = ".") {
  single_analysis(design, list(analysis_type = "prop_select_multiple", analysis_var = analysis_var,
                               group_var = as.character(group_var), level = level),
                  sm_separator = sm_separator)
}

#' @rdname single_analysis_fast
#' @export
create_analysis_ratio_fast <- function(design,
                                       group_var = NA,
                                       analysis_var_numerator,
                                       analysis_var_denominator,
                                       numerator_NA_to_0 = TRUE,
                                       filter_denominator_0 = TRUE,
                                       level = .95,
                                       ratio_legacy = FALSE) {
  single_analysis(design, list(analysis_type = "ratio", analysis_var = NA_character_,
                               group_var = as.character(group_var), level = level,
                               analysis_var_numerator = analysis_var_numerator,
                               analysis_var_denominator = analysis_var_denominator,
                               numerator_NA_to_0 = numerator_NA_to_0,
                               filter_denominator_0 = filter_denominator_0),
                  ratio_legacy = ratio_legacy)
}
