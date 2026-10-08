#' Create an analysis from a design and a list of analysis (fast)
#'
#' Drop-in, faster replacement for `analysistools::create_analysis()`.
#' Same inputs (a srvyr design and a *list of analysis*, LOA), same output
#' (a list with `results_table`, `dataset`, `loa`) and the same results table
#' (columns, analysis keys, estimates, confidence intervals and sample sizes).
#' In addition the list contains `skipped`, the log of analyses that could not
#' be computed (as the `create_analysis_safe()` helper of the NGA scripts).
#'
#' The *loa* must contain `analysis_type`, `analysis_var`, `group_var` and
#' optionally `level` (default 0.95). Ratios also need
#' `analysis_var_numerator`, `analysis_var_denominator` and optionally
#' `numerator_NA_to_0` and `filter_denominator_0` (default `TRUE`).
#' `group_var` follows analysistools: `"admin1, admin2"` groups by the two
#' variables together; `NA` means no disaggregation.
#'
#' @param design A srvyr/survey design (`srvyr::as_survey_design()`), a
#'   `fastmsna_design` from [as_fast_design()] (fastest when the same data is
#'   analysed several times), or a data frame (then pass `weights`, `strata`,
#'   `ids`, `fpc` through `...`).
#' @param loa List of analysis. If `NULL`, one is created from the column types
#'   as in `analysistools::create_loa()`.
#' @param group_var Used only when `loa` is `NULL` (see
#'   `analysistools::create_loa()`).
#' @param sm_separator Separator between a select-multiple question and its
#'   choices. Default ".".
#' @param on_error `"stop"` (default, like `create_analysis()`) stops at the
#'   first analysis that cannot be computed; `"skip"` logs it in `$skipped`
#'   and continues (like `create_analysis_safe()`).
#' @param ratio_legacy `FALSE` (default) applies the documented ratio filters
#'   (rows with a missing or - if `filter_denominator_0` - zero denominator are
#'   always removed). `TRUE` reproduces an analysistools bug where these
#'   denominator filters are dropped when `numerator_NA_to_0 = FALSE`.
#' @param verbose Print progress / skipped messages.
#' @param ... Passed to [as_fast_design()] when `design` is a data frame
#'   (`weights`, `strata`, `ids`, `fpc`, `nest`) and to override lonely PSU
#'   options (`lonely_psu`, `adjust_domain_lonely`).
#'
#' @return A list with
#'   * `results_table`: long results table (analysistools format);
#'   * `dataset`: the dataset used;
#'   * `loa`: the list of analysis used;
#'   * `skipped`: data frame of analyses that failed (empty if none).
#' @export
#'
#' @examples
#' df <- data.frame(
#'   admin1 = rep(c("a", "b"), each = 6),
#'   income = c(10, 20, 30, NA, 50, 60, 15, 25, 35, 45, 55, 65),
#'   water = rep(c("tap", "well", "river"), 4),
#'   w = rep(c(1, 2), 6)
#' )
#' loa <- data.frame(
#'   analysis_type = c("mean", "median", "prop_select_one"),
#'   analysis_var = c("income", "income", "water"),
#'   group_var = c(NA, "admin1", "admin1"),
#'   level = 0.95
#' )
#' res <- create_analysis_fast(df, loa = loa, weights = "w", strata = "admin1")
#' res$results_table
create_analysis_fast <- function(design,
                                 loa = NULL,
                                 group_var = NULL,
                                 sm_separator = ".",
                                 on_error = c("stop", "skip"),
                                 ratio_legacy = FALSE,
                                 verbose = FALSE,
                                 ...) {
  on_error <- match.arg(on_error)
  des <- as_fast_design(design, ...)
  if (!is.null(loa) && !is.null(group_var)) {
    warning("You have provided a list of analysis and group variable, group variable will be ignored")
  }
  loa <- if (!is.null(loa)) {
    check_loa_fast(loa, des)
  } else {
    create_loa_fast(des, group_var = group_var, sm_separator = sm_separator)
  }
  run <- run_loa(des, loa, sm_separator = sm_separator, on_error = on_error,
                 ratio_legacy = ratio_legacy, verbose = verbose)
  list(
    results_table = run$results_table,
    dataset = des$variables,
    loa = loa,
    skipped = run$skipped
  )
}

# Internal: run a checked LOA on a fast design. Returns results with the
# internal `.row` column (LOA row of every result row) when keep_row = TRUE.
run_loa <- function(des, loa, sm_separator = ".", on_error = "stop", ratio_legacy = FALSE,
                    verbose = FALSE, keep_row = FALSE, qrule = "school") {
  L <- as.data.table(loa)
  L[, .row := .I]
  L[, .gs := normalise_group_var(group_var)]
  L[, .level := suppressWarnings(as.numeric(level))]

  ctx <- new.env(parent = emptyenv())
  ctx$des <- des
  ctx$gs <- new.env(parent = emptyenv())
  ctx$lonely <- lonely_options(des)
  ctx$on_error <- on_error
  ctx$ratio_legacy <- ratio_legacy
  ctx$sm_separator <- sm_separator
  ctx$verbose <- verbose
  ctx$qrule <- qrule
  ctx$skipped <- list()
  ctx$max_cells <- getOption("fastmsna.max_cells", 4e6)

  runners <- list(
    prop_select_one = run_select_one,
    prop_select_multiple = run_select_multiple,
    mean = run_mean,
    median = run_median,
    ratio = run_ratio
  )
  out <- list()
  for (tp in names(runners)) {
    rows <- L[analysis_type == tp]
    if (!nrow(rows)) next
    t0 <- proc.time()[["elapsed"]]
    res <- runners[[tp]](rows, ctx)
    if (verbose) {
      message(sprintf("  %-21s %5d loa rows  %8.2f s", tp, nrow(rows), proc.time()[["elapsed"]] - t0))
    }
    out[[tp]] <- res[, intersect(result_columns, names(res)), with = FALSE]
  }
  res <- if (length(out)) rbindlist(out, use.names = TRUE, fill = TRUE) else empty_result()
  setorderv(res, c(".row", ".ord"))
  res[, analysis_type := L$analysis_type[.row]]

  is_ratio <- res$analysis_type == "ratio"
  key <- character(nrow(res))
  if (any(!is_ratio)) {
    i <- which(!is_ratio)
    key[i] <- paste0(res$analysis_type[i], " @/@ ", res$analysis_var[i], " %/% ",
                     res$analysis_var_value[i], " @/@ ", res$.key_group[i])
  }
  if (any(is_ratio)) {
    i <- which(is_ratio)
    av <- str_split_keep(res$analysis_var[i], " %/% ")
    vv <- str_split_keep(res$analysis_var_value[i], " %/% ")
    part <- vapply(seq_along(i), function(j) {
      paste(paste(av[[j]], vv[[j]], sep = " %/% "), collapse = " -/- ")
    }, character(1))
    key[i] <- paste(res$analysis_type[i], "@/@", part, "@/@", res$.key_group[i])
  }
  res[, analysis_key := key]

  final_cols <- c("analysis_type", "analysis_var", "analysis_var_value", "group_var",
                  "group_var_value", "stat", "stat_low", "stat_upp", "n", "n_total",
                  "n_w", "n_w_total", "analysis_key")
  if (keep_row) final_cols <- c(final_cols, ".row")
  res <- res[, final_cols, with = FALSE]
  for (cc in c("n", "n_total", "n_w", "n_w_total", "stat", "stat_low", "stat_upp")) {
    set(res, j = cc, value = as.numeric(res[[cc]]))
  }

  skipped <- if (length(ctx$skipped)) rbindlist(ctx$skipped) else
    data.table(row_index = integer(), analysis_type = character(), analysis_var = character(),
               group_var = character(), n_obs = integer(), error_message = character())
  setorderv(skipped, "row_index")
  list(results_table = as_output(res), skipped = as.data.frame(skipped))
}

# Results are returned as tibbles when tibble is installed (as analysistools).
as_output <- function(dt) {
  df <- as.data.frame(dt)
  if (requireNamespace("tibble", quietly = TRUE)) df <- tibble::as_tibble(df)
  df
}

#' @importFrom data.table set
NULL

#' Check a list of analysis (fast design version)
#'
#' Same checks as `analysistools::check_loa()`: required columns, implemented
#' analysis types, existing variables, and defaults for `level`,
#' `numerator_NA_to_0` and `filter_denominator_0`. In addition, missing
#' (`NA`) `numerator_NA_to_0` / `filter_denominator_0` values on ratio rows
#' are set to the default `TRUE` with a warning (analysistools would stop with
#' "missing value where TRUE/FALSE needed").
#'
#' @param loa A list of analysis.
#' @param design A design accepted by [as_fast_design()].
#' @return The checked list of analysis (a data frame).
#' @export
check_loa_fast <- function(loa, design) {
  des <- as_fast_design(design)
  loa <- as.data.frame(loa, stringsAsFactors = FALSE)
  analysis_type_dictionary <- c("prop_select_one", "prop_select_multiple", "mean", "median", "ratio")

  if (!"analysis_type" %in% names(loa)) stop("Make sure to have analysis_type in your loa")
  if (any(analysis_type_dictionary[analysis_type_dictionary != "ratio"] %in% loa[["analysis_type"]])) {
    if (!all(c("group_var", "analysis_var") %in% names(loa))) {
      stop("Make sure to have group_var, analysis_var in your loa")
    }
  }
  if (!"group_var" %in% names(loa)) loa[["group_var"]] <- NA_character_
  if (!"analysis_var" %in% names(loa)) loa[["analysis_var"]] <- NA_character_
  if (!"level" %in% names(loa)) {
    loa[["level"]] <- 0.95
    warning("No column level identified, set to 0.95 as default value.")
  }
  if ("ratio" %in% loa[["analysis_type"]]) {
    if (!all(c("analysis_var_numerator", "analysis_var_denominator") %in% names(loa))) {
      stop("You have ratio, you need analysis_var_numerator, analysis_var_denominator columns in your loa")
    }
    is_ratio <- loa[["analysis_type"]] %in% "ratio"
    for (flag in c("numerator_NA_to_0", "filter_denominator_0")) {
      if (!flag %in% names(loa)) {
        loa[[flag]] <- ifelse(is_ratio, TRUE, NA)
        warning("No column level ", flag, ", set to TRUE as default value.")
      } else {
        val <- as.logical(loa[[flag]])
        if (any(is_ratio & is.na(val))) {
          warning(flag, " is missing for ", sum(is_ratio & is.na(val)),
                  " ratio row(s), set to TRUE as default value.")
          val[is_ratio & is.na(val)] <- TRUE
        }
        loa[[flag]] <- val
      }
    }
  }
  verify_in(loa[["analysis_type"]], analysis_type_dictionary,
            "The following analysis type are not yet implemented or check for typo: ")
  gv <- normalise_group_var(loa[["group_var"]])
  gv_vars <- unique(unlist(lapply(gv[!is.na(gv)], split_group_var)))
  verify_in(gv_vars, names(des$variables), "The following group variables are not present in the dataset: ")
  verify_in(loa[["analysis_var"]][!loa[["analysis_type"]] %in% "ratio"], names(des$variables),
            "The following analysis variables are not present in the dataset: ")
  if ("ratio" %in% loa[["analysis_type"]]) {
    is_ratio <- loa[["analysis_type"]] %in% "ratio"
    verify_in(loa[["analysis_var_numerator"]][is_ratio], names(des$variables),
              "The following analysis numerator variables are not present in the dataset: ")
    verify_in(loa[["analysis_var_denominator"]][is_ratio], names(des$variables),
              "The following analysis denominator variables are not present in the dataset: ")
  }
  loa
}

verify_in <- function(a, b, msg) {
  ua <- unique(a)
  ua <- ua[!is.na(ua)]
  miss <- ua[!ua %in% b]
  if (length(miss)) stop(msg, paste(miss, collapse = ", "), call. = FALSE)
}

#' Create a list of analysis from a design (fast design version)
#'
#' Same rules as `analysistools::create_loa()`: character columns ->
#' `prop_select_one`, double/integer -> `mean` and `median`, logical and
#' select-multiple parents -> `prop_select_multiple`; select-multiple choice
#' columns, columns starting with `_` / `X_` and `start`, `end`, `today`,
#' `uuid` are ignored. Each `group_var` is added after the overall analysis.
#'
#' @param design A design accepted by [as_fast_design()].
#' @param group_var Grouping variable(s), see `analysistools::create_loa()`.
#' @param sm_separator Select-multiple separator.
#' @return A list of analysis (data frame).
#' @export
create_loa_fast <- function(design, group_var = NULL, sm_separator = ".") {
  des <- as_fast_design(design)
  vars <- names(des$variables)
  types <- vapply(des$variables, typeof, character(1))
  dict <- data.frame(
    type = c("character", "double", "double", "logical", "integer", "integer"),
    analysis_type = c("prop_select_one", "mean", "median", "prop_select_multiple", "mean", "median"),
    stringsAsFactors = FALSE
  )
  # cleaningtools::auto_detect_sm_parents() / auto_sm_parent_children()
  stripped <- sub(paste0(".[^\\", sm_separator, "]*$"), "", vars)
  stripped <- stripped[stripped != ""]
  tab <- table(stripped)
  parents <- names(tab)[tab > 1]
  children <- if (length(parents)) {
    pre <- tolower(paste0(parents, sm_separator))
    vars[vapply(tolower(vars), function(v) any(startsWith(v, pre)), logical(1))]
  } else character(0)

  rows <- lapply(seq_along(vars), function(i) {
    at <- dict$analysis_type[dict$type == types[i]]
    if (!length(at)) return(NULL)
    data.frame(analysis_var = vars[i], analysis_type = at, stringsAsFactors = FALSE)
  })
  loa <- do.call(rbind, rows)
  if (is.null(loa)) loa <- data.frame(analysis_var = character(), analysis_type = character())
  loa <- loa[!grepl("^(X_|_)", loa$analysis_var), , drop = FALSE]
  loa <- loa[!loa$analysis_var %in% c("start", "end", "today", "uuid"), , drop = FALSE]
  loa <- loa[!loa$analysis_var %in% children, , drop = FALSE]
  loa$analysis_type[loa$analysis_var %in% parents] <- "prop_select_multiple"

  if (is.null(group_var)) {
    loa$group_var <- NA_character_
  } else {
    loa <- do.call(rbind, lapply(c(NA, group_var), function(x) {
      loa$group_var <- x
      loa
    }))
    keep <- is.na(loa$group_var) |
      !mapply(function(p, s) grepl(p, s, perl = TRUE), loa$analysis_var, loa$group_var)
    loa <- loa[keep, , drop = FALSE]
  }
  rownames(loa) <- NULL
  out <- loa[, c("analysis_type", "analysis_var", "group_var")]
  out$level <- 0.95
  out
}
