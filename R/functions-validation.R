#' Compare two results tables row by row
#'
#' Matches two analysistools-format results tables on `analysis_key` (and the
#' occurrence number when a key is repeated) and compares the numeric columns.
#' Use it to validate a migration: run the old and the new tool on the same
#' data and LOA, then `compare_analysis_results(old, new)`.
#'
#' Two values are considered equal when both are missing (`NA`/`NaN`), or when
#' `|old - new| <= tolerance * max(1, |old|)`.
#'
#' @param old,new Results tables, or the lists returned by `create_analysis()`
#'   / [create_analysis_fast()] (their `results_table` is used).
#' @param tolerance Relative tolerance (default `1e-8`).
#' @param columns Numeric columns to compare.
#' @param key Name of the key column.
#'
#' @return A list of class `fastmsna_comparison` with
#'   * `summary`: one row per compared column (mismatches, max abs/rel diff);
#'   * `mismatches`: the rows/columns that differ (long format);
#'   * `only_in_old`, `only_in_new`: keys present in one table only;
#'   * `label_mismatches`: rows whose label columns differ;
#'   * `same_row_order`: whether both tables list the keys in the same order;
#'   * `equivalent`: `TRUE` if nothing differs.
#' @export
compare_analysis_results <- function(old, new,
                                     tolerance = 1e-8,
                                     columns = c("stat", "stat_low", "stat_upp", "n",
                                                 "n_total", "n_w", "n_w_total"),
                                     key = "analysis_key") {
  if (!is.data.frame(old) && is.list(old)) old <- old$results_table
  if (!is.data.frame(new) && is.list(new)) new <- new$results_table
  o <- as.data.table(old)
  w <- as.data.table(new)
  o[, .occ := seq_len(.N), by = key]
  w[, .occ := seq_len(.N), by = key]
  o[, .pos := .I]
  w[, .pos := .I]
  m <- merge(o, w, by = c(key, ".occ"), all = TRUE, suffixes = c(".old", ".new"), sort = FALSE)
  in_old <- !is.na(m$.pos.old)
  in_new <- !is.na(m$.pos.new)
  both <- m[in_old & in_new]

  summ <- list()
  mism <- list()
  for (cc in columns) {
    a <- as.numeric(both[[paste0(cc, ".old")]])
    b <- as.numeric(both[[paste0(cc, ".new")]])
    na_a <- is.na(a)
    na_b <- is.na(b)
    d <- abs(a - b)
    tol <- tolerance * pmax(1, abs(a))
    bad <- (na_a != na_b) | (!na_a & !na_b & !(d <= tol))
    bad[is.na(bad)] <- TRUE
    ok <- !na_a & !na_b
    summ[[cc]] <- data.frame(
      column = cc,
      n_compared = length(a),
      n_both_missing = sum(na_a & na_b),
      n_mismatch = sum(bad),
      max_abs_diff = if (any(ok)) max(d[ok]) else 0,
      max_rel_diff = if (any(ok)) max(d[ok] / pmax(1e-300, abs(a[ok]))) else 0
    )
    if (any(bad)) {
      mism[[cc]] <- data.frame(analysis_key = both[[key]][bad], column = cc,
                               old = a[bad], new = b[bad])
    }
  }
  label_cols <- intersect(c("analysis_type", "analysis_var", "analysis_var_value",
                            "group_var", "group_var_value"),
                          intersect(names(o), names(w)))
  lab_bad <- rep(FALSE, nrow(both))
  for (cc in label_cols) {
    x <- as.character(both[[paste0(cc, ".old")]])
    y <- as.character(both[[paste0(cc, ".new")]])
    lab_bad <- lab_bad | !((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y))
  }
  same_order <- nrow(o) == nrow(w) && identical(as.character(o[[key]]), as.character(w[[key]]))
  out <- list(
    summary = do.call(rbind, summ),
    mismatches = if (length(mism)) do.call(rbind, mism) else
      data.frame(analysis_key = character(), column = character(), old = numeric(), new = numeric()),
    only_in_old = m[[key]][in_old & !in_new],
    only_in_new = m[[key]][!in_old & in_new],
    label_mismatches = both[[key]][lab_bad],
    same_row_order = same_order,
    n_old = nrow(o),
    n_new = nrow(w)
  )
  out$equivalent <- nrow(out$mismatches) == 0 && !length(out$only_in_old) &&
    !length(out$only_in_new) && !length(out$label_mismatches)
  class(out) <- "fastmsna_comparison"
  out
}

#' @export
print.fastmsna_comparison <- function(x, ...) {
  cat("<fastmsna_comparison> old rows: ", x$n_old, " | new rows: ", x$n_new,
      " | same row order: ", x$same_row_order, "\n", sep = "")
  cat("  only in old: ", length(x$only_in_old), " | only in new: ", length(x$only_in_new),
      " | label mismatches: ", length(x$label_mismatches), "\n", sep = "")
  print(x$summary, row.names = FALSE)
  cat(if (x$equivalent) "=> EQUIVALENT\n" else "=> DIFFERENCES FOUND (see $mismatches)\n")
  invisible(x)
}
