#' Run the analysis one disaggregation at a time, saving after each one (fast)
#'
#' Faster drop-in for the `run_disaggregated_analysis()` helper of the NGA
#' MSNA scripts (`src/utils_data_analysis.R`). Same arguments, same output
#' files (`<file_prefix>_disagg_<NN>_<disaggregation>.<ext>`, with the same
#' Windows-safe file names), so `Step 3 - merge analysis.R` keeps working.
#'
#' Differences with the original helper:
#' * the design is prepared once ([as_fast_design()]) and converted columns
#'   are cached across disaggregations;
#' * analyses that cannot be computed are skipped and logged (CSV in
#'   `log_dir`) instead of failing the whole disaggregation (set
#'   `stop_on_error = TRUE` to stop instead);
#' * LOA rows without `analysis_var` (ratios) are kept. The original helper
#'   dropped them silently through `filter(analysis_var != group_var)`
#'   (`NA != x` is `NA`). Use `keep_ratio_rows = FALSE` for the old behaviour.
#'
#' @param sh Name of the sheet being analysed (messages / file names).
#' @param loa_sheet LOA for this sheet, *before* expansion over
#'   disaggregations (one row per indicator; `group_var` is overwritten).
#' @param disagg_vars Ordered vector of disaggregation variables.
#' @param disagg_index Which elements of `disagg_vars` to run, in order.
#' @param design Design for this sheet (srvyr design, `fastmsna_design`, or a
#'   data frame + `...` design arguments).
#' @param sm_separator Select-multiple separator.
#' @param output_dir Folder for the per-disaggregation results (created).
#' @param file_prefix File name prefix, defaults to `sh`.
#' @param save_format `"rds"`, `"xlsx"` (needs writexl) or `"csv"`.
#' @param stop_on_error Stop at the first failing analysis (`TRUE`) or skip
#'   and log it (`FALSE`, default).
#' @param keep_ratio_rows Keep LOA rows with a missing `analysis_var`
#'   (ratios). See Details.
#' @param log_dir Folder for the skipped-analysis logs (CSV). It is a
#'   sub-folder of `output_dir` by default so that Step 3 (which merges every
#'   `.xlsx` of the folder) does not pick the logs up.
#' @param verbose Print progress.
#' @param ... Passed to [as_fast_design()] when `design` is a data frame.
#'
#' @return Invisibly, a named list of results tables (one per disaggregation).
#' @export
run_disaggregated_analysis_fast <- function(sh,
                                            loa_sheet,
                                            disagg_vars,
                                            disagg_index = seq_along(disagg_vars),
                                            design,
                                            sm_separator = ".",
                                            output_dir = "./outputs/results",
                                            file_prefix = sh,
                                            save_format = c("rds", "xlsx", "csv"),
                                            stop_on_error = FALSE,
                                            keep_ratio_rows = TRUE,
                                            log_dir = file.path(output_dir, "_logs"),
                                            verbose = TRUE,
                                            ...) {
  save_format <- match.arg(save_format)
  if (save_format == "xlsx" && !requireNamespace("writexl", quietly = TRUE)) {
    stop("save_format = 'xlsx' needs the writexl package.", call. = FALSE)
  }
  if (any(disagg_index < 1 | disagg_index > length(disagg_vars))) {
    stop(sprintf("disagg_index must contain values between 1 and %d (length of disagg_vars). Got: %s",
                 length(disagg_vars), paste(disagg_index, collapse = ", ")), call. = FALSE)
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  t_prep <- proc.time()[["elapsed"]]
  des <- as_fast_design(design, ...)
  if (verbose) message(sprintf("[%s] design prepared in %.2f s (%d rows, %d strata, %d PSUs)",
                               sh, proc.time()[["elapsed"]] - t_prep, des$n, des$H, des$P))
  # check the LOA (and fill its defaults) once, not once per disaggregation
  loa_sheet <- as.data.frame(loa_sheet, stringsAsFactors = FALSE)
  loa_sheet$group_var <- NA_character_
  loa_sheet <- check_loa_fast(loa_sheet, des)

  n_total <- length(disagg_index)
  results_list <- list()
  for (i in seq_along(disagg_index)) {
    idx <- disagg_index[i]
    gvar <- disagg_vars[idx]
    t0 <- proc.time()[["elapsed"]]
    if (verbose) message(sprintf("[%s] (%d/%d) Running disaggregation index %d: '%s'", sh, i, n_total, idx, gvar))

    loa_sub <- as.data.frame(loa_sheet, stringsAsFactors = FALSE)
    loa_sub$group_var <- gvar
    av <- loa_sub$analysis_var
    drop <- if (keep_ratio_rows) !is.na(av) & av == gvar else is.na(av) | av == gvar
    loa_sub <- loa_sub[!drop, , drop = FALSE]

    run_result <- tryCatch(
      create_analysis_fast(des, loa = loa_sub, sm_separator = sm_separator,
                           on_error = if (stop_on_error) "stop" else "skip"),
      error = function(e) {
        message(sprintf("  ! Failed on '%s' (index %d): %s", gvar, idx, conditionMessage(e)))
        if (stop_on_error) stop(e)
        NULL
      }
    )
    if (is.null(run_result)) next

    results_table <- run_result$results_table
    out_file <- disagg_file_name(output_dir, file_prefix, idx, gvar, save_format)
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    save_status <- tryCatch({
      switch(save_format,
             rds = saveRDS(results_table, out_file),
             xlsx = writexl::write_xlsx(results_table, out_file),
             csv = utils::write.csv(results_table, out_file, row.names = FALSE, na = ""))
      TRUE
    }, error = function(e) {
      message(sprintf("  ! Failed to save '%s': %s", out_file, conditionMessage(e)))
      FALSE
    })
    if (!save_status) {
      if (stop_on_error) stop("Failed to save output file: ", out_file, call. = FALSE)
      next
    }
    if (nrow(run_result$skipped)) {
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
      log_file <- file.path(log_dir, sub("\\.[a-z]+$", "_skipped.csv", basename(out_file)))
      utils::write.csv(run_result$skipped, log_file, row.names = FALSE)
      if (verbose) message(sprintf("  ! %d analysis row(s) skipped, see %s", nrow(run_result$skipped), log_file))
    }
    if (verbose) message(sprintf("  -> saved: %s (%d rows, %.2f s)", out_file, nrow(results_table),
                                 proc.time()[["elapsed"]] - t0))
    results_list[[gvar]] <- results_table
  }
  invisible(results_list)
}

# Same file naming as the original helper (Windows-safe, MAX_PATH aware).
disagg_file_name <- function(output_dir, file_prefix, idx, gvar, save_format) {
  safe_gvar <- gsub("&", "and", gvar, fixed = TRUE)
  safe_gvar <- gsub('[<>:"/\\\\|?*]', "_", safe_gvar)
  full_len <- nchar(file.path(
    normalizePath(output_dir, winslash = "/", mustWork = FALSE),
    sprintf("%s_disagg_%02d_%s.%s", file_prefix, idx, safe_gvar, save_format)
  ))
  if (full_len > 259) {
    safe_gvar <- trimws(substr(safe_gvar, 1, nchar(safe_gvar) - (full_len - 259)))
    message(sprintf("  ! Path too long for Windows; filename shortened to '%s'", safe_gvar))
  }
  file.path(output_dir, sprintf("%s_disagg_%02d_%s.%s", file_prefix, idx, safe_gvar, save_format))
}
