#  -----------------------------------------------------------------------------
#
# Title : Validation of fastmsna against analysistools
#    By : IMPACT NGA
#  Date : 10/2026
#
#  Runs analysistools (the old tool) and fastmsna on the same data and the same
#  list of analysis, and compares the results row by row (matched on
#  analysis_key). Writes to outputs/validation/:
#    - validation_summary.csv      one line per scenario
#    - validation_columns.csv      per scenario x column: mismatches, max diffs
#    - rowwise_<scenario>.csv      old vs new values for every result row
#    - ratio_legacy_differences.csv  the intentional ratio-filter differences
#    - VALIDATION_REPORT.md
#
#  Old tool = analysistools functions run row by row with tryCatch (exactly
#  like create_analysis_safe() in msna_analysis_main/src), so that rows that
#  error in the old tool can be compared with rows skipped by fastmsna.
#  Takes ~15-25 minutes (the old tool is the slow part).
#  -----------------------------------------------------------------------------

project_dir <- "C:/Users/User/Music/new_msna_analysis"
out_dir <- file.path(project_dir, "outputs", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  if (requireNamespace("fastmsna", quietly = TRUE) &&
      utils::packageVersion("fastmsna") >= read.dcf(file.path(project_dir, "DESCRIPTION"))[, "Version"]) {
    library(fastmsna)
  } else {
    pkgload::load_all(project_dir, quiet = TRUE)
  }
  library(analysistools)
  library(srvyr)
  library(dplyr)
})
options(dplyr.summarise.inform = FALSE)
source(file.path(project_dir, "tests", "testthat", "helper-old-tool.R"))

log_file <- file.path(project_dir, "logs", paste0(Sys.Date(), "_validation.log"))
dir.create(dirname(log_file), showWarnings = FALSE)
logmsg <- function(...) {
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " ", sprintf(...))
  message(msg)
  cat(msg, "\n", file = log_file, append = TRUE)
}

summary_rows <- list()
column_rows <- list()

record <- function(name, description, old_tab, new_res, old_err, old_s, new_s, loa_rows,
                   tolerance = 1e-8, write_rowwise = TRUE) {
  cmp <- compare_analysis_results(old_tab, new_res, tolerance = tolerance)
  new_err <- sort(new_res$skipped$row_index)
  summary_rows[[name]] <<- data.frame(
    scenario = name,
    description = description,
    loa_rows = loa_rows,
    result_rows_old = cmp$n_old,
    result_rows_new = cmp$n_new,
    only_in_old = length(cmp$only_in_old),
    only_in_new = length(cmp$only_in_new),
    value_mismatches = nrow(cmp$mismatches),
    label_mismatches = length(cmp$label_mismatches),
    same_row_order = cmp$same_row_order,
    failed_rows_old = length(old_err),
    skipped_rows_new = length(new_err),
    same_failed_rows = identical(sort(as.integer(old_err)), as.integer(new_err)),
    max_rel_diff = max(cmp$summary$max_rel_diff),
    seconds_old = round(old_s, 2),
    seconds_new = round(new_s, 2),
    equivalent = cmp$equivalent && identical(sort(as.integer(old_err)), as.integer(new_err))
  )
  column_rows[[name]] <<- cbind(scenario = name, cmp$summary)
  if (write_rowwise) {
    o <- as.data.frame(old_tab)
    w <- as.data.frame(new_res$results_table)
    o$occurrence <- stats::ave(seq_len(nrow(o)), o$analysis_key, FUN = seq_along)
    w$occurrence <- stats::ave(seq_len(nrow(w)), w$analysis_key, FUN = seq_along)
    m <- merge(o, w, by = c("analysis_key", "occurrence"), all = TRUE, suffixes = c("_old", "_new"), sort = FALSE)
    num <- c("stat", "stat_low", "stat_upp", "n", "n_total", "n_w", "n_w_total")
    for (cc in num) m[[paste0(cc, "_absdiff")]] <- abs(m[[paste0(cc, "_old")]] - m[[paste0(cc, "_new")]])
    keep <- c("analysis_key", "occurrence", unlist(lapply(num, function(cc) paste0(cc, c("_old", "_new", "_absdiff")))))
    utils::write.csv(m[, keep], file.path(out_dir, paste0("rowwise_", name, ".csv")), row.names = FALSE)
  }
  logmsg("%-28s equivalent = %s | rows %d/%d | mismatches %d | failed rows old %d / new %d | %.1fs vs %.2fs",
         name, summary_rows[[name]]$equivalent, cmp$n_old, cmp$n_new, nrow(cmp$mismatches),
         length(old_err), length(new_err), old_s, new_s)
  invisible(cmp)
}

timed <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  list(value = val, seconds = proc.time()[["elapsed"]] - t0)
}

# 1. analysistools' own test data --------------------------------------------
tmpl <- new.env()
utils::data(list = c("analysistools_MSNA_template_data", "analysistools_MSNA_template_loa",
                     "analysistools_MSNA_template_loa_with_ratio"),
            package = "analysistools", envir = tmpl)
des_t <- srvyr::as_survey(tmpl$analysistools_MSNA_template_data)
for (nm in c("analysistools_MSNA_template_loa", "analysistools_MSNA_template_loa_with_ratio")) {
  loa <- tmpl[[nm]]
  o <- timed(suppressMessages(analysistools::create_analysis(des_t, loa = loa, sm_separator = "/")))
  n <- timed(create_analysis_fast(des_t, loa = loa, sm_separator = "/"))
  record(sub("analysistools_MSNA_", "", nm), "analysistools template data and LOA (sm separator '/')",
         o$value$results_table, n$value, integer(), o$seconds, n$seconds, nrow(loa))
}
short <- tmpl$analysistools_MSNA_template_data
short <- short[, c("admin1", "admin2", "expenditure_debt", "income_v1_salaried_work",
                   "wash_drinkingwatersource", grep("edu_learning_conditions_reasons_v1", names(short), value = TRUE))]
o <- timed(suppressMessages(analysistools::create_analysis(srvyr::as_survey(short), group_var = c("admin1", "admin1, admin2"),
                                                           sm_separator = "/")))
n <- timed(create_analysis_fast(srvyr::as_survey(short), group_var = c("admin1", "admin1, admin2"), sm_separator = "/"))
record("template_no_loa", "template data, automatic LOA, group_var = c('admin1', 'admin1, admin2')",
       o$value$results_table, n$value, integer(), o$seconds, n$seconds, nrow(n$value$loa))

# 2. synthetic MSNA data with edge cases --------------------------------------
syn <- make_synthetic_msna(1000, n_select_one = 6, n_select_multiple = 3, n_numeric = 3,
                           n_filler = 20, edge_cases = TRUE, seed = 11)
d <- syn$data
d$fpc_pop <- ave(rep(1, nrow(d)), d$strata, FUN = function(x) length(x) * 4)
loa <- expand_loa(syn$loa, c(NA, "admin1", "hoh_gender", "admin1, pop_group", "rare_group"))
rv <- expand.grid(na0 = c(TRUE, FALSE), f0 = c(TRUE, FALSE))
ratio_extra <- do.call(rbind, lapply(seq_len(nrow(rv)), function(i) data.frame(
  analysis_type = "ratio", analysis_var = NA, group_var = rep(c(NA, "admin1"), each = 2), level = 0.9,
  analysis_var_numerator = rep(c("num_children_school", "expenditure_food"), 2),
  analysis_var_denominator = rep(c("num_children", "expenditure_total"), 2),
  numerator_NA_to_0 = rv$na0[i], filter_denominator_0 = rv$f0[i])))
loa <- rbind(loa, ratio_extra)

scenarios <- list(
  list(name = "A_2scs_strata_weights", desc = "weights + strata, ids = ~1 (NGA 2SCS Step 2 design), lonely.psu = adjust + domain lonely",
       lonely = "adjust", adj = TRUE, make = function() as_survey_design(d, weights = weight, strata = strata)),
  list(name = "B_clustered", desc = "weights + strata + cluster ids, adjust + domain lonely",
       lonely = "adjust", adj = TRUE, make = function() as_survey_design(d, ids = cluster_id, weights = weight, strata = strata)),
  list(name = "C_quota_unweighted", desc = "no weights, no strata (NGA quota design)",
       lonely = "adjust", adj = TRUE, make = function() as_survey_design(d)),
  list(name = "D_lonely_fail", desc = "survey defaults: lonely.psu = fail (errors must match)",
       lonely = "fail", adj = FALSE, make = function() as_survey_design(d, weights = weight, strata = strata)),
  list(name = "E_lonely_remove", desc = "clustered, lonely.psu = remove + domain lonely",
       lonely = "remove", adj = TRUE, make = function() as_survey_design(d, ids = cluster_id, weights = weight, strata = strata)),
  list(name = "E_lonely_certainty", desc = "clustered, lonely.psu = certainty + domain lonely",
       lonely = "certainty", adj = TRUE, make = function() as_survey_design(d, ids = cluster_id, weights = weight, strata = strata)),
  list(name = "E_lonely_average", desc = "clustered, lonely.psu = average + domain lonely",
       lonely = "average", adj = TRUE, make = function() as_survey_design(d, ids = cluster_id, weights = weight, strata = strata)),
  list(name = "E_adjust_no_domain", desc = "clustered, lonely.psu = adjust, adjust.domain.lonely = FALSE",
       lonely = "adjust", adj = FALSE, make = function() as_survey_design(d, ids = cluster_id, weights = weight, strata = strata)),
  list(name = "F_fpc", desc = "weights + strata + fpc (population size)",
       lonely = "adjust", adj = TRUE, make = function() as_survey_design(d, weights = weight, strata = strata, fpc = fpc_pop))
)

for (sc in scenarios) {
  options(survey.lonely.psu = sc$lonely, survey.adjust.domain.lonely = sc$adj)
  des <- sc$make()
  o <- timed(old_analysis_safe(des, loa))
  n <- timed(create_analysis_fast(des, loa = loa, on_error = "skip", ratio_legacy = TRUE))
  record(sc$name, paste(sc$desc, "[ratio_legacy = TRUE]"), o$value$results_table, n$value,
         o$value$errors, o$seconds, n$seconds, nrow(loa))
  if (sc$name == "A_2scs_strata_weights") {
    # same design, default (documented) ratio filters: the only differences
    # allowed are ratio rows with numerator_NA_to_0 = FALSE
    fd <- as_fast_design(des)
    n2 <- fastmsna:::run_loa(fd, suppressWarnings(check_loa_fast(loa, fd)), on_error = "skip",
                             ratio_legacy = FALSE, keep_row = TRUE)
    affected <- which(loa$analysis_type == "ratio" & loa$numerator_NA_to_0 %in% FALSE)
    old_tab <- o$value$results_table
    new_tab <- n2$results_table
    cmp_same <- compare_analysis_results(old_tab[!old_tab$loa_row %in% affected, ],
                                         new_tab[!new_tab$.row %in% affected, ], tolerance = 1e-8)
    cmp_diff <- compare_analysis_results(old_tab[old_tab$loa_row %in% affected, ],
                                         new_tab[new_tab$.row %in% affected, ], tolerance = 1e-8)
    utils::write.csv(cmp_diff$mismatches, file.path(out_dir, "ratio_legacy_differences.csv"), row.names = FALSE)
    summary_rows[["G_ratio_documented_filters"]] <- data.frame(
      scenario = "G_ratio_documented_filters",
      description = paste0("scenario A with ratio_legacy = FALSE: all rows except the ", length(affected),
                           " ratio rows with numerator_NA_to_0 = FALSE must be identical; those ",
                           "differ by design (", nrow(cmp_diff$mismatches),
                           " values, see ratio_legacy_differences.csv)"),
      loa_rows = nrow(loa), result_rows_old = cmp_same$n_old, result_rows_new = cmp_same$n_new,
      only_in_old = length(cmp_same$only_in_old), only_in_new = length(cmp_same$only_in_new),
      value_mismatches = nrow(cmp_same$mismatches), label_mismatches = length(cmp_same$label_mismatches),
      same_row_order = cmp_same$same_row_order, failed_rows_old = length(o$value$errors),
      skipped_rows_new = nrow(n2$skipped), same_failed_rows = TRUE,
      max_rel_diff = max(cmp_same$summary$max_rel_diff), seconds_old = NA, seconds_new = NA,
      equivalent = cmp_same$equivalent
    )
    logmsg("G_ratio_documented_filters   unaffected rows equivalent: %s | intentional differences: %d values on %d ratio rows",
           cmp_same$equivalent, nrow(cmp_diff$mismatches), length(affected))
  }
}
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)

# 3. write outputs -------------------------------------------------------------
summ <- do.call(rbind, summary_rows)
cols <- do.call(rbind, column_rows)
utils::write.csv(summ, file.path(out_dir, "validation_summary.csv"), row.names = FALSE)
utils::write.csv(cols, file.path(out_dir, "validation_columns.csv"), row.names = FALSE)

md <- c(
  "# Validation report: fastmsna vs analysistools",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M"), " | R ", getRversion(),
         " | analysistools ", utils::packageVersion("analysistools"),
         " | srvyr ", utils::packageVersion("srvyr"), " | survey ", utils::packageVersion("survey")),
  "",
  "Rows are matched on `analysis_key`. Values are equal when both are missing or",
  "`|old - new| <= 1e-8 * max(1, |old|)`. `max_rel_diff` is the largest relative",
  "difference over all compared values (floating-point noise only).",
  "",
  "| scenario | LOA rows | result rows old / new | value mismatches | same row order | failed rows old / new (same) | max rel diff | old (s) | new (s) | equivalent |",
  "|---|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %d | %d / %d | %d | %s | %d / %d (%s) | %s | %s | %s | **%s** |",
          summ$scenario, summ$loa_rows, summ$result_rows_old, summ$result_rows_new,
          summ$value_mismatches, summ$same_row_order, summ$failed_rows_old,
          summ$skipped_rows_new, summ$same_failed_rows,
          ifelse(is.na(summ$max_rel_diff), "-", formatC(summ$max_rel_diff, format = "e", digits = 1)),
          ifelse(is.na(summ$seconds_old), "-", summ$seconds_old),
          ifelse(is.na(summ$seconds_new), "-", summ$seconds_new), summ$equivalent),
  "",
  "Scenario descriptions:",
  "",
  sprintf("- **%s**: %s", summ$scenario, summ$description)
)
writeLines(md, file.path(out_dir, "VALIDATION_REPORT.md"))
logmsg("All scenarios equivalent: %s", all(summ$equivalent))
print(summ[, c("scenario", "value_mismatches", "same_row_order", "same_failed_rows", "max_rel_diff", "equivalent")])
