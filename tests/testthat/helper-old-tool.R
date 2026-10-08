# Run the analysistools functions row by row (like create_analysis_safe()),
# collecting results and errors, to compare with fastmsna.
old_analysis_safe <- function(design, loa, sm_separator = ".") {
  loa <- suppressWarnings(analysistools:::check_loa(loa, design))
  res <- list()
  err <- integer()
  for (i in seq_len(nrow(loa))) {
    row <- loa[i, ]
    r <- tryCatch(suppressWarnings(suppressMessages(switch(
      row$analysis_type,
      mean = analysistools::create_analysis_mean(design, group_var = row$group_var,
                                                 analysis_var = row$analysis_var, level = row$level),
      median = analysistools::create_analysis_median(design, group_var = row$group_var,
                                                     analysis_var = row$analysis_var, level = row$level),
      prop_select_one = analysistools::create_analysis_prop_select_one(
        design, group_var = row$group_var, analysis_var = row$analysis_var, level = row$level),
      prop_select_multiple = analysistools::create_analysis_prop_select_multiple(
        design, group_var = row$group_var, analysis_var = row$analysis_var, level = row$level,
        sm_separator = sm_separator),
      ratio = analysistools::create_analysis_ratio(
        design, group_var = row$group_var,
        analysis_var_numerator = row$analysis_var_numerator,
        analysis_var_denominator = row$analysis_var_denominator,
        numerator_NA_to_0 = row$numerator_NA_to_0,
        filter_denominator_0 = row$filter_denominator_0, level = row$level)
    ))), error = function(e) e)
    if (inherits(r, "error")) {
      err <- c(err, i)
    } else {
      r$loa_row <- rep(i, nrow(r))
      res[[length(res) + 1]] <- r
    }
  }
  list(results_table = do.call(rbind, res), errors = err)
}

expand_loa <- function(loa, group_vars, levels = c(0.95, 0.9)) {
  do.call(rbind, lapply(seq_along(group_vars), function(i) {
    l <- loa
    l$group_var <- group_vars[i]
    l$level <- levels[(i %% length(levels)) + 1]
    l
  }))
}

small_synthetic <- function(seed = 11, n = 300) {
  make_synthetic_msna(n, n_select_one = 3, n_select_multiple = 2, n_numeric = 2,
                      n_filler = 3, edge_cases = TRUE, seed = seed)
}

with_survey_options <- function(lonely, adj, code) {
  old <- options(survey.lonely.psu = lonely, survey.adjust.domain.lonely = adj)
  on.exit(options(old))
  force(code)
}
