test_that("create_loa_fast matches analysistools::create_loa", {
  skip_if_not_installed("analysistools")
  skip_if_not_installed("srvyr")
  data_env <- new.env()
  utils::data("analysistools_MSNA_template_data", package = "analysistools", envir = data_env)
  des <- srvyr::as_survey(data_env$analysistools_MSNA_template_data)
  for (gv in list(NULL, "admin1", c("admin1", "admin1, admin2"))) {
    old <- suppressMessages(analysistools::create_loa(des, group_var = gv, sm_separator = "/"))
    new <- create_loa_fast(des, group_var = gv, sm_separator = "/")
    expect_equal(new, as.data.frame(old), ignore_attr = TRUE)
  }
})

test_that("check_loa_fast keeps analysistools errors and fills ratio defaults", {
  d <- data.frame(a = 1:4, g = c("x", "x", "y", "y"), num = c(1, NA, 2, 3), den = c(2, 2, 0, 1))
  expect_error(check_loa_fast(data.frame(analysis = "mean"), d), "analysis_type")
  expect_error(check_loa_fast(data.frame(analysis_type = "quantile", analysis_var = "a", group_var = NA, level = .9), d),
               "not yet implemented")
  expect_error(check_loa_fast(data.frame(analysis_type = "mean", analysis_var = "zz", group_var = NA, level = .9), d),
               "analysis variables are not present")
  expect_error(check_loa_fast(data.frame(analysis_type = "mean", analysis_var = "a", group_var = "zz", level = .9), d),
               "group variables are not present")
  loa <- data.frame(analysis_type = "ratio", analysis_var = NA, group_var = NA, level = .9,
                    analysis_var_numerator = "num", analysis_var_denominator = "den",
                    numerator_NA_to_0 = NA, filter_denominator_0 = NA)
  warns <- character()
  out <- withCallingHandlers(check_loa_fast(loa, d), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(warns, 2)
  expect_true(all(grepl("set to TRUE", warns)))
  expect_true(out$numerator_NA_to_0)
  expect_true(out$filter_denominator_0)
})

test_that("ratio filters follow the documented rules; ratio_legacy reproduces the old bug", {
  d <- data.frame(num = c(1, NA, 0, 4, 2), den = c(3, 0, 2, NA, 0), g = "a")
  # documented: drop NA / 0 denominators, NA numerators -> 0 (default)
  r <- create_analysis_ratio_fast(d, analysis_var_numerator = "num", analysis_var_denominator = "den")
  expect_equal(r$stat, 1 / 5)
  expect_equal(r$n, 2)
  # numerator_NA_to_0 = FALSE: rows with NA numerator also removed, denominator filter kept
  r2 <- create_analysis_ratio_fast(d, analysis_var_numerator = "num", analysis_var_denominator = "den",
                                   numerator_NA_to_0 = FALSE)
  expect_equal(r2$stat, 1 / 5)
  expect_equal(r2$n, 2)
  # legacy: the denominator filter is dropped (zero denominators kept in n and in the ratio)
  r3 <- create_analysis_ratio_fast(d, analysis_var_numerator = "num", analysis_var_denominator = "den",
                                   numerator_NA_to_0 = FALSE, ratio_legacy = TRUE)
  expect_equal(r3$stat, (1 + 0 + 2) / (3 + 2 + 0))
  expect_equal(r3$n, 4)
})

test_that("add_weights_fast matches analysistools::add_weights", {
  skip_if_not_installed("analysistools")
  clean_data <- data.frame(uuid = 1:8, strata = c("s1", "s2", "s1", "s2", "s1", "s2", "s1", "s1"))
  sample <- data.frame(strata = c("s1", "s2"), population = c(30000, 50000))
  expect_equal(add_weights_fast(clean_data, sample, "strata", "strata", "population"),
               analysistools::add_weights(clean_data, sample, "strata", "strata", "population"),
               ignore_attr = TRUE)
  expect_error(add_weights_fast(clean_data, sample[1, ], "strata", "strata", "population"),
               "Not all strata from dataset")
})

test_that("run_disaggregated_analysis_fast writes one file per disaggregation and keeps ratio rows", {
  syn <- make_synthetic_msna(200, n_select_one = 2, n_select_multiple = 1, n_numeric = 1, n_filler = 0)
  out_dir <- file.path(tempdir(), "fastmsna_disagg_test")
  unlink(out_dir, recursive = TRUE)
  res <- run_disaggregated_analysis_fast(
    sh = "hh data", loa_sheet = syn$loa, disagg_vars = c("admin1", "Age & Gender"),
    disagg_index = 1, design = syn$data, weights = "weight", strata = "strata",
    lonely_psu = "adjust", adjust_domain_lonely = TRUE,
    output_dir = out_dir, save_format = "rds", verbose = FALSE
  )
  expect_true(file.exists(file.path(out_dir, "hh data_disagg_01_admin1.rds")))
  expect_true(any(res$admin1$analysis_type == "ratio"))
  expect_equal(fastmsna:::disagg_file_name(out_dir, "hh data", 2, "Age & Gender", "xlsx"),
               file.path(out_dir, "hh data_disagg_02_Age and Gender.xlsx"))
  old_style <- run_disaggregated_analysis_fast(
    sh = "hh data", loa_sheet = syn$loa, disagg_vars = "admin1", design = syn$data,
    weights = "weight", strata = "strata", output_dir = out_dir, keep_ratio_rows = FALSE,
    lonely_psu = "adjust", adjust_domain_lonely = TRUE, verbose = FALSE
  )
  expect_true(nrow(old_style$admin1) > 0)
  expect_false(any(old_style$admin1$analysis_type == "ratio"))
})

test_that("compare_analysis_results detects differences", {
  a <- data.frame(analysis_key = c("k1", "k2"), stat = c(1, 2), stat_low = 0, stat_upp = 3,
                  n = 1, n_total = 1, n_w = 1, n_w_total = 1)
  b <- a
  b$stat[2] <- 2.5
  cmp <- compare_analysis_results(a, b)
  expect_false(cmp$equivalent)
  expect_equal(nrow(cmp$mismatches), 1)
  expect_true(compare_analysis_results(a, a)$equivalent)
})
