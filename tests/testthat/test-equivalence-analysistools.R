# Equivalence with analysistools on its own test data (template data / LOA /
# expected results shipped with analysistools).

test_that("create_analysis_fast reproduces analysistools on the template data", {
  skip_if_not_installed("analysistools")
  skip_if_not_installed("srvyr")
  data_env <- new.env()
  utils::data(list = c("analysistools_MSNA_template_data", "analysistools_MSNA_template_loa",
                       "analysistools_MSNA_template_loa_with_ratio"),
              package = "analysistools", envir = data_env)
  des <- srvyr::as_survey(data_env$analysistools_MSNA_template_data)

  for (loa in list(data_env$analysistools_MSNA_template_loa,
                   data_env$analysistools_MSNA_template_loa_with_ratio)) {
    old <- suppressMessages(analysistools::create_analysis(des, loa = loa, sm_separator = "/"))
    new <- create_analysis_fast(des, loa = loa, sm_separator = "/")
    cmp <- compare_analysis_results(old, new, tolerance = 1e-10)
    expect_true(cmp$equivalent)
    expect_true(cmp$same_row_order)
    expect_identical(names(new$results_table), names(old$results_table))
  }
})

test_that("create_analysis_fast without a LOA reproduces create_analysis(group_var =)", {
  skip_if_not_installed("analysistools")
  skip_if_not_installed("srvyr")
  data_env <- new.env()
  utils::data("analysistools_MSNA_template_data", package = "analysistools", envir = data_env)
  d <- data_env$analysistools_MSNA_template_data
  short <- d[, c("admin1", "admin2", "expenditure_debt", "income_v1_salaried_work",
                 "wash_drinkingwatersource",
                 grep("edu_learning_conditions_reasons_v1", names(d), value = TRUE))]
  des <- srvyr::as_survey(short)
  old <- suppressMessages(analysistools::create_analysis(des, group_var = c("admin1", "admin1, admin2"),
                                                         sm_separator = "/"))
  new <- create_analysis_fast(des, group_var = c("admin1", "admin1, admin2"), sm_separator = "/")
  expect_equal(new$loa, as.data.frame(old$loa), ignore_attr = TRUE)
  cmp <- compare_analysis_results(old, new, tolerance = 1e-10)
  expect_true(cmp$equivalent)
})

test_that("weighted, stratified design with edge cases matches analysistools row by row", {
  skip_if_not_installed("analysistools")
  skip_if_not_installed("srvyr")
  syn <- small_synthetic()
  d <- syn$data
  loa <- expand_loa(syn$loa, c(NA, "admin1", "hoh_gender", "admin1, pop_group", "rare_group"))
  with_survey_options("adjust", TRUE, {
    des <- srvyr::as_survey_design(d, weights = weight, strata = strata)
    old <- old_analysis_safe(des, loa)
    new <- create_analysis_fast(des, loa = loa, on_error = "skip")
  })
  cmp <- compare_analysis_results(old, new, tolerance = 1e-8)
  expect_true(cmp$equivalent)
  expect_identical(sort(new$skipped$row_index), sort(old$errors))
})

test_that("clustered design and lonely.psu = 'fail' give the same results and the same failures", {
  skip_if_not_installed("analysistools")
  skip_if_not_installed("srvyr")
  syn <- small_synthetic(seed = 5)
  d <- syn$data
  loa <- expand_loa(syn$loa, c(NA, "admin1"))
  with_survey_options("fail", FALSE, {
    des <- srvyr::as_survey_design(d, ids = cluster_id, weights = weight, strata = strata)
    old <- old_analysis_safe(des, loa)
    new <- create_analysis_fast(des, loa = loa, on_error = "skip")
  })
  expect_true(length(old$errors) > 0)
  expect_identical(sort(new$skipped$row_index), sort(old$errors))
  expect_true(compare_analysis_results(old, new, tolerance = 1e-8)$equivalent)
})

test_that("single-indicator functions match analysistools", {
  skip_if_not_installed("analysistools")
  skip_if_not_installed("srvyr")
  syn <- small_synthetic(seed = 3)
  with_survey_options("adjust", TRUE, {
    des <- srvyr::as_survey_design(syn$data, weights = weight, strata = strata)
    pairs <- list(
      list(suppressWarnings(analysistools::create_analysis_mean(des, "admin1", "income", .9)),
           create_analysis_mean_fast(des, "admin1", "income", .9)),
      list(suppressWarnings(suppressMessages(analysistools::create_analysis_median(des, "pop_group", "hh_size", .95))),
           create_analysis_median_fast(des, "pop_group", "hh_size", .95)),
      list(suppressWarnings(suppressMessages(analysistools::create_analysis_prop_select_one(des, "admin1", "so_001", .95))),
           create_analysis_prop_select_one_fast(des, "admin1", "so_001", .95)),
      list(suppressWarnings(analysistools::create_analysis_prop_select_multiple(des, NA, "sm_001", .95)),
           create_analysis_prop_select_multiple_fast(des, NA, "sm_001", .95)),
      list(suppressWarnings(analysistools::create_analysis_ratio(des, "admin1", "num_children_school", "num_children")),
           create_analysis_ratio_fast(des, "admin1", "num_children_school", "num_children"))
    )
  })
  for (p in pairs) expect_true(compare_analysis_results(p[[1]], p[[2]])$equivalent)
})
