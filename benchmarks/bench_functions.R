# Benchmark helpers ----------------------------------------------------------
#
# Every benchmark case runs in a fresh R session (callr) so that timings and
# memory are not polluted by earlier runs. Memory = peak R heap ("max used"
# of gc() after gc(reset = TRUE)) minus the heap in use before the analysis.

bench_project_dir <- function() {
  d <- getwd()
  if (file.exists(file.path(d, "DESCRIPTION"))) return(normalizePath(d, winslash = "/"))
  if (file.exists(file.path(d, "..", "DESCRIPTION"))) return(normalizePath(file.path(d, ".."), winslash = "/"))
  stop("Run the benchmarks from the project root (new_msna_analysis).")
}

# Synthetic dataset for a given size, cached on disk (generated once).
bench_dataset <- function(n, n_filler, cache_dir, seed = 2026) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(cache_dir, sprintf("synthetic_n%d_f%d_s%d.rds", n, n_filler, seed))
  if (!file.exists(f)) {
    syn <- fastmsna::make_synthetic_msna(
      n = n, n_select_one = 10, n_select_multiple = 4, n_numeric = 3,
      n_filler = n_filler, seed = seed
    )
    saveRDS(syn, f, compress = FALSE)
  }
  f
}

# Expand the base LOA over disaggregations.
bench_loa <- function(base_loa, group_vars, types = NULL) {
  if (!is.null(types)) base_loa <- base_loa[base_loa$analysis_type %in% types, , drop = FALSE]
  out <- do.call(rbind, lapply(group_vars, function(g) {
    l <- base_loa
    l$group_var <- g
    # as run_disaggregated_analysis(): never group an indicator by itself
    if (!is.na(g)) l <- l[is.na(l$analysis_var) | !l$analysis_var %in% trimws(strsplit(g, ",")[[1]]), ]
    l
  }))
  out$level <- 0.9
  out
}

# Run one case in a fresh session. tool = "old" (analysistools) or "new".
run_case <- function(tool, data_file, group_vars, types = NULL, design_type = "2scs",
                     timeout = 3600, project_dir = bench_project_dir()) {
  f <- function(tool, data_file, group_vars, types, design_type, project_dir, bench_file) {
    suppressPackageStartupMessages({
      library(srvyr)
      pkgload::load_all(project_dir, quiet = TRUE, export_all = FALSE)
      if (tool == "old") library(analysistools)
    })
    source(bench_file)
    options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE,
            dplyr.summarise.inform = FALSE)
    syn <- readRDS(data_file)
    loa <- bench_loa(syn$loa, group_vars, types)
    d <- syn$data
    design <- if (design_type == "2scs") {
      srvyr::as_survey_design(d, weights = weight, strata = strata)
    } else {
      srvyr::as_survey_design(d)
    }
    invisible(gc(reset = TRUE))
    base_mb <- sum(gc()[, 2])
    invisible(gc(reset = TRUE))
    t <- system.time({
      res <- if (tool == "old") {
        suppressWarnings(suppressMessages(analysistools::create_analysis(design, loa = loa)))
      } else {
        fastmsna::create_analysis_fast(design, loa = loa, on_error = "skip")
      }
    })
    peak_mb <- sum(gc()[, 6])
    list(
      seconds = unname(t[["elapsed"]]),
      peak_mb = peak_mb - base_mb,
      n_rows = nrow(d),
      n_cols = ncol(d),
      loa_rows = nrow(loa),
      result_rows = nrow(res$results_table),
      results = res$results_table
    )
  }
  out <- tryCatch(
    callr::r(f, args = list(tool, data_file, group_vars, types, design_type, project_dir,
                            file.path(project_dir, "benchmarks", "bench_functions.R")),
             timeout = timeout),
    error = function(e) {
      list(seconds = NA_real_, peak_mb = NA_real_, error = conditionMessage(e))
    }
  )
  out
}
