#  -----------------------------------------------------------------------------
#
# Title : Benchmarks - analysistools (old) vs fastmsna (new)
#    By : IMPACT NGA
#  Date : 10/2026
#
#  Synthetic MSNA-like datasets (make_synthetic_msna(): strata = admin2 x
#  population group, stratum weights, select-one / select-multiple / numeric
#  questions with skip-logic NAs, ratio pairs, 300 unused text columns to mimic
#  the width of real MSNA datasets).
#
#  Part A - full list of analysis (31 indicators x 5 disaggregations = 155 LOA
#           rows) at increasing dataset sizes.
#  Part B - one analysis type at a time (disaggregation "admin1, pop_group",
#           18 groups) to see which functions gain the most.
#
#  Every case runs in a fresh R session (callr). Memory = peak R heap during
#  the analysis. Results of old and new are compared for every case.
#  Run from the project root:  Rscript benchmarks/run_benchmarks.R
#  Env vars: FASTMSNA_BENCH_SIZES (default "1000,5000,20000,100000"),
#            FASTMSNA_BENCH_OLD_MAX_N (largest n for the old full-LOA run, default 20000)
#  -----------------------------------------------------------------------------

source(file.path("benchmarks", "bench_functions.R"))
project_dir <- bench_project_dir()
suppressPackageStartupMessages(pkgload::load_all(project_dir, quiet = TRUE))

cache_dir <- file.path(project_dir, "benchmarks", "cache")
out_dir <- file.path(project_dir, "outputs", "benchmarks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sizes <- as.integer(strsplit(Sys.getenv("FASTMSNA_BENCH_SIZES", "1000,5000,20000,100000"), ",")[[1]])
old_max_n <- as.integer(Sys.getenv("FASTMSNA_BENCH_OLD_MAX_N", "20000"))
n_filler <- 300
group_vars <- c(NA, "admin1", "pop_group", "hoh_gender", "admin1, pop_group")
types <- c("prop_select_one", "prop_select_multiple", "mean", "median", "ratio")

stamp <- function(...) message(format(Sys.time(), "%H:%M:%S "), sprintf(...))
check_equal <- function(old, new) {
  if (is.null(old$results) || is.null(new$results)) return(NA)
  compare_analysis_results(old$results, new$results, tolerance = 1e-8)$equivalent
}
row_of <- function(part, n, what, old, new, g) {
  data.frame(
    part = part, n_rows = n, n_cols = new$n_cols, analysis = what, disaggregations = g,
    loa_rows = new$loa_rows, result_rows = new$result_rows,
    seconds_old = round(old$seconds, 2), seconds_new = round(new$seconds, 3),
    speedup = round(old$seconds / new$seconds, 1),
    peak_mb_old = round(old$peak_mb, 1), peak_mb_new = round(new$peak_mb, 1),
    row_indicators_per_s_old = round(n * new$loa_rows / old$seconds),
    row_indicators_per_s_new = round(n * new$loa_rows / new$seconds),
    results_equivalent = check_equal(old, new),
    old_status = if (!is.null(old$error)) gsub("\\s+", " ", old$error) else if (is.na(old$seconds)) "not run" else "ok"
  )
}

rows <- list()
not_run <- list(seconds = NA_real_, peak_mb = NA_real_, error = "not run (too slow, see part B)")

# ---- Part A: full LOA -------------------------------------------------------
for (n in sizes) {
  f <- bench_dataset(n, n_filler, cache_dir)
  stamp("A n=%d: new", n)
  new <- run_case("new", f, group_vars)
  old <- if (n <= old_max_n) {
    stamp("A n=%d: old (slow)", n)
    run_case("old", f, group_vars, timeout = 4 * 3600)
  } else not_run
  rows[[length(rows) + 1]] <- row_of("A_full_loa", n, "all types", old, new, length(group_vars))
  print(rows[[length(rows)]][, c("n_rows", "seconds_old", "seconds_new", "speedup", "peak_mb_old",
                                  "peak_mb_new", "results_equivalent")])
  utils::write.csv(do.call(rbind, rows), file.path(out_dir, "benchmark_results.csv"), row.names = FALSE)
}

# ---- Part B: per analysis type ----------------------------------------------
for (n in intersect(c(20000L, 100000L), sizes)) {
  f <- bench_dataset(n, n_filler, cache_dir)
  for (tp in types) {
    stamp("B n=%d %s: new", n, tp)
    new <- run_case("new", f, "admin1, pop_group", types = tp)
    stamp("B n=%d %s: old", n, tp)
    old <- run_case("old", f, "admin1, pop_group", types = tp, timeout = 3 * 3600)
    rows[[length(rows) + 1]] <- row_of("B_per_function", n, tp, old, new, 1)
    print(rows[[length(rows)]][, c("n_rows", "analysis", "seconds_old", "seconds_new", "speedup",
                                    "results_equivalent")])
    utils::write.csv(do.call(rbind, rows), file.path(out_dir, "benchmark_results.csv"), row.names = FALSE)
  }
}

res <- do.call(rbind, rows)
utils::write.csv(res, file.path(out_dir, "benchmark_results.csv"), row.names = FALSE)

# Estimated old full-LOA time where it was not run: per-row cost of part B at the
# same n (sum of per-type seconds / LOA rows) x full LOA rows.
for (i in which(res$part == "A_full_loa" & is.na(res$seconds_old))) {
  b <- res[res$part == "B_per_function" & res$n_rows == res$n_rows[i] & !is.na(res$seconds_old), ]
  if (nrow(b)) {
    est <- sum(b$seconds_old) / sum(b$loa_rows) * res$loa_rows[i]
    res$old_status[i] <- sprintf("not run; estimated from part B: %.0f s (%.1f h)", est, est / 3600)
  }
}
utils::write.csv(res, file.path(out_dir, "benchmark_results.csv"), row.names = FALSE)

md <- c(
  "# Benchmark results: analysistools (old) vs fastmsna (new)",
  "",
  paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " on ", Sys.info()[["sysname"]], ", R ",
         getRversion(), ", ", parallel::detectCores(), " cores. Each case in a fresh R session."),
  "",
  "`row-indicators/s` = dataset rows x LOA rows / seconds. Memory = peak R heap during the analysis.",
  "",
  "| part | rows | analysis | LOA rows | old (s) | new (s) | speed-up | old peak MB | new peak MB | old row-ind/s | new row-ind/s | results equal | note |",
  "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %s | %s | %d | %s | %.2f | %s | %s | %.0f | %s | %s | %s | %s |",
          res$part, format(res$n_rows, big.mark = ","), res$analysis, res$loa_rows,
          ifelse(is.na(res$seconds_old), "-", format(res$seconds_old, nsmall = 1)), res$seconds_new,
          ifelse(is.na(res$speedup), "-", paste0(res$speedup, "x")),
          ifelse(is.na(res$peak_mb_old), "-", res$peak_mb_old), res$peak_mb_new,
          ifelse(is.na(res$row_indicators_per_s_old), "-", format(res$row_indicators_per_s_old, big.mark = ",")),
          format(res$row_indicators_per_s_new, big.mark = ","),
          ifelse(is.na(res$results_equivalent), "-", res$results_equivalent),
          ifelse(res$old_status == "ok", "", res$old_status))
)
writeLines(md, file.path(out_dir, "BENCHMARKS.md"))
print(res)
