#  -----------------------------------------------------------------------------
#
# Title : Real-data benchmark and validation (NGA MSNA 2026 Round 1, 2SCS)
#    By : IMPACT NGA
#  Date : 10/2026
#
#  Reproduces the data preparation of "Step 2 - main analysis - v2.R" for one
#  sheet, then runs analysistools::create_analysis() (old) and
#  fastmsna::create_analysis_fast() (new) for the selected disaggregations,
#  each in a fresh R session, and compares the results row by row.
#
#  Only timings, memory and comparison statistics are written (to
#  outputs/benchmarks/); no survey data is copied into this project.
#  The old tool takes ~10-20 minutes per disaggregation on the hh sheet.
#  -----------------------------------------------------------------------------

project_dir <- "C:/Users/User/Music/new_msna_analysis"
# Local path of the msna_analysis_main project: set it here or with the
# environment variable MSNA_ROOT (e.g. in ~/.Renviron).
msna_root <- Sys.getenv("MSNA_ROOT", "path/to/msna_analysis_main")
if (!dir.exists(msna_root)) {
  stop("Set msna_root (or the MSNA_ROOT environment variable) to your local msna_analysis_main folder.")
}
data_file <- file.path(msna_root, "outputs/datasets/2026-10-04_pre_analysis_dataset_2-STAGE_CLUSTER_SAMPLING.xlsx")
dap_file <- file.path(msna_root, "resources/data_analysis_plan.xlsx")

sheet_index <- 1                 # 1 = "hh data" (order of sheets_to_analyze in Step 2)
# disaggregations to benchmark (hh: 1 = MSNA-wide, 3 = State); override with FASTMSNA_DISAGG_INDEX="1,3"
disagg_index <- as.integer(strsplit(Sys.getenv("FASTMSNA_DISAGG_INDEX", "1,3"), ",")[[1]])
run_old <- Sys.getenv("FASTMSNA_RUN_OLD", "TRUE") == "TRUE"   # FALSE = only time the new tool
old_timeout <- 3 * 3600          # seconds per old run
cached_sheet_rds <- Sys.getenv("FASTMSNA_REAL_SHEET_RDS", "")  # optional pre-read sheet (text) as .rds

out_dir <- file.path(project_dir, "outputs", "benchmarks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(project_dir, "benchmarks", "bench_functions.R"))
suppressPackageStartupMessages(pkgload::load_all(project_dir, quiet = TRUE))

sheets_to_analyze <- c("hh data", "roster data", "health data", "edu data", "nut data", "prot data")
sh <- sheets_to_analyze[sheet_index]

# readxl fails on very long SharePoint paths: work on temporary copies
tmp <- file.path(tempdir(), "fastmsna_real")
dir.create(tmp, showWarnings = FALSE)
file.copy(dap_file, file.path(tmp, "dap.xlsx"), overwrite = TRUE)

prep_file <- file.path(tmp, paste0("prepared_", gsub(" ", "_", sh), ".rds"))
suppressPackageStartupMessages({ library(readxl); library(dplyr) })
if (nzchar(cached_sheet_rds) && file.exists(cached_sheet_rds)) {
  raw <- readRDS(cached_sheet_rds)
} else {
  file.copy(data_file, file.path(tmp, "data.xlsx"), overwrite = TRUE)
  sheet_map <- c("hh data" = "^hh data", "roster data" = "^roster data", "edu data" = "^education data",
                 "health data" = "^health data", "nut data" = "^nutrition data", "prot data" = "^protection data")
  xl_sheet <- grep(sheet_map[[sh]], excel_sheets(file.path(tmp, "data.xlsx")), value = TRUE)[1]
  raw <- read_excel(file.path(tmp, "data.xlsx"), sheet = xl_sheet, col_types = "text")
}

# ---- same preparation as Step 2 v2 ------------------------------------------
my_loa <- read_excel(file.path(tmp, "dap.xlsx"), sheet = "final") %>% mutate(level = as.numeric(level))
disagg <- read_excel(file.path(tmp, "dap.xlsx"), sheet = "disaggregation_2scs") %>%
  pull(sheet_index) %>% na.omit() %>% as.character()
data_with_indicators <- raw
if (sh == "nut data") {
  data_with_indicators <- data_with_indicators %>%
    filter(!is.na(as.numeric(nut_ind_under5_age_months)) & as.numeric(nut_ind_under5_age_months) >= 0 &
             as.numeric(nut_ind_under5_age_months) <= 59) %>%
    mutate(exclusive_breastfeeding = as.numeric(exclusive_breastfeeding),
           under6_month_age = as.numeric(under6_month_age))
}
if ("hesper_clean_women" %in% names(data_with_indicators)) {
  data_with_indicators <- data_with_indicators %>%
    mutate(hesper_clean_women = if_else(hesper_clean_women == "na", "not_applicable", hesper_clean_women))
}
only_nas_lgl <- vapply(data_with_indicators, function(x) all(is.na(x)), logical(1))
data_main <- data_with_indicators[, !only_nas_lgl]
only_nas <- names(only_nas_lgl)[only_nas_lgl]
my_loa_sh <- my_loa %>% filter(sheet == sh, !analysis_var %in% only_nas)
numeric_vars <- my_loa_sh %>% filter(analysis_type %in% c("mean", "median")) %>% pull(analysis_var)
data_main <- data_main %>%
  mutate(across(all_of(numeric_vars), as.numeric)) %>%
  mutate(weight = as.numeric(weight))
saveRDS(list(data = data_main, loa = my_loa_sh, disagg = disagg), prep_file, compress = FALSE)
rm(raw, data_with_indicators)

cat(sprintf("Sheet '%s': %d rows x %d columns, %d LOA rows, disaggregations: %s\n", sh, nrow(data_main),
            ncol(data_main), nrow(my_loa_sh), paste(disagg[disagg_index], collapse = " | ")))

# ---- one run in a fresh session ---------------------------------------------
real_case <- function(tool, gvar) {
  f <- function(tool, gvar, prep_file, project_dir) {
    suppressPackageStartupMessages({
      library(srvyr)
      pkgload::load_all(project_dir, quiet = TRUE, export_all = FALSE)
      if (tool == "old") library(analysistools)
    })
    options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE,
            dplyr.summarise.inform = FALSE, scipen = 999)
    p <- readRDS(prep_file)
    # same expansion as run_disaggregated_analysis(): group_var = gvar, analysis_var != group_var
    loa_sub <- p$loa
    loa_sub$group_var <- gvar
    loa_sub <- loa_sub[!is.na(loa_sub$analysis_var) & loa_sub$analysis_var != gvar, ]
    design <- srvyr::as_survey_design(p$data, weights = "weight", strata = "Strata")
    invisible(gc(reset = TRUE))
    base_mb <- sum(gc()[, 2])
    invisible(gc(reset = TRUE))
    err <- NULL
    t <- system.time({
      res <- if (tool == "old") {
        tryCatch(suppressWarnings(suppressMessages(analysistools::create_analysis(design, loa = loa_sub))),
                 error = function(e) { err <<- conditionMessage(e); NULL })
      } else {
        fastmsna::create_analysis_fast(design, loa = loa_sub, on_error = "skip")
      }
    })
    peak <- sum(gc()[, 6]) - base_mb
    list(seconds = unname(t[["elapsed"]]), peak_mb = peak, loa_rows = nrow(loa_sub),
         n_rows = nrow(p$data), n_cols = ncol(p$data), error = err,
         results = if (is.null(res)) NULL else res$results_table,
         skipped = if (tool == "new") res$skipped else NULL)
  }
  tryCatch(callr::r(f, args = list(tool, gvar, prep_file, project_dir),
                    timeout = if (tool == "old") old_timeout else 3600),
           error = function(e) list(seconds = NA, peak_mb = NA, error = conditionMessage(e)))
}

rows <- list()
for (idx in disagg_index) {
  gvar <- disagg[idx]
  message(sprintf("[%s] disaggregation %d '%s': new tool ...", format(Sys.time(), "%H:%M"), idx, gvar))
  new <- real_case("new", gvar)
  old <- if (run_old) {
    message(sprintf("[%s] disaggregation %d '%s': old tool (slow) ...", format(Sys.time(), "%H:%M"), idx, gvar))
    real_case("old", gvar)
  } else list(seconds = NA, peak_mb = NA)
  cmp <- if (!is.null(old$results) && !is.null(new$results)) {
    compare_analysis_results(old$results, new$results, tolerance = 1e-8)
  } else NULL
  rows[[length(rows) + 1]] <- data.frame(
    sheet = sh, disagg_index = idx, disaggregation = gvar,
    n_rows = new$n_rows, n_cols = new$n_cols, loa_rows = new$loa_rows,
    result_rows_old = if (is.null(cmp)) NA else cmp$n_old,
    result_rows_new = if (is.null(new$results)) NA else nrow(new$results),
    seconds_old = round(old$seconds, 2), seconds_new = round(new$seconds, 2),
    speedup = round(old$seconds / new$seconds, 1),
    peak_mb_old = round(old$peak_mb, 1), peak_mb_new = round(new$peak_mb, 1),
    value_mismatches = if (is.null(cmp)) NA else nrow(cmp$mismatches),
    keys_only_old = if (is.null(cmp)) NA else length(cmp$only_in_old),
    keys_only_new = if (is.null(cmp)) NA else length(cmp$only_in_new),
    same_row_order = if (is.null(cmp)) NA else cmp$same_row_order,
    max_rel_diff = if (is.null(cmp)) NA else max(cmp$summary$max_rel_diff),
    skipped_new = if (is.null(new$skipped)) NA else nrow(new$skipped),
    old_error = if (is.null(old$error)) "" else old$error
  )
  print(rows[[length(rows)]])
  if (!is.null(cmp) && !cmp$equivalent) {
    utils::write.csv(cmp$mismatches, file.path(out_dir, sprintf("real_data_mismatches_%s_%02d.csv", gsub(" ", "_", sh), idx)),
                     row.names = FALSE)
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(out_dir, sprintf("real_data_benchmark_%s.csv", gsub(" ", "_", sh))), row.names = FALSE)
print(res)
