#  -----------------------------------------------------------------------------
#
# Title : Real-data full-sheet timing (new tool only)
#    By : IMPACT NGA
#  Date : 10/2026
#
#  Runs run_disaggregated_analysis_fast() over ALL disaggregations of one sheet,
#  exactly as "Step 2 - main analysis - fast.R" does, saving the results to a
#  temporary folder (nothing is written to msna_analysis_main). Reports the
#  time per disaggregation and the total.
#  Needs the prepared sheet written by benchmark_real_data.R, or set
#  FASTMSNA_REAL_PREPARED to an .rds with list(data, loa, disagg).
#  -----------------------------------------------------------------------------

project_dir <- "C:/Users/User/Music/new_msna_analysis"
prepared <- Sys.getenv("FASTMSNA_REAL_PREPARED", "")
sheet_name <- Sys.getenv("FASTMSNA_SHEET", "hh data")
suppressPackageStartupMessages({
  library(srvyr)
  pkgload::load_all(project_dir, quiet = TRUE)
})
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)

p <- readRDS(prepared)
out <- file.path(tempdir(), "fastmsna_full_sheet")
unlink(out, recursive = TRUE)

design <- srvyr::as_survey_design(p$data, weights = "weight", strata = "Strata")
timings <- list()
t_all <- system.time({
  t_des <- system.time(fd <- as_fast_design(design))
  for (i in seq_along(p$disagg)) {
    t <- system.time(res <- run_disaggregated_analysis_fast(
      sh = sheet_name, loa_sheet = p$loa, disagg_vars = p$disagg, disagg_index = i,
      design = fd, output_dir = out, save_format = "rds", verbose = FALSE
    ))
    timings[[i]] <- data.frame(
      sheet = sheet_name, disagg_index = i, disaggregation = p$disagg[i],
      groups = length(unique(res[[1]]$group_var_value)),
      result_rows = nrow(res[[1]]), seconds = round(t[["elapsed"]], 2)
    )
    message(sprintf("%2d %-62s %6d rows %7.2f s", i, p$disagg[i], nrow(res[[1]]), t[["elapsed"]]))
  }
})
tab <- do.call(rbind, timings)
peak <- sum(gc()[, 6])
cat(sprintf("\nSheet '%s': %d rows x %d cols, %d LOA rows, %d disaggregations\n", sheet_name, nrow(p$data),
            ncol(p$data), nrow(p$loa), length(p$disagg)))
cat(sprintf("design preparation: %.2f s | total incl. saving: %.1f s (%.1f min) | peak R memory: %.0f MB\n",
            t_des[["elapsed"]], t_all[["elapsed"]], t_all[["elapsed"]] / 60, peak))
utils::write.csv(tab, file.path(project_dir, "outputs", "benchmarks",
                                sprintf("real_data_full_sheet_%s.csv", gsub(" ", "_", sheet_name))),
                 row.names = FALSE)
