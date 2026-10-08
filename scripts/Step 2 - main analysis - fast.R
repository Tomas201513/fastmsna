#  -----------------------------------------------------------------------------
#
# Title : Analysis Generator (fastmsna version of "Step 2 - main analysis - v2.R")
#    By : IMPACT NGA
#  Date : 10/2026
#
#  Same inputs, same outputs (one results file per disaggregation in
#  ./outputs/analysis/<sampling>/results table/<sheet>/), so "Step 3 - merge
#  analysis.R" and "Step 4 - formatted analysis.R" run unchanged.
#
#  What changed compared with Step 2 v2 (search for "# FAST:"):
#    1. fastmsna is loaded (after src/init.R).
#    2. the srvyr design is converted ONCE with as_fast_design().
#    3. run_disaggregated_analysis() -> run_disaggregated_analysis_fast()
#       (same arguments). Failing indicators are skipped and logged in
#       <output_dir>/_logs/ instead of failing a whole disaggregation.
#    4. because a sheet now takes seconds/minutes instead of hours, the script
#       can loop over all sheets and all disaggregations in one run.
#
#  Run it from the msna_analysis_main project (working directory = project root)
#  like the other Step scripts.
#  -----------------------------------------------------------------------------

rm(list=ls())

# Resources path -------------------------------------------------------------
kobo_tool <- "./resources/kobo_for_analysis.xlsx"
dap <- "resources/data_analysis_plan.xlsx"
samping_type <- "2scs"  # "quota"  "2scs"

if (samping_type == "quota") {
  clean_data_path <- "./outputs/datasets/2026-10-04_pre_analysis_dataset_REFUGEE_QUOTA_SAMPLING.xlsx"
} else if (samping_type == "2scs") {
  clean_data_path <- "./outputs/datasets/2026-10-04_pre_analysis_dataset_2-STAGE_CLUSTER_SAMPLING.xlsx"
} else {
  stop("Invalid sampling type. Please choose either 'quota' or '2scs'.")
}

# Load src scripts -------------------------------------------------------------
source("src/init.R")

# FAST: load fastmsna. Install it once from GitHub:
#   remotes::install_github("Tomas201513/fastmsna", upgrade = "never")
if (!requireNamespace("fastmsna", quietly = TRUE)) {
  remotes::install_github("Tomas201513/fastmsna", upgrade = "never")
}
library(fastmsna)
# or, without installing: pkgload::load_all("C:/Users/User/Music/new_msna_analysis")

#############################################################################################################
### LOAD DATASET & ANALYSIS PLAN
#############################################################################################################
clean_sheets <- grep("data", excel_sheets(clean_data_path), value = TRUE, ignore.case = TRUE)

label_map <- c(
  "^hh data"     = "hh data",
  "^roster data" = "roster data",
  "^education data"    = "edu data",
  "^health data" = "health data",
  "^nutrition data"    = "nut data",
  "^protection data"   = "prot data"
)

labels <- vapply(clean_sheets, \(x) {
  hit <- names(label_map)[sapply(names(label_map), grepl, x = x)]
  if (length(hit) == 0) stop("No label match for sheet: ", x)
  label_map[[hit]]
}, character(1))

data <- setNames(
  lapply(clean_sheets, \(x) read_excel(clean_data_path, sheet = x, col_types = "text")),
  labels
)

# --- Load DAP ---
my_loa <- read_excel(dap, sheet = "final") %>%
  mutate(level = as.numeric(level))

missing <- find_missing_indicators(kobo_tool, dap, return_details = TRUE)
print(missing)

sheets_to_analyze <- c("hh data", "roster data", "health data", "edu data", "nut data", "prot data")

# FAST: run one, several or all sheets in one go
sheets_to_run <- 1:6          #' e.g. 3, c(1, 3), 1:6
disagg_index_to_run <- NULL   #' NULL = all disaggregations of the sheet, or e.g. c(24, 1)

for (sheet in sheets_to_run) {

  sh <- sheets_to_analyze[sheet]

  if (samping_type == "quota") {
    disagg <- read_excel(dap, sheet = "disaggregation_quota") %>%
      pull(sheet) %>% na.omit() %>% as.character()
  } else {
    disagg <- read_excel(dap, sheet = "disaggregation_2scs") %>%
      pull(sheet) %>% na.omit() %>% as.character()
  }

  data_with_indicators <- data[[sh]]

  #' *the nut data sheet has a special age filter for under-5 indicators*
  if (sh == "nut data") {
    data_with_indicators <- data_with_indicators %>%
      filter(!is.na(as.numeric(nut_ind_under5_age_months)) & as.numeric(nut_ind_under5_age_months) >= 0 & as.numeric(nut_ind_under5_age_months) <= 59) %>%
      mutate(
        exclusive_breastfeeding = as.numeric(exclusive_breastfeeding),
        under6_month_age        = as.numeric(under6_month_age)
      )
  }

  #' *recode the "na" choice to "not_applicable" for hesper_clean_women*
  if ("hesper_clean_women" %in% names(data_with_indicators)) {
    data_with_indicators <- data_with_indicators %>%
      mutate(hesper_clean_women = if_else(hesper_clean_women == "na", "not_applicable", hesper_clean_women))
  }

  only_nas_lgl <- data_with_indicators %>%
    summarise(across(everything(), ~ sum(is.na(.)) == nrow(data_with_indicators))) %>%
    unlist()

  data_main <- data_with_indicators[, !only_nas_lgl]
  only_nas  <- names(only_nas_lgl)[only_nas_lgl]

  # base LOA for this sheet (group_var is set per disaggregation by the runner)
  my_loa_sh <- my_loa %>%
    filter(sheet == sh, !analysis_var %in% only_nas)

  numeric_vars <- my_loa_sh %>%
    filter(analysis_type %in% c("mean", "median")) %>%
    pull(analysis_var)

  data_main <- data_main %>%
    mutate(across(all_of(numeric_vars), as.numeric))

  if (samping_type == "quota") {
    my_design <- srvyr::as_survey_design(data_main)
  } else {
    data_main <- data_main %>%
      mutate(weight = as.numeric(weight))
    my_design <- srvyr::as_survey_design(data_main, weights = "weight", strata = "Strata")
  }

  # FAST: prepare the design once (weights, strata, PSU counts, cached columns)
  my_fast_design <- as_fast_design(my_design)
  print(my_fast_design)

  output_dir <- paste0("./outputs/analysis/", samping_type, "/results table/", sh, "/")

  cat(green(paste0("Running disaggregated analysis for sheet: ", sh, "\n\n",
                   "Number of analysis variables: ", nrow(my_loa_sh), "\n\n",
                   "Disaggregation levels: ", paste(disagg, collapse = ", "), "\n\n",
                   "Output directory: ", output_dir, "\n\n",
                   "Missing indicators: ", paste(only_nas, collapse = ", "), "\n")))

  # FAST: same arguments as run_disaggregated_analysis()
  run_disaggregated_analysis_fast(
    sh           = sh,
    loa_sheet    = my_loa_sh,
    disagg_vars  = disagg,
    disagg_index = if (is.null(disagg_index_to_run)) seq_along(disagg) else disagg_index_to_run,
    design       = my_fast_design,
    output_dir   = output_dir,
    save_format  = "xlsx"
  )
}

##################################################################################
# Column/indicator-specific re-run (one indicator had an error)
#
# vars_to_rerun <- c("ind_age")
# my_loa_fix <- my_loa %>% filter(sheet == sh, analysis_var %in% vars_to_rerun)
# run_disaggregated_analysis_fast(
#   sh = sh, loa_sheet = my_loa_fix, disagg_vars = disagg, design = my_fast_design,
#   output_dir = paste0(output_dir, "specific_rerun/"), save_format = "xlsx"
# )
