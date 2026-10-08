#  -----------------------------------------------------------------------------
#
# Title : fastmsna - how to use (examples)
#    By : IMPACT NGA
#  Date : 10/2026
#  -----------------------------------------------------------------------------

library(fastmsna)          # or pkgload::load_all("C:/Users/User/Music/new_msna_analysis")
library(srvyr)

# Same options as msna_analysis_main/src/init.R (fastmsna reads them like survey)
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)

# 1. Some data --------------------------------------------------------------------
syn <- make_synthetic_msna(n = 2000, n_filler = 20)   # synthetic MSNA-like household data
df  <- syn$data
loa <- syn$loa                                         # analysis_type / analysis_var / group_var / level / ratio columns
head(loa)

# 2. Exactly as with analysistools: a srvyr design + a LOA --------------------------
my_design <- as_survey_design(df, weights = weight, strata = strata)

res <- create_analysis_fast(my_design, loa = loa, sm_separator = ".")
res$results_table        # same columns / analysis keys as analysistools::create_analysis()
res$skipped              # analyses that could not be computed (empty here)

# 3. Disaggregations: group_var exactly as in analysistools -------------------------
loa_admin1 <- transform(loa, group_var = "admin1")                 # by admin1
loa_cross  <- transform(loa, group_var = "admin1, pop_group")      # admin1 x population group
res2 <- create_analysis_fast(my_design, loa = rbind(loa, loa_admin1, loa_cross), on_error = "skip")

# 4. Prepare the design once when you run many analyses on the same data -----------
fd <- as_fast_design(my_design)       # weights / strata / PSU counts as plain vectors
fd
res3 <- create_analysis_fast(fd, loa = loa_admin1)

# A data frame works too (no srvyr needed): same arguments as as_survey_design()
res4 <- create_analysis_fast(df, loa = loa, weights = "weight", strata = "strata")

# 5. Single indicators (same signatures as analysistools::create_analysis_*) ------
create_analysis_mean_fast(my_design, group_var = "admin1", analysis_var = "income", level = .9)
create_analysis_median_fast(my_design, group_var = NA, analysis_var = "hh_size")
create_analysis_prop_select_one_fast(my_design, group_var = "pop_group", analysis_var = "setting")
create_analysis_prop_select_multiple_fast(my_design, group_var = NA, analysis_var = "sm_001", sm_separator = ".")
create_analysis_ratio_fast(my_design, group_var = "admin1",
                           analysis_var_numerator = "num_children_school",
                           analysis_var_denominator = "num_children")

# 6. Many disaggregations, one file each (drop-in for run_disaggregated_analysis) --
run_disaggregated_analysis_fast(
  sh = "hh data",
  loa_sheet = loa,
  disagg_vars = c("admin1", "pop_group", "admin1, pop_group"),
  design = fd,
  output_dir = file.path(tempdir(), "fastmsna_example"),
  save_format = "csv"
)

# 7. Compare with the old tool (validation of a migration) --------------------------
# old <- analysistools::create_analysis(my_design, loa = loa_admin1)
# compare_analysis_results(old, res3)
