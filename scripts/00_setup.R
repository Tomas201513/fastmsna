#  -----------------------------------------------------------------------------
#
# Title : fastmsna - setup
#    By : IMPACT NGA
#  Date : 10/2026
#
#  Installs (or loads) the fastmsna package from this project folder and checks
#  the dependencies. fastmsna only needs data.table (already installed with the
#  msna_analysis_main environment); analysistools / srvyr / survey are only
#  needed for validation and benchmarks, writexl for xlsx outputs.
#  -----------------------------------------------------------------------------

fastmsna_path <- "C:/Users/User/Music/new_msna_analysis"

# 1. Dependencies --------------------------------------------------------------
required  <- c("data.table")
suggested <- c("srvyr", "survey", "analysistools", "writexl", "readxl", "tibble", "callr",
               "testthat", "pkgload")
have <- rownames(installed.packages())
missing_required <- setdiff(required, have)
if (length(missing_required)) install.packages(missing_required)
missing_suggested <- setdiff(suggested, have)
if (length(missing_suggested)) {
  message("Optional packages not installed (only needed for validation / benchmarks / xlsx): ",
          paste(missing_suggested, collapse = ", "))
}

# 2. Install fastmsna into your library (once; re-run after updating the code) --
#    upgrade = "never" makes sure no other package is re-installed or updated.
if (!requireNamespace("fastmsna", quietly = TRUE) ||
    utils::packageVersion("fastmsna") < read.dcf(file.path(fastmsna_path, "DESCRIPTION"))[, "Version"]) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::install(fastmsna_path, upgrade = "never", dependencies = FALSE, quiet = TRUE)
  } else {
    install.packages(fastmsna_path, repos = NULL, type = "source")
  }
}
library(fastmsna)

# Alternative without installing (development): pkgload::load_all(fastmsna_path)

packageVersion("fastmsna")
