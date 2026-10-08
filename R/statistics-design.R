# Design layer -------------------------------------------------------------
#
# A `fastmsna_design` holds only what the variance engine needs, as plain
# integer/double vectors: weights, stratum codes, PSU codes, the full-design
# number of PSUs per stratum (survey's fpc$sampsize) and the finite population
# correction. The dataset itself is kept by reference (never copied); columns
# are pulled one at a time when an indicator needs them.

#' Prepare a survey design for fast analysis
#'
#' Converts a `srvyr`/`survey` design object, or a data frame plus the names of
#' its design columns, into a light-weight `fastmsna_design`. Build it **once**
#' per dataset and re-use it for every call to [create_analysis_fast()]: this
#' is where most of the old per-call overhead disappears.
#'
#' Supported designs are the linearisation designs used in MSNAs: weights,
#' optional strata, optional (first-stage) clusters and optional fpc. Replicate
#' weight, two-phase, PPS and calibrated/post-stratified designs are rejected
#' with an error (use `analysistools` for those).
#'
#' @param design A `tbl_svy`/`survey.design2` object (from
#'   `srvyr::as_survey_design()` or `survey::svydesign()`), a data frame, or an
#'   existing `fastmsna_design` (returned unchanged).
#' @param weights,strata,ids,fpc Only used when `design` is a data frame:
#'   column names (strings) for the sampling weights, strata, first-stage
#'   cluster ids and finite population correction. `NULL` means "none"
#'   (weights = 1, one stratum, every row its own PSU, no fpc), exactly like
#'   `srvyr::as_survey_design()`. `fpc` follows `survey`: values > 1 are
#'   population sizes (number of PSUs), values <= 1 are sampling fractions.
#' @param nest Only for data frames: if `TRUE`, cluster ids are re-coded to be
#'   unique within strata (as `survey::svydesign(nest = TRUE)`).
#' @param lonely_psu,adjust_domain_lonely Override the `survey.lonely.psu` and
#'   `survey.adjust.domain.lonely` options for this design. `NULL` (default)
#'   reads the options at analysis time, like `survey` does.
#'
#' @return An object of class `fastmsna_design`.
#' @export
#'
#' @examples
#' df <- data.frame(y = c(1, 4, 2, 8, 5, 7), s = rep(c("a", "b"), 3),
#'                  w = c(1, 2, 1, 2, 1, 2))
#' des <- as_fast_design(df, weights = "w", strata = "s")
#' des
as_fast_design <- function(design,
                           weights = NULL,
                           strata = NULL,
                           ids = NULL,
                           fpc = NULL,
                           nest = FALSE,
                           lonely_psu = NULL,
                           adjust_domain_lonely = NULL) {
  if (inherits(design, "fastmsna_design")) {
    if (!is.null(lonely_psu)) design$lonely_psu <- lonely_psu
    if (!is.null(adjust_domain_lonely)) design$adjust_domain_lonely <- adjust_domain_lonely
    return(design)
  }
  unsupported <- c("svyrep.design", "twophase", "twophase2", "pps",
                   "DBIsvydesign", "ODBCsvydesign", "svyimputationList")
  if (inherits(design, unsupported)) {
    stop("fastmsna supports linearisation designs (weights, strata, clusters, fpc) only. ",
         "Replicate-weight, two-phase, PPS and database-backed designs are not supported: ",
         "use analysistools/srvyr for these.", call. = FALSE)
  }
  out <- if (inherits(design, "survey.design2")) {
    design_from_survey(design)
  } else if (is.data.frame(design)) {
    design_from_data(design, weights = weights, strata = strata, ids = ids, fpc = fpc, nest = nest)
  } else {
    stop("`design` must be a srvyr/survey design, a data frame or a fastmsna_design.", call. = FALSE)
  }
  out$lonely_psu <- lonely_psu
  out$adjust_domain_lonely <- adjust_domain_lonely
  out
}

design_from_survey <- function(design) {
  if (!is.null(design$postStrata)) {
    stop("Calibrated / post-stratified / raked designs are not supported by fastmsna ",
         "(their variance needs the calibration residuals). Use analysistools for these.",
         call. = FALSE)
  }
  if (!is.null(design$pps) && !isFALSE(design$pps)) {
    stop("PPS designs (svydesign(pps = ...)) are not supported by fastmsna.", call. = FALSE)
  }
  popsize <- if (is.null(design$fpc$popsize)) NULL else design$fpc$popsize[, 1]
  if (NCOL(design$cluster) > 1 && !is.null(popsize) &&
      !isTRUE(getOption("survey.ultimate.cluster"))) {
    stop("Multi-stage designs with a finite population correction are not supported ",
         "(later-stage variance contributions). Drop the fpc or use analysistools.",
         call. = FALSE)
  }
  new_fast_design(
    variables = design$variables,
    w = 1 / as.numeric(design$prob),
    strata = design$strata[[1]],
    cluster = design$cluster[[1]],
    sampsize = design$fpc$sampsize[, 1],
    popsize = popsize,
    source = "survey design"
  )
}

design_from_data <- function(data, weights, strata, ids, fpc, nest) {
  n <- nrow(data)
  get_col <- function(x, what) {
    if (is.null(x)) return(NULL)
    if (is.character(x) && length(x) == 1L) {
      if (!x %in% names(data)) stop("Cannot find the ", what, " column '", x, "' in the data.", call. = FALSE)
      return(data[[x]])
    }
    if (length(x) != n) stop("`", what, "` must be a column name or a vector of length nrow(data).", call. = FALSE)
    x
  }
  w <- get_col(weights, "weights")
  if (is.null(w)) {
    w <- rep(1, n)
  } else {
    w <- suppressWarnings(as.numeric(w))
    if (anyNA(w)) stop("missing (or non-numeric) values in weights.", call. = FALSE)
    # survey stores probabilities (1 / w) and uses 1 / prob as weights: mirror
    # that round trip so results are bit-for-bit comparable with srvyr designs.
    w <- 1 / (1 / w)
  }
  s <- get_col(strata, "strata")
  if (is.null(s)) s <- rep(1L, n)
  if (anyNA(s)) stop("missing values in strata.", call. = FALSE)
  cl <- get_col(ids, "ids")
  if (is.null(cl)) cl <- seq_len(n)
  if (anyNA(cl)) stop("missing values in ids.", call. = FALSE)
  if (isTRUE(nest) && !is.null(strata)) cl <- paste(s, cl, sep = "\r")
  popsize <- NULL
  f <- get_col(fpc, "fpc")
  if (!is.null(f)) {
    f <- as.numeric(f)
    if (anyNA(f)) stop("missing values in fpc.", call. = FALSE)
    # sampsize is filled in by new_fast_design(); fractions are converted there.
    popsize <- f
  }
  new_fast_design(variables = data, w = w, strata = s, cluster = cl,
                  sampsize = NULL, popsize = popsize, fpc_is_raw = !is.null(f),
                  source = "data frame")
}

new_fast_design <- function(variables, w, strata, cluster, sampsize, popsize,
                            fpc_is_raw = FALSE, source = "") {
  n <- length(w)
  s <- frankv(strata, ties.method = "dense")
  H <- if (n) max(s) else 0L
  cl <- frankv(cluster, ties.method = "dense")
  P <- if (n) max(cl) else 0L

  first_psu <- !duplicated(cl)
  psu_stratum <- integer(P)
  psu_stratum[cl[first_psu]] <- s[first_psu]
  if (any(psu_stratum[cl] != s)) {
    stop("Clusters (PSUs) are not nested within strata. Re-code the ids or use nest = TRUE.",
         call. = FALSE)
  }

  first_str <- !duplicated(s)
  if (is.null(sampsize)) {
    nPSU <- tabulate(psu_stratum, nbins = H)
  } else {
    nPSU <- integer(H)
    nPSU[s[first_str]] <- as.integer(sampsize[first_str])
  }

  fpc_f <- rep(1, H)
  if (!is.null(popsize)) {
    pop <- numeric(H)
    pop[s[first_str]] <- popsize[first_str]
    if (fpc_is_raw && !any(pop > 1)) pop <- nPSU / pop # sampling fractions -> pop size
    if (any(pop < nPSU)) stop("FPC implies >100% sampling in some strata.", call. = FALSE)
    fpc_f <- ifelse(pop == Inf, 1, (pop - nPSU) / pop)
  }
  scale <- ifelse(nPSU > 1, fpc_f * nPSU / (nPSU - 1), fpc_f)

  structure(
    list(
      variables = variables,
      n = n,
      w = as.numeric(w),
      strata = s,
      strata_labels = as.character(strata[first_str][order(s[first_str])]),
      H = H,
      cluster = cl,
      P = P,
      psu_stratum = psu_stratum,
      has_clusters = P < n,
      nPSU = nPSU,
      fpc_f = fpc_f,
      scale = scale,
      cache = new.env(parent = emptyenv()),
      source = source,
      lonely_psu = NULL,
      adjust_domain_lonely = NULL
    ),
    class = "fastmsna_design"
  )
}

#' @export
print.fastmsna_design <- function(x, ...) {
  cat("<fastmsna_design> from ", x$source, "\n", sep = "")
  cat("  rows: ", x$n, " | columns: ", length(x$variables), "\n", sep = "")
  cat("  strata: ", x$H, " | PSUs: ", x$P,
      if (x$has_clusters) " (clustered)" else " (one row per PSU)", "\n", sep = "")
  cat("  weights: ", if (all(x$w == 1)) "none (all 1)" else
    sprintf("min %.4g / max %.4g / sum %.6g", min(x$w), max(x$w), sum(x$w)), "\n", sep = "")
  cat("  fpc: ", if (all(x$fpc_f == 1)) "none" else "yes", "\n", sep = "")
  lp <- lonely_options(x)
  cat("  lonely PSU: ", lp$psu, " | adjust domain lonely: ", lp$adj, "\n", sep = "")
  invisible(x)
}

# Resolve the lonely-PSU options exactly like survey does (options read at
# computation time unless overridden on the design).
lonely_options <- function(des) {
  psu <- des$lonely_psu
  if (is.null(psu)) psu <- getOption("survey.lonely.psu")
  if (is.null(psu)) psu <- "fail"
  psu <- match.arg(psu, c("fail", "adjust", "remove", "certainty", "average"))
  adj <- des$adjust_domain_lonely
  if (is.null(adj)) adj <- isTRUE(getOption("survey.adjust.domain.lonely"))
  list(psu = psu, adj = isTRUE(adj))
}

#' Clear the column cache of a fast design
#'
#' Converted columns (numeric versions of text columns, select-one codes,...)
#' are cached on the design so that running many disaggregations does not
#' convert the same column again. Call this to free that memory.
#'
#' @param design A `fastmsna_design`.
#' @return The design, invisibly.
#' @export
clear_design_cache <- function(design) {
  stopifnot(inherits(design, "fastmsna_design"))
  rm(list = ls(design$cache, all.names = TRUE), envir = design$cache)
  invisible(design)
}
