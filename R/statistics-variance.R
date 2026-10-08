# Variance engine -----------------------------------------------------------
#
# Taylor-linearisation variance of domain estimators, for all domains (groups)
# and all indicators (matrix columns) at once. It reproduces survey's
# svyrecvar() -> onestage() -> onestrat() for a one-stage (ultimate cluster)
# stratified design, including:
#   * domain estimation on a physically subset design that keeps the
#     full-design number of PSUs per stratum (missing PSUs padded with 0);
#   * finite population correction per stratum;
#   * survey.lonely.psu = "fail" / "remove" / "certainty" / "adjust" /
#     "average", and survey.adjust.domain.lonely.
#
# Inputs
#   Z : n x k influence values (w * linearised residual / domain total);
#       must be 0 for rows outside the domain or invalid for the column.
#   V : n x k 0/1 matrix, 1 if the row belongs to the (sub)design used for
#       the column (i.e. is in the domain and not missing).
#   gs: group spec (see build_group_spec()).
# Output
#   G x k matrix of variances; attribute "failed" = list(column -> stratum)
#   for columns that hit a lonely PSU under survey.lonely.psu = "fail".

domain_variance <- function(Z, V, gs, des, lonely = lonely_options(des)) {
  k <- ncol(Z)
  if (des$has_clusters) {
    Zu <- rowsum(Z, gs$pc, reorder = TRUE)          # PSU totals within domain
    Vu <- rowsum(V, gs$pc, reorder = TRUE) > 0       # PSU present in subset
    storage.mode(Vu) <- "integer"
    ucell <- gs$pc_cell
  } else {
    Zu <- Z
    Vu <- V
    ucell <- gs$cell
  }
  S1 <- rowsum(Zu, ucell, reorder = TRUE)            # C x k
  M <- rowsum(Vu, ucell, reorder = TRUE)             # PSUs present (nsubset)
  nP <- gs$cell_nPSU
  present <- M > 0

  failed <- NULL
  if (lonely$psu == "fail") {
    bad <- present & (nP == 1) & (gs$cell_f >= 1e-07)
    if (any(bad)) {
      bad_cols <- which(colSums(bad) > 0)
      failed <- lapply(bad_cols, function(j) {
        h <- gs$cell_h[which(bad[, j])[1]]
        des$strata_labels[h]
      })
      names(failed) <- bad_cols
    }
  }

  # centring: stratum mean over the nPSU (padded) PSUs, or - for lonely PSUs
  # under "adjust" - the grand mean of all PSU totals of the domain.
  center <- S1 / nP
  if (lonely$psu == "adjust") {
    tot <- rowsum(S1, gs$cell_g, reorder = TRUE)
    np_sum <- rowsum(present * nP, gs$cell_g, reorder = TRUE)
    rec <- (tot / np_sum)[gs$cell_g, , drop = FALSE]
    use_rec <- !((M > 1) | (nP > 1 & !lonely$adj))
    center[use_rec] <- rec[use_rec]
  }
  dev <- Zu - center[ucell, , drop = FALSE]
  ss <- rowsum(Vu * dev * dev, ucell, reorder = TRUE) + (nP - M) * center * center
  contrib <- gs$cell_scale * ss
  contrib[gs$cell_f < 1e-07, ] <- 0

  if (lonely$psu == "average") {
    na_cell <- present & ((nP == 1) | (M == 1 & nP > 1 & lonely$adj))
    contrib[na_cell] <- NA
  }
  contrib[!present] <- 0

  if (lonely$psu == "average") {
    v <- rowsum(contrib, gs$cell_g, reorder = TRUE, na.rm = TRUE)
    nstrat <- rowsum(present + 0, gs$cell_g, reorder = TRUE)
    nok <- rowsum((present & !is.na(contrib)) + 0, gs$cell_g, reorder = TRUE)
    v <- v * nstrat / nok
  } else {
    v <- rowsum(contrib, gs$cell_g, reorder = TRUE)
  }
  v <- complete_groups(v, gs$G)
  attr(v, "failed") <- failed
  v
}

# rowsum() only returns rows for groups that occur; make sure we always have
# G rows (all groups occur in gid by construction, but be defensive).
complete_groups <- function(m, G) {
  if (nrow(m) == G) {
    rownames(m) <- NULL
    return(m)
  }
  out <- matrix(0, G, ncol(m))
  out[as.integer(rownames(m)), ] <- m
  out
}

# survey::degf() of the (filtered) design: number of PSUs with non-zero
# weight minus number of strata, for each column of R (n x k 0/1 matrix).
design_df <- function(R, des) {
  Rw <- R * (des$w != 0)
  npsu <- if (des$has_clusters) {
    colSums(rowsum(Rw, des$cluster, reorder = TRUE) > 0)
  } else {
    colSums(Rw)
  }
  nstr <- colSums(rowsum(Rw, des$strata, reorder = TRUE) > 0)
  as.numeric(npsu - nstr)
}

# survey::degf() of each domain subset (used by svyquantile): G x k matrix.
group_df <- function(V, gs, des) {
  Vw <- V * (des$w != 0)
  npsu <- if (des$has_clusters) {
    rowsum((rowsum(Vw, gs$pc, reorder = TRUE) > 0) + 0, gs$pc_g, reorder = TRUE)
  } else {
    rowsum(Vw, gs$gid, reorder = TRUE)
  }
  nstr <- rowsum((rowsum(Vw, gs$cell, reorder = TRUE) > 0) + 0, gs$cell_g, reorder = TRUE)
  complete_groups(npsu, gs$G) - complete_groups(nstr, gs$G)
}

# confint() of a svystat with a t distribution: est + se * qt(a, df).
ci_bounds <- function(est, se, level, df) {
  a <- (1 - level) / 2
  lo_fac <- suppressWarnings(stats::qt(a, df = df))
  up_fac <- suppressWarnings(stats::qt(1 - a, df = df))
  list(low = est + se * lo_fac, upp = est + se * up_fac)
}
