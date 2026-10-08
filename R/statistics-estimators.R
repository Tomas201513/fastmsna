# Estimator kernels ---------------------------------------------------------
#
# Every kernel works on n x k matrices (k = indicators / choices / levels) and
# returns G x k matrices (G = groups of the disaggregation). One call replaces
# G x k separate srvyr::summarise() evaluations of the old tool.

# Weighted domain mean (also proportions: Y = 0/1 indicator).
#   Y : n x k numeric, MUST be 0 where V == 0 (see zero_invalid())
#   V : n x k integer 0/1 validity (row in the domain subset for that column)
# survey::svymean(): mean = sum(w y) / sum(w); influence = w (y - mean) / sum(w).
kernel_mean <- function(Y, V, gs, des, lonely) {
  gid <- gs$gid
  WV <- des$w * V
  Wg <- complete_groups(rowsum(WV, gid, reorder = TRUE), gs$G)
  SY <- complete_groups(rowsum(WV * Y, gid, reorder = TRUE), gs$G)
  est <- SY / Wg
  # rows outside the domain have WV = 0, so they get Z = 0 (empty domains:
  # estimate NaN, use 0 for the expansion to avoid NaN * 0)
  e0 <- est
  e0[!is.finite(e0)] <- 0
  iw <- 1 / Wg
  iw[!is.finite(iw)] <- 0
  Z <- WV * (Y - e0[gid, , drop = FALSE]) * iw[gid, , drop = FALSE]
  v <- domain_variance(Z, V, gs, des, lonely)
  list(
    est = est,
    var = v,
    n = complete_groups(rowsum(V, gid, reorder = TRUE), gs$G),
    nw = Wg,
    sum_y = complete_groups(rowsum(Y, gid, reorder = TRUE), gs$G),
    sum_wy = SY,
    failed = attr(v, "failed")
  )
}

# Replace missing values by 0 (validity is carried by V).
zero_na <- function(y) {
  y[is.na(y)] <- 0
  y
}

# Weighted domain ratio sum(w y) / sum(w x) (survey::svyratio()).
# influence = w (y - R x) / sum(w x).
kernel_ratio <- function(Yn, Yd, V, gs, des, lonely) {
  w <- des$w
  Yn[V == 0L] <- 0
  Yd[V == 0L] <- 0
  WV <- w * V
  SN <- complete_groups(rowsum(WV * Yn, gs$gid, reorder = TRUE), gs$G)
  SD <- complete_groups(rowsum(WV * Yd, gs$gid, reorder = TRUE), gs$G)
  R <- SN / SD
  r0 <- R
  r0[!is.finite(r0)] <- 0
  isd <- 1 / SD
  isd[!is.finite(isd)] <- 0
  Z <- WV * (Yn - r0[gs$gid, , drop = FALSE] * Yd) * isd[gs$gid, , drop = FALSE]
  v <- domain_variance(Z, V, gs, des, lonely)
  list(est = R, var = v, failed = attr(v, "failed"))
}

# survey:::qrule_school() on x already sorted ascending (stable) with weights w.
# `sw` is the total weight summed in the ORIGINAL row order, exactly as
# survey:::qs() does (sum(w) before ordering), so knife-edge comparisons
# (cumw <= p * sum(w)) are decided identically.
school_quantile <- function(x, w, p, sw = sum(w)) {
  if (any(zero <- w == 0)) {
    w <- w[!zero]
    x <- x[!zero]
  }
  cumw <- cumsum(w)
  sel <- cumw <= p * sw
  pos <- if (any(sel)) max(which(sel)) else 1L
  posnext <- if (pos == length(x)) pos else pos + 1L
  wlow <- p - cumw[pos] / sw
  if (wlow <= 0) (x[pos] + x[posnext]) / 2 else x[posnext]
}

# survey:::qrule_math() (lower value), offered for completeness.
math_quantile <- function(x, w, p, sw = sum(w)) {
  if (any(zero <- w == 0)) {
    w <- w[!zero]
    x <- x[!zero]
  }
  cumw <- cumsum(w)
  sel <- cumw <= p * sw
  pos <- if (any(sel)) max(which(sel)) else 1L
  posnext <- if (pos == length(x)) pos else pos + 1L
  wlow <- p - cumw[pos] / sw
  if (wlow <= 0) x[pos] else x[posnext]
}

# Weighted median per group with Woodruff confidence intervals, as
# srvyr::survey_median(qrule = "school", vartype = "ci") -> survey::svyquantile(
# interval.type = "mean", df = degf(group design)).
#   y      : numeric vector (NA = missing)
#   levels : confidence level(s) needed for this variable x group spec
# Returns list(est, n, nw, ci = list(level -> list(low, upp)), failed)
kernel_median <- function(y, gs, des, lonely, levels, qrule = "school") {
  G <- gs$G
  qfun <- if (qrule == "school") school_quantile else math_quantile
  valid <- !is.na(y)
  V <- matrix(as.integer(valid), ncol = 1L)
  n_g <- complete_groups(rowsum(V, gs$gid, reorder = TRUE), G)[, 1]
  nw_g <- complete_groups(rowsum(des$w * V, gs$gid, reorder = TRUE), G)[, 1]
  est <- rep(NA_real_, G)
  ci <- stats::setNames(lapply(levels, function(l) list(low = rep(NA_real_, G), upp = rep(NA_real_, G))),
                        as.character(levels))
  if (!any(valid)) {
    return(list(est = est, n = n_g, nw = nw_g, ci = ci, failed = NULL))
  }
  idx <- which(valid)
  gg <- gs$gid[idx]
  yy <- y[idx]
  ww <- des$w[idx]
  # total weight per group, summed in original row order (survey:::qs / svymean)
  # (base sum(), not data.table's GForce sum: it must round exactly like survey)
  psum <- vapply(split(ww, factor(gg, levels = seq_len(G))), sum, numeric(1))
  psum[!seq_len(G) %in% gg] <- NA_real_

  o <- order(gg, yy, method = "radix")
  dt <- data.table(g = gg[o], y = yy[o], w = ww[o])
  qd <- dt[, list(q = qfun(y, w, 0.5, psum[.BY[[1]]])), by = g]
  est[qd$g] <- qd$q

  # Woodruff: CDF at the median and its SE, i.e. svymean(x <= qhat) on the
  # group design. The proportion and its influence values are computed with
  # survey's arithmetic (colSums(I * w / sum(w)), original row order) because
  # the CI branch "p + t * se > 1" is decided on the last bits when the median
  # is the largest value of the group.
  ind <- as.numeric(yy <= est[gg])
  pd <- data.table(g = gg, w = ww, I = ind)[, {
    s <- psum[.BY[[1]]]
    list(p = sum((I * w) / s))
  }, by = g]
  phat <- rep(NaN, G)
  phat[pd$g] <- pd$p
  Z <- matrix(0, des$n, 1L)
  Z[idx, 1] <- ((ind - phat[gg]) * ww) / psum[gg]
  v <- domain_variance(Z, V, gs, des, lonely)
  se <- sqrt(v[, 1])
  km <- list(failed = attr(v, "failed"))
  dfg <- group_df(V, gs, des)[, 1]

  for (li in seq_along(levels)) {
    alpha <- round(1 - levels[li], 7)          # srvyr::survey_quantile()
    pc <- ci_bounds(phat, se, 1 - alpha, dfg)  # survey:::woodruffCI(): confint(m, level = 1 - alpha, df)
    pl <- pc$low
    pu <- pc$upp
    cd <- dt[, {
      j <- .BY[[1]]
      list(
        lo = if (is.nan(pl[j]) || pl[j] < 0) NaN else qfun(y, w, pl[j], psum[j]),
        up = if (is.nan(pu[j]) || pu[j] > 1) NaN else qfun(y, w, pu[j], psum[j])
      )
    }, by = g]
    ci[[li]]$low[cd$g] <- cd$lo
    ci[[li]]$upp[cd$g] <- cd$up
  }
  list(est = est, n = n_g, nw = nw_g, ci = ci, failed = km$failed)
}
