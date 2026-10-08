# Analysis runners ----------------------------------------------------------
#
# One runner per analysis type. Each runner
#   1. prepares the columns of all its indicators once,
#   2. packs them into n x k matrices (chunked to bound memory),
#   3. for every disaggregation (group spec) requested, runs ONE kernel call
#      for all indicators of the chunk and all groups at once,
#   4. expands the G x k results into the analysistools long format.
#
# All runners return a data.table with columns
#   .row .ord analysis_var analysis_var_value group_var group_var_value
#   .key_group stat stat_low stat_upp n n_total n_w n_w_total

result_columns <- c(".row", ".ord", "analysis_var", "analysis_var_value", "group_var",
                    "group_var_value", ".key_group", "stat", "stat_low", "stat_upp",
                    "n", "n_total", "n_w", "n_w_total")

empty_result <- function() {
  data.table(.row = integer(), .ord = numeric(), analysis_var = character(),
             analysis_var_value = character(), group_var = character(),
             group_var_value = character(), .key_group = character(),
             stat = numeric(), stat_low = numeric(), stat_upp = numeric(),
             n = numeric(), n_total = numeric(), n_w = numeric(), n_w_total = numeric())
}

# Log (or raise) an error for a set of LOA rows.
record_error <- function(ctx, rows, msg) {
  if (!nrow(rows)) return(invisible(NULL))
  if (ctx$on_error == "stop") {
    r <- rows[1]
    stop(sprintf("analysis failed (loa row %d: %s / %s, group_var = %s): %s",
                 r$.row, r$analysis_type, r$.unit_label, r$.gs, msg), call. = FALSE)
  }
  vars_cols <- intersect(c("analysis_var", "analysis_var_numerator", "analysis_var_denominator"), names(rows))
  n_obs <- vapply(seq_len(nrow(rows)), function(i) {
    n_complete(ctx$des, unlist(rows[i, vars_cols, with = FALSE], use.names = FALSE))
  }, integer(1))
  ctx$skipped[[length(ctx$skipped) + 1L]] <- data.table(
    row_index = rows$.row,
    analysis_type = rows$analysis_type,
    analysis_var = ifelse(is.na(rows$.unit_label), "NA", rows$.unit_label),
    group_var = ifelse(is.na(rows$.gs), "NA", rows$.gs),
    n_obs = n_obs,
    error_message = msg
  )
  if (isTRUE(ctx$verbose)) {
    message(sprintf("  ! Skipped %d row(s) (%s / %s): %s", nrow(rows), rows$analysis_type[1],
                    rows$.unit_label[1], msg))
  }
  invisible(NULL)
}

# Attach LOA rows (one or several per unit x group spec) to the long block
# `base` (one row per group x value) and compute confidence intervals.
attach_rows <- function(base, map, gs, ci = TRUE) {
  out <- map[base, on = ".unit", allow.cartesian = TRUE, nomatch = NULL]
  if (ci) {
    b <- ci_bounds(out$stat, out$se, out$.level, out$df)
    out[, `:=`(stat_low = b$low, stat_upp = b$upp)]
  }
  out[, `:=`(group_var = gs$label,
             group_var_value = gs$values[.g],
             .key_group = gs$key_group[.g])]
  out
}

# analysistools::correct_nan_total_is_0()
nan_if_total_0 <- function(dt) {
  dt[n_total == 0, `:=`(stat = NaN, stat_low = NaN, stat_upp = NaN,
                        n_w_total = NaN, n_total = NaN, n_w = NaN)]
  dt
}

unit_map <- function(rows, units_in, gs_key) {
  sel <- rows$.unit %in% units_in & (if (is.na(gs_key)) is.na(rows$.gs) else (!is.na(rows$.gs) & rows$.gs == gs_key))
  rows[sel, list(.unit, .row, .level)]
}

gs_keys_for <- function(rows, units_in) unique(rows$.gs[rows$.unit %in% units_in])

failed_units <- function(failed, units_cols) {
  if (is.null(failed)) return(list())
  cols <- as.integer(names(failed))
  out <- list()
  for (i in seq_along(cols)) {
    u <- units_cols[cols[i]]
    if (is.null(out[[u]])) out[[u]] <- paste0("Stratum (", failed[[i]], ") has only one PSU at stage 1")
  }
  out
}

# ---- mean -------------------------------------------------------------------
run_mean <- function(rows, ctx) {
  des <- ctx$des
  rows[, `:=`(.unit = analysis_var, .unit_label = analysis_var)]
  units <- unique(rows$.unit)
  prep <- list()
  for (u in units) {
    res <- tryCatch(get_numeric(des, u), error = function(e) e)
    if (inherits(res, "error")) record_error(ctx, rows[.unit == u], conditionMessage(res)) else prep[[u]] <- res
  }
  units <- names(prep)
  if (!length(units)) return(empty_result())
  chunks <- chunk_units(rep(1, length(units)), des$n, ctx$max_cells)
  out <- list()
  for (ch in unique(chunks)) {
    cu <- units[chunks == ch]
    k <- length(cu)
    Y <- matrix(0, des$n, k)
    V <- matrix(0L, des$n, k)
    for (j in seq_len(k)) {
      y <- prep[[cu[j]]]
      Y[, j] <- zero_na(y)
      V[, j] <- as.integer(!is.na(y))
    }
    dfs <- design_df(V, des)
    for (gk in gs_keys_for(rows, cu)) {
      gs <- get_group_spec(ctx, gk)
      map <- unit_map(rows, cu, gk)
      cols <- which(cu %in% map$.unit)
      km <- kernel_mean(Y[, cols, drop = FALSE], V[, cols, drop = FALSE], gs, des, ctx$lonely)
      fu <- failed_units(km$failed, cu[cols])
      for (u in names(fu)) {
        record_error(ctx, rows[.unit == u & .row %in% map[.unit == u, .row]], fu[[u]])
        map <- map[.unit != u]
      }
      G <- gs$G
      kk <- length(cols)
      base <- data.table(
        .unit = rep(cu[cols], each = G),
        .g = rep(seq_len(G), kk),
        analysis_var_value = NA_character_,
        stat = as.vector(km$est),
        se = sqrt(as.vector(km$var)),
        df = rep(dfs[cols], each = G),
        n = as.vector(km$n),
        n_w = as.vector(km$nw)
      )
      base[, `:=`(n_total = n, n_w_total = n_w, .ord = .g)]
      out[[length(out) + 1L]] <- attach_rows(base, map, gs)
    }
  }
  res <- rbindlist(out, use.names = TRUE, fill = TRUE)
  res[, analysis_var := .unit]
  nan_if_total_0(res)
}

# ---- median -----------------------------------------------------------------
run_median <- function(rows, ctx) {
  des <- ctx$des
  rows[, `:=`(.unit = analysis_var, .unit_label = analysis_var)]
  units <- unique(rows$.unit)
  out <- list()
  for (u in units) {
    y <- tryCatch(get_numeric(des, u), error = function(e) e)
    if (inherits(y, "error")) {
      record_error(ctx, rows[.unit == u], conditionMessage(y))
      next
    }
    for (gk in gs_keys_for(rows, u)) {
      gs <- get_group_spec(ctx, gk)
      map <- unit_map(rows, u, gk)
      levels <- unique(map$.level)
      km <- kernel_median(y, gs, des, ctx$lonely, levels, qrule = ctx$qrule)
      if (!is.null(km$failed)) {
        record_error(ctx, rows[.row %in% map$.row], paste0("Stratum (", km$failed[[1]], ") has only one PSU at stage 1"))
        next
      }
      G <- gs$G
      base <- data.table(.unit = u, .g = seq_len(G), analysis_var_value = NA_character_,
                         stat = km$est, n = km$n, n_w = km$nw)
      base[, `:=`(n_total = n, n_w_total = n_w, .ord = .g)]
      blk <- attach_rows(base, map, gs, ci = FALSE)
      lv <- as.character(blk$.level)
      gi <- blk$.g
      blk[, `:=`(stat_low = NA_real_, stat_upp = NA_real_)]
      for (l in unique(lv)) {
        idx <- which(lv == l)
        blk[idx, `:=`(stat_low = km$ci[[l]]$low[gi[idx]], stat_upp = km$ci[[l]]$upp[gi[idx]])]
      }
      out[[length(out) + 1L]] <- blk
    }
  }
  if (!length(out)) return(empty_result())
  res <- rbindlist(out, use.names = TRUE, fill = TRUE)
  res[, analysis_var := .unit]
  nan_if_total_0(res)
}

# ---- select one -------------------------------------------------------------
run_select_one <- function(rows, ctx) {
  des <- ctx$des
  rows[, `:=`(.unit = analysis_var, .unit_label = analysis_var)]
  # analysistools cannot group a select-one by itself (the grouping column is
  # renamed to analysis_var_value): keep the same behaviour (error / skip).
  self_grouped <- vapply(seq_len(nrow(rows)), function(i) {
    !is.na(rows$.gs[i]) && rows$analysis_var[i] %in% split_group_var(rows$.gs[i])
  }, logical(1))
  if (any(self_grouped)) {
    record_error(ctx, rows[self_grouped],
                 "The analysis variable is also a grouping variable (not supported for prop_select_one).")
    rows <- rows[!self_grouped]
    if (!nrow(rows)) return(empty_result())
  }
  units <- unique(rows$.unit)
  prep <- list()
  for (u in units) {
    res <- tryCatch(get_codes(des, u), error = function(e) e)
    if (inherits(res, "error")) record_error(ctx, rows[.unit == u], conditionMessage(res)) else prep[[u]] <- res
  }
  units <- names(prep)
  if (!length(units)) return(empty_result())
  ncols <- vapply(units, function(u) max(1L, prep[[u]]$L), numeric(1))
  chunks <- chunk_units(ncols, des$n, ctx$max_cells)
  out <- list()
  for (ch in unique(chunks)) {
    cu <- units[chunks == ch]
    Ls <- vapply(cu, function(u) prep[[u]]$L, integer(1))
    k <- sum(Ls)
    col_unit <- rep(cu, Ls)
    col_level <- unlist(lapply(Ls, seq_len), use.names = FALSE)
    Y <- matrix(0, des$n, k)
    V <- matrix(0L, des$n, k)
    Vu <- matrix(0L, des$n, length(cu))
    NAm <- matrix(0L, des$n, length(cu))
    off <- 0L
    for (j in seq_along(cu)) {
      cd <- prep[[cu[j]]]$codes
      ok <- !is.na(cd)
      Vu[, j] <- as.integer(ok)
      NAm[, j] <- as.integer(!ok)
      if (Ls[j] > 0) {
        cols <- off + seq_len(Ls[j])
        Y[cbind(which(ok), off + cd[ok])] <- 1
        V[, cols] <- as.integer(ok)
        off <- off + Ls[j]
      }
    }
    dfs <- design_df(Vu, des)
    names(dfs) <- cu
    for (gk in gs_keys_for(rows, cu)) {
      gs <- get_group_spec(ctx, gk)
      map <- unit_map(rows, cu, gk)
      um <- cu[cu %in% map$.unit]
      cols <- which(col_unit %in% um)
      G <- gs$G
      na_cnt <- complete_groups(rowsum(NAm[, match(um, cu), drop = FALSE], gs$gid, reorder = TRUE), G)
      colnames(na_cnt) <- um
      blocks <- list()
      if (length(cols)) {
        km <- kernel_mean(Y[, cols, drop = FALSE], V[, cols, drop = FALSE], gs, des, ctx$lonely)
        fu <- failed_units(km$failed, col_unit[cols])
        for (u in names(fu)) {
          record_error(ctx, rows[.row %in% map[.unit == u, .row]], fu[[u]])
          map <- map[.unit != u]
        }
        cu_c <- col_unit[cols]
        nwt <- vapply(um, function(u) {
          cc <- which(cu_c == u)
          if (length(cc)) rowSums(km$sum_wy[, cc, drop = FALSE]) else rep(0, G)
        }, numeric(G))
        nwt <- matrix(nwt, nrow = G)
        colnames(nwt) <- um
        kk <- length(cols)
        lvl <- col_level[cols]
        base <- data.table(
          .unit = rep(cu_c, each = G),
          .g = rep(seq_len(G), kk),
          .l = rep(lvl, each = G),
          stat = as.vector(km$est),
          se = sqrt(as.vector(km$var)),
          n = as.vector(km$sum_y),
          n_total = as.vector(km$n),
          n_w = as.vector(km$sum_wy)
        )
        base[, df := dfs[.unit]]
        base[, n_w_total := nwt[cbind(.g, match(.unit, um))]]
        base <- base[n > 0]
        base[, analysis_var_value := prep_labels(prep, .unit, .l)]
        base[, .na_row := FALSE]
        blocks[[1]] <- base
      }
      # rows for missing values of the analysis variable (n = number of NA)
      na_base <- data.table(.unit = rep(um, each = G), .g = rep(seq_len(G), length(um)),
                            n = as.vector(na_cnt))
      na_base <- na_base[n > 0]
      if (nrow(na_base)) {
        na_base[, `:=`(.l = 0L, analysis_var_value = NA_character_, stat = NaN, se = NaN, df = NaN,
                       n_total = NaN, n_w = NaN, n_w_total = NaN, .na_row = TRUE)]
        blocks[[length(blocks) + 1L]] <- na_base
      }
      if (!length(blocks)) next
      base <- rbindlist(blocks, use.names = TRUE, fill = TRUE)
      setorderv(base, c(".unit", ".g", ".na_row", ".l"))
      base[, .ord := seq_len(.N), by = .unit]
      blk <- attach_rows(base, map, gs)
      blk[.na_row == TRUE, `:=`(stat_low = NaN, stat_upp = NaN)]
      out[[length(out) + 1L]] <- blk
    }
  }
  if (!length(out)) return(empty_result())
  res <- rbindlist(out, use.names = TRUE, fill = TRUE)
  res[, analysis_var := .unit]
  res
}

prep_labels <- function(prep, unit, level) {
  out <- character(length(unit))
  for (u in unique(unit)) {
    i <- which(unit == u)
    out[i] <- prep[[u]]$labels[level[i]]
  }
  out
}

# ---- select multiple --------------------------------------------------------
run_select_multiple <- function(rows, ctx) {
  des <- ctx$des
  sep <- ctx$sm_separator
  rows[, `:=`(.unit = analysis_var, .unit_label = analysis_var)]
  units <- unique(rows$.unit)
  prep <- list()
  for (u in units) {
    res <- tryCatch({
      kids <- sm_children(des, u, sep)
      if (!length(kids)) {
        stop("No select multiple choice columns found for '", u, "' with separator '", sep, "'.", call. = FALSE)
      }
      parent_ok <- !is.na(get_column(des, u))
      vals <- lapply(kids, function(k) get_sm_child(des, k))
      labels <- gsub(paste0(u, sep), "", kids)
      list(kids = kids, labels = labels, parent_ok = parent_ok, vals = vals)
    }, error = function(e) e)
    if (inherits(res, "error")) record_error(ctx, rows[.unit == u], conditionMessage(res)) else prep[[u]] <- res
  }
  units <- names(prep)
  if (!length(units)) return(empty_result())
  ncols <- vapply(units, function(u) length(prep[[u]]$kids), numeric(1))
  chunks <- chunk_units(ncols, des$n, ctx$max_cells)
  out <- list()
  for (ch in unique(chunks)) {
    cu <- units[chunks == ch]
    Ks <- vapply(cu, function(u) length(prep[[u]]$kids), integer(1))
    k <- sum(Ks)
    col_unit <- rep(cu, Ks)
    col_child <- unlist(lapply(Ks, seq_len), use.names = FALSE)
    Y <- matrix(0, des$n, k)
    V <- matrix(0L, des$n, k)
    Pm <- matrix(0L, des$n, length(cu))
    off <- 0L
    for (j in seq_along(cu)) {
      p <- prep[[cu[j]]]
      Pm[, j] <- as.integer(p$parent_ok)
      for (c in seq_along(p$kids)) {
        y <- p$vals[[c]]
        ok <- p$parent_ok & !is.na(y)
        y[!ok] <- 0
        Y[, off + c] <- y
        V[, off + c] <- as.integer(ok)
      }
      off <- off + Ks[j]
    }
    dfs <- design_df(Pm, des)
    names(dfs) <- cu
    for (gk in gs_keys_for(rows, cu)) {
      gs <- get_group_spec(ctx, gk)
      map <- unit_map(rows, cu, gk)
      um <- cu[cu %in% map$.unit]
      cols <- which(col_unit %in% um)
      G <- gs$G
      km <- kernel_mean(Y[, cols, drop = FALSE], V[, cols, drop = FALSE], gs, des, ctx$lonely)
      fu <- failed_units(km$failed, col_unit[cols])
      for (u in names(fu)) {
        record_error(ctx, rows[.row %in% map[.unit == u, .row]], fu[[u]])
        map <- map[.unit != u]
      }
      kk <- length(cols)
      cu_c <- col_unit[cols]
      base <- data.table(
        .unit = rep(cu_c, each = G),
        .g = rep(seq_len(G), kk),
        .c = rep(col_child[cols], each = G),
        stat = as.vector(km$est),
        se = sqrt(as.vector(km$var)),
        n = as.vector(km$sum_y),
        n_total = as.vector(km$n),
        n_w = as.vector(km$sum_wy),
        n_w_total = as.vector(km$nw)
      )
      base[, df := dfs[.unit]]
      base[, analysis_var_value := sm_labels(prep, .unit, .c)]
      base <- nan_if_total_0(base)
      # "NA" rows: number of rows with a missing parent in the group
      na_cnt <- complete_groups(rowsum(1L - Pm[, match(um, cu), drop = FALSE], gs$gid, reorder = TRUE), G)
      na_base <- data.table(.unit = rep(um, each = G), .g = rep(seq_len(G), length(um)),
                            n = as.vector(na_cnt))
      na_base <- na_base[n > 0]
      if (nrow(na_base)) {
        na_base[, `:=`(.c = .Machine$integer.max, analysis_var_value = "NA",
                       stat = NA_real_, se = NA_real_, df = NA_real_,
                       n_total = NA_real_, n_w = NA_real_, n_w_total = NA_real_)]
        base <- rbindlist(list(base, na_base), use.names = TRUE, fill = TRUE)
      }
      setorderv(base, c(".unit", ".g", ".c"))
      base[, .ord := seq_len(.N), by = .unit]
      out[[length(out) + 1L]] <- attach_rows(base, map, gs)
    }
  }
  if (!length(out)) return(empty_result())
  res <- rbindlist(out, use.names = TRUE, fill = TRUE)
  res[, analysis_var := .unit]
  res
}

sm_labels <- function(prep, unit, child) {
  out <- character(length(unit))
  for (u in unique(unit)) {
    i <- which(unit == u)
    out[i] <- prep[[u]]$labels[child[i]]
  }
  out
}

# ---- ratio ------------------------------------------------------------------
run_ratio <- function(rows, ctx) {
  des <- ctx$des
  rows[, .unit := paste(analysis_var_numerator, analysis_var_denominator,
                        numerator_NA_to_0, filter_denominator_0, sep = "\r")]
  rows[, .unit_label := paste(analysis_var_numerator, "%/%", analysis_var_denominator)]
  units <- unique(rows$.unit)
  prep <- list()
  for (u in units) {
    r1 <- rows[.unit == u][1]
    res <- tryCatch({
      num <- get_numeric(des, r1$analysis_var_numerator)
      den <- get_numeric(des, r1$analysis_var_denominator)
      na0 <- as.logical(r1$numerator_NA_to_0)
      f0 <- as.logical(r1$filter_denominator_0)
      if (isTRUE(na0)) num[is.na(num)] <- 0
      keep <- !is.na(den)
      if (isTRUE(f0)) keep <- keep & den != 0
      if (!isTRUE(na0)) {
        keep <- if (isTRUE(ctx$ratio_legacy)) !is.na(num) else keep & !is.na(num)
      }
      valid <- keep & !is.na(num) & !is.na(den)
      list(num = num, den = den, keep = keep, valid = valid,
           num_name = r1$analysis_var_numerator, den_name = r1$analysis_var_denominator)
    }, error = function(e) e)
    if (inherits(res, "error")) record_error(ctx, rows[.unit == u], conditionMessage(res)) else prep[[u]] <- res
  }
  units <- names(prep)
  if (!length(units)) return(empty_result())
  chunks <- chunk_units(rep(2, length(units)), des$n, ctx$max_cells)
  out <- list()
  for (ch in unique(chunks)) {
    cu <- units[chunks == ch]
    k <- length(cu)
    Yn <- matrix(0, des$n, k)
    Yd <- matrix(0, des$n, k)
    V <- matrix(0L, des$n, k)
    Rm <- matrix(0L, des$n, k)
    for (j in seq_len(k)) {
      p <- prep[[cu[j]]]
      Yn[, j] <- p$num
      Yd[, j] <- p$den
      V[, j] <- as.integer(p$valid)
      Rm[, j] <- as.integer(p$keep)
    }
    dfs <- design_df(Rm, des)
    for (gk in gs_keys_for(rows, cu)) {
      gs <- get_group_spec(ctx, gk)
      map <- unit_map(rows, cu, gk)
      cols <- which(cu %in% map$.unit)
      kr <- kernel_ratio(Yn[, cols, drop = FALSE], Yd[, cols, drop = FALSE], V[, cols, drop = FALSE],
                         gs, des, ctx$lonely)
      fu <- failed_units(kr$failed, cu[cols])
      for (u in names(fu)) {
        record_error(ctx, rows[.row %in% map[.unit == u, .row]], fu[[u]])
        map <- map[.unit != u]
      }
      G <- gs$G
      kk <- length(cols)
      Rk <- Rm[, cols, drop = FALSE]
      nn <- complete_groups(rowsum(Rk, gs$gid, reorder = TRUE), G)
      nw <- complete_groups(rowsum(des$w * Rk, gs$gid, reorder = TRUE), G)
      base <- data.table(
        .unit = rep(cu[cols], each = G),
        .g = rep(seq_len(G), kk),
        analysis_var_value = "NA %/% NA",
        stat = as.vector(kr$est),
        se = sqrt(as.vector(kr$var)),
        df = rep(dfs[cols], each = G),
        n = as.vector(nn),
        n_w = as.vector(nw)
      )
      base[, `:=`(n_total = n, n_w_total = n_w, .ord = .g)]
      out[[length(out) + 1L]] <- attach_rows(base, map, gs)
    }
  }
  if (!length(out)) return(empty_result())
  res <- rbindlist(out, use.names = TRUE, fill = TRUE)
  lab <- vapply(prep, function(p) paste(p$num_name, "%/%", p$den_name), character(1))
  res[, analysis_var := lab[.unit]]
  nan_if_total_0(res)
}
