# Grouping layer -----------------------------------------------------------
#
# A "group spec" is computed once per distinct `group_var` string and re-used
# by every indicator that is disaggregated by it. It holds:
#   * gid    : group (domain) id of every row, numbered in dplyr::group_by()
#              order (each variable ascending, C locale, NA last);
#   * labels : group_var / group_var_value strings and the group part of the
#              analysis key;
#   * cells  : domain x stratum (and domain x PSU) indexes used by the
#              variance engine.

# analysistools::char_to_vector() equivalent.
split_group_var <- function(string) {
  if (length(string) > 1) stop("The group_var to be turned into a vector is already a vector.")
  if (is.na(string)) return(character(0))
  out <- trimws(strsplit(string, ",", fixed = TRUE)[[1]])
  if (!length(out)) out <- ""
  if (any(out == "") && length(out) > 1) {
    stop("The group_var seems to have empty value, please check the inputs values")
  }
  if (length(out) == 1 && out == "") return(character(0))
  out
}

# analysistools::create_group_var() equivalent (stringr::str_squish semantics).
group_var_label <- function(group_var) {
  if (is.na(group_var)) return(NA_character_)
  squish(gsub(",", " %/% ", group_var, fixed = TRUE))
}

squish <- function(x) gsub("\\s+", " ", gsub("^\\s+|\\s+$", "", x, perl = TRUE), perl = TRUE)

# stringr::str_split(fixed) semantics: keeps trailing empty pieces.
str_split_keep <- function(x, sep) {
  lapply(strsplit(paste0(x, sep, "\001"), sep, fixed = TRUE), function(v) v[-length(v)])
}

# Normalise LOA group_var entries: NA / "" / whitespace -> NA.
normalise_group_var <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & !nzchar(trimws(x))] <- NA_character_
  x
}

build_group_spec <- function(des, group_var) {
  vars <- split_group_var(group_var)
  n <- des$n
  label <- if (length(vars)) group_var_label(group_var) else NA_character_

  if (!length(vars)) {
    gid <- rep.int(1L, n)
    G <- 1L
    values <- NA_character_
    key_group <- "NA %/% NA"
  } else {
    missing_vars <- setdiff(vars, names(des$variables))
    if (length(missing_vars)) {
      stop("The following group variables are not present in the dataset: ",
           paste(missing_vars, collapse = ", "), call. = FALSE)
    }
    cols <- lapply(vars, function(v) des$variables[[v]])
    gid <- if (length(cols) == 1L) {
      frankv(cols[[1]], ties.method = "dense", na.last = TRUE)
    } else {
      frankv(cols, ties.method = "dense", na.last = TRUE)
    }
    gid <- as.integer(gid)
    G <- if (n) max(gid) else 0L
    first <- which(!duplicated(gid))
    first <- first[order(gid[first])]
    # tidyr::unite(sep = " %/% ") semantics: paste(), NA -> "NA"
    values <- do.call(paste, c(lapply(cols, function(cc) cc[first]), sep = " %/% "))
    # analysistools::adding_analysis_key(): split label and value on " %/% ",
    # paste pairwise, collapse with " -/- ".
    key_names <- str_split_keep(label, " %/% ")[[1]]
    key_group <- vapply(str_split_keep(values, " %/% "), function(v) {
      paste(paste(key_names, v, sep = " %/% "), collapse = " -/- ")
    }, character(1))
  }

  # domain x stratum cells (sorted by group, then stratum)
  cell <- if (G > 1L) frankv(list(gid, des$strata), ties.method = "dense") else des$strata
  cell <- as.integer(cell)
  C <- if (n) max(cell) else 0L
  fc <- which(!duplicated(cell))
  cell_g <- integer(C)
  cell_h <- integer(C)
  cell_g[cell[fc]] <- gid[fc]
  cell_h[cell[fc]] <- des$strata[fc]

  spec <- list(
    group_var = group_var,
    vars = vars,
    label = label,
    gid = gid,
    G = G,
    values = values,
    key_group = key_group,
    cell = cell,
    C = C,
    cell_g = cell_g,
    cell_h = cell_h,
    cell_nPSU = des$nPSU[cell_h],
    cell_f = des$fpc_f[cell_h],
    cell_scale = des$scale[cell_h]
  )

  if (des$has_clusters) {
    pc <- as.integer(frankv(list(gid, des$cluster), ties.method = "dense"))
    fp <- which(!duplicated(pc))
    pc_cell <- integer(max(pc))
    pc_cell[pc[fp]] <- cell[fp]
    pc_g <- integer(max(pc))
    pc_g[pc[fp]] <- gid[fp]
    spec$pc <- pc
    spec$pc_cell <- pc_cell
    spec$pc_g <- pc_g
  }
  spec
}

# Fetch (and cache on the run context) the group spec for a group_var string.
get_group_spec <- function(ctx, group_var) {
  key <- if (is.na(group_var)) "\001NA" else group_var
  spec <- ctx$gs[[key]]
  if (is.null(spec)) {
    spec <- build_group_spec(ctx$des, group_var)
    ctx$gs[[key]] <- spec
  }
  spec
}
