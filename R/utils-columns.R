# Column preparation (cached on the design) ---------------------------------
#
# The old workflow converted columns inside every srvyr pipeline (and copied
# the whole dataset for every group). Here every column is converted at most
# once per design, and only the columns an indicator needs are touched.

cache_get <- function(des, key, fun) {
  val <- des$cache[[key]]
  if (is.null(val)) {
    val <- fun()
    assign(key, val, envir = des$cache)
  }
  val
}

get_column <- function(des, var) {
  x <- des$variables[[var]]
  if (is.null(x)) stop("The variable '", var, "' is not present in the dataset.", call. = FALSE)
  x
}

# Numeric version of a column for mean / median / ratio.
# Text columns that hold numbers (e.g. data read with col_types = "text") are
# converted; anything else non-numeric is an error, as in srvyr.
get_numeric <- function(des, var) {
  cache_get(des, paste0("num\r", var), function() {
    x <- get_column(des, var)
    if (is.factor(x)) {
      stop("Factor columns are not supported for numeric analyses ('", var, "'). ",
           "Convert it with as.numeric(as.character(x)).", call. = FALSE)
    }
    if (is.numeric(x) || is.logical(x)) return(as.numeric(x))
    if (is.character(x)) {
      y <- suppressWarnings(as.numeric(x))
      bad <- is.na(y) & !is.na(x) & nzchar(trimws(x))
      if (any(bad)) {
        stop("The variable '", var, "' is not numeric (e.g. value '", x[which(bad)[1]],
             "'). Convert it before running a mean / median / ratio.", call. = FALSE)
      }
      y[is.na(x) | !nzchar(trimws(x))] <- NA_real_
      return(y)
    }
    stop("The variable '", var, "' has an unsupported type (", class(x)[1], ").", call. = FALSE)
  })
}

# Select-multiple dummy columns: analysistools uses as.numeric(); in addition,
# text "TRUE"/"FALSE" (logical columns written to Excel and read back as text)
# are mapped to 1/0 instead of silently becoming NA.
get_sm_child <- function(des, var) {
  cache_get(des, paste0("smc\r", var), function() {
    x <- get_column(des, var)
    if (is.character(x)) {
      y <- suppressWarnings(as.numeric(x))
      bad <- which(is.na(y) & !is.na(x))
      if (length(bad)) {
        up <- toupper(trimws(x[bad]))
        tf <- up %in% c("TRUE", "FALSE")
        y[bad[tf]] <- as.numeric(up[tf] == "TRUE")
      }
      return(y)
    }
    if (is.factor(x)) return(suppressWarnings(as.numeric(as.character(x))))
    as.numeric(x)
  })
}

# Select-one codes: integer codes in dplyr::group_by() order (C locale,
# factors by level), NA kept, plus the level labels.
get_codes <- function(des, var) {
  cache_get(des, paste0("cod\r", var), function() {
    x <- get_column(des, var)
    codes <- as.integer(frankv(x, ties.method = "dense", na.last = "keep"))
    L <- if (all(is.na(codes))) 0L else max(codes, na.rm = TRUE)
    first <- which(!duplicated(codes) & !is.na(codes))
    first <- first[order(codes[first])]
    list(codes = codes, labels = as.character(x[first]), L = L)
  })
}

# select-multiple children of a parent (analysistools: dplyr::starts_with(),
# which is case-insensitive), in dataset column order.
sm_children <- function(des, parent, sm_separator) {
  prefix <- paste0(parent, sm_separator)
  nms <- names(des$variables)
  nms[startsWith(tolower(nms), tolower(prefix))]
}

# Number of complete observations for the variables of a LOA row (for the
# skipped-analysis log, as create_analysis_safe()).
n_complete <- function(des, vars) {
  vars <- unique(vars[!is.na(vars) & nzchar(vars)])
  vars <- intersect(vars, names(des$variables))
  if (!length(vars)) return(NA_integer_)
  ok <- rep(TRUE, des$n)
  for (v in vars) ok <- ok & !is.na(des$variables[[v]])
  sum(ok)
}

# Split units into chunks so that n x (columns in chunk) stays below
# `max_cells` (bounds peak memory of the n x k matrices).
chunk_units <- function(ncols, n, max_cells) {
  per_chunk <- max(1, floor(max_cells / max(1, n)))
  chunk <- integer(length(ncols))
  cur <- 1L
  used <- 0
  for (i in seq_along(ncols)) {
    if (used > 0 && used + ncols[i] > per_chunk) {
      cur <- cur + 1L
      used <- 0
    }
    chunk[i] <- cur
    used <- used + ncols[i]
  }
  chunk
}
