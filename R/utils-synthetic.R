#' Generate a synthetic MSNA-like household dataset
#'
#' Produces a dataset with the features that matter for analysis speed and for
#' the edge cases of the estimators: admin1/admin2/population-group strata,
#' optional clusters, stratum weights, select-one questions with skip-logic
#' missingness, select-multiple questions (parent + 0/1 choice columns),
#' numeric questions (with ties, for medians), ratio numerators/denominators
#' with zeros and missing values, and wide "filler" columns (real MSNA
#' datasets have hundreds or thousands of columns, which is what made the
#' old tool slow).
#'
#' @param n Number of households (rows).
#' @param n_select_one,n_select_multiple,n_numeric Number of generated
#'   questions of each type (in addition to a few named ones).
#' @param n_filler Number of extra text columns not used in the analysis.
#' @param n_admin1,n_admin2_per_admin1 Size of the admin hierarchy.
#' @param cluster_size Average households per cluster (0 = no clusters).
#' @param sm_separator Separator for select-multiple choice columns.
#' @param edge_cases If `TRUE`, inject edge cases: a stratum with a single
#'   PSU, missing group values, a group where a variable is entirely missing,
#'   an all-missing question, single-observation groups.
#' @param seed Random seed.
#'
#' @return A list with `data` (data frame), `loa` (a list of analysis for all
#'   generated questions, `group_var = NA`), and `disaggregations` (suggested
#'   group variables).
#' @export
#' @examples
#' syn <- make_synthetic_msna(500, n_filler = 10)
#' dim(syn$data)
make_synthetic_msna <- function(n = 5000,
                                n_select_one = 20,
                                n_select_multiple = 8,
                                n_numeric = 6,
                                n_filler = 200,
                                n_admin1 = 6,
                                n_admin2_per_admin1 = 8,
                                cluster_size = 12,
                                sm_separator = ".",
                                edge_cases = FALSE,
                                seed = 2026) {
  set.seed(seed)
  admin1 <- sprintf("state_%02d", sample.int(n_admin1, n, replace = TRUE))
  admin2 <- paste0(admin1, "_lga_", sprintf("%02d", sample.int(n_admin2_per_admin1, n, replace = TRUE)))
  pop_group <- sample(c("host", "idp", "returnee"), n, replace = TRUE, prob = c(.55, .35, .10))
  strata <- paste(admin2, pop_group, sep = "_")
  # stratum weights: population / interviews, rescaled to mean 1
  st <- unique(strata)
  pop <- stats::setNames(round(stats::runif(length(st), 2000, 60000)), st)
  cnt <- table(strata)
  w_st <- pop[st] / as.numeric(cnt[st])
  weight <- as.numeric(w_st[strata])
  weight <- weight / mean(weight)
  cluster <- if (cluster_size > 0) {
    paste(strata, sprintf("c%03d", stats::ave(seq_len(n), strata, FUN = function(i) {
      ceiling(seq_along(i) / cluster_size)
    })), sep = "_")
  } else as.character(seq_len(n))

  d <- data.frame(
    uuid = sprintf("uuid_%07d", seq_len(n)),
    admin1 = admin1,
    admin2 = admin2,
    pop_group = pop_group,
    strata = strata,
    cluster_id = cluster,
    weight = weight,
    hoh_gender = sample(c("female", "male"), n, replace = TRUE, prob = c(.35, .65)),
    setting = sample(c("urban", "rural", "camp"), n, replace = TRUE, prob = c(.3, .6, .1)),
    hoh_age_group = sample(c("18-29", "30-59", "60+"), n, replace = TRUE, prob = c(.3, .55, .15)),
    stringsAsFactors = FALSE
  )
  d$hh_size <- pmax(1L, as.integer(stats::rpois(n, 6)))
  d$income <- round(stats::rlnorm(n, 10, 1), -2)
  d$income[stats::runif(n) < .08] <- NA
  d$expenditure_food <- round(stats::rlnorm(n, 9.5, .8), -2)
  d$expenditure_total <- d$expenditure_food + round(stats::rlnorm(n, 9, 1), -2)
  d$expenditure_total[stats::runif(n) < .05] <- NA
  d$num_children <- stats::rpois(n, 2.2)
  d$num_children[stats::runif(n) < .03] <- NA
  size <- ifelse(is.na(d$num_children), 0L, d$num_children)
  d$num_children_school <- ifelse(!is.na(d$num_children) & d$num_children > 0,
                                  stats::rbinom(n, size, .6), NA)

  loa <- list(
    data.frame(analysis_type = "prop_select_one", analysis_var = c("hoh_gender", "setting", "hoh_age_group")),
    data.frame(analysis_type = c("mean", "median"), analysis_var = rep(c("hh_size", "income", "expenditure_food"), each = 2))
  )

  for (i in seq_len(n_select_one)) {
    L <- sample(2:8, 1)
    v <- sample(paste0("opt_", letters[seq_len(L)]), n, replace = TRUE, prob = stats::runif(L, .2, 1))
    v[stats::runif(n) < stats::runif(1, 0, .4)] <- NA
    nm <- sprintf("so_%03d", i)
    d[[nm]] <- v
    loa[[length(loa) + 1L]] <- data.frame(analysis_type = "prop_select_one", analysis_var = nm)
  }
  for (i in seq_len(n_select_multiple)) {
    K <- sample(3:10, 1)
    nm <- sprintf("sm_%03d", i)
    ch <- matrix(as.integer(stats::runif(n * K) < rep(stats::runif(K, .05, .7), each = n)), n, K)
    ch[cbind(seq_len(n), sample.int(K, n, replace = TRUE))] <- 1L
    miss <- stats::runif(n) < stats::runif(1, 0, .3)
    parent <- apply(ch, 1, function(r) paste(paste0("c", which(r == 1)), collapse = " "))
    parent[miss] <- NA
    d[[nm]] <- parent
    for (k in seq_len(K)) {
      x <- ch[, k]
      x[miss] <- NA
      d[[paste0(nm, sm_separator, "c", k)]] <- x
    }
    loa[[length(loa) + 1L]] <- data.frame(analysis_type = "prop_select_multiple", analysis_var = nm)
  }
  for (i in seq_len(n_numeric)) {
    nm <- sprintf("num_%03d", i)
    x <- if (i %% 2) round(stats::rgamma(n, 2, .1)) else round(stats::rnorm(n, 50, 15), 1)
    x[stats::runif(n) < stats::runif(1, 0, .3)] <- NA
    d[[nm]] <- x
    loa[[length(loa) + 1L]] <- data.frame(analysis_type = c("mean", "median"), analysis_var = nm)
  }
  for (i in seq_len(n_filler)) {
    d[[sprintf("filler_%04d", i)]] <- sample(c("x", "y", "z", NA), n, replace = TRUE)
  }

  if (edge_cases) {
    # single-PSU stratum (lonely PSU) and a 2-row stratum
    d$strata[1] <- "lonely_stratum"
    d$cluster_id[1] <- "lonely_cluster"
    d$strata[2:3] <- "tiny_stratum"
    d$cluster_id[2:3] <- c("tiny_c1", "tiny_c2")
    # missing group values
    d$hoh_gender[sample.int(n, max(2, n %/% 100))] <- NA
    # a group where a variable is entirely missing
    d$income[d$setting == "camp"] <- NA
    d$so_001[d$admin1 == unique(d$admin1)[1]] <- NA
    # an all-missing select one / numeric
    d$all_missing_text <- NA_character_
    d$all_missing_num <- NA_real_
    # single-observation groups
    d$rare_group <- ifelse(seq_len(n) %in% 4:5, c("only_4", "only_5")[match(seq_len(n), 4:5)], "common")
    # numeric with heavy ties and a constant group
    d$ties <- sample(c(1, 2, 2, 3, 3, 3), n, replace = TRUE)
    d$ties[d$pop_group == "returnee"] <- 7
    loa[[length(loa) + 1L]] <- data.frame(
      analysis_type = c("prop_select_one", "mean", "median", "mean", "median"),
      analysis_var = c("all_missing_text", "all_missing_num", "all_missing_num", "ties", "ties"))
  }

  loa <- do.call(rbind, loa)
  loa$group_var <- NA_character_
  loa$level <- 0.95
  loa$analysis_var_numerator <- NA_character_
  loa$analysis_var_denominator <- NA_character_
  loa$numerator_NA_to_0 <- NA
  loa$filter_denominator_0 <- NA
  ratios <- data.frame(
    analysis_type = "ratio", analysis_var = NA_character_, group_var = NA_character_, level = 0.95,
    analysis_var_numerator = c("num_children_school", "expenditure_food"),
    analysis_var_denominator = c("num_children", "expenditure_total"),
    numerator_NA_to_0 = c(TRUE, TRUE),
    filter_denominator_0 = c(TRUE, TRUE)
  )
  loa <- rbind(loa, ratios)
  list(
    data = d,
    loa = loa,
    disaggregations = c("admin1", "admin2", "pop_group", "hoh_gender", "setting",
                        "admin1, pop_group", "admin1, hoh_gender")
  )
}
