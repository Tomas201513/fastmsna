# The variance engine against the survey package directly (no analysistools).

engine_mean_se <- function(des_fast, y, domain) {
  gs <- fastmsna:::build_group_spec(des_fast, NA)
  V <- matrix(as.integer(!is.na(y) & domain), ncol = 1)
  Y <- matrix(ifelse(V[, 1] == 1L, y, 0), ncol = 1) # kernels expect 0 outside the domain
  k <- fastmsna:::kernel_mean(Y, V, gs, des_fast, fastmsna:::lonely_options(des_fast))
  c(est = k$est[1, 1], se = sqrt(k$var[1, 1]))
}

survey_mean_se <- function(des, y, domain) {
  des$variables$.y <- y
  sub <- subset(des, domain & !is.na(y))
  m <- survey::svymean(~.y, sub)
  c(est = unname(stats::coef(m)), se = unname(survey::SE(m)))
}

test_that("domain means and SEs match survey for all lonely PSU options and designs", {
  skip_if_not_installed("survey")
  set.seed(42)
  n <- 240
  d <- data.frame(
    s = rep(paste0("s", 1:12), each = 20),
    psu = paste0("p", rep(1:60, each = 4)),
    w = stats::runif(n, 0.5, 3),
    y = stats::rnorm(n, 10, 3),
    dom = stats::runif(n) < 0.6
  )
  d$y[sample.int(n, 30)] <- NA
  # a stratum with a single PSU and a domain that leaves single PSUs
  d$s[1:4] <- "lonely"
  d$fpc <- 50 # population size (PSUs) per stratum

  designs <- list(
    rows = function() survey::svydesign(ids = ~1, strata = ~s, weights = ~w, data = d),
    clustered = function() survey::svydesign(ids = ~psu, strata = ~s, weights = ~w, data = d),
    fpc = function() survey::svydesign(ids = ~psu, strata = ~s, weights = ~w, fpc = ~fpc, data = d)
  )
  for (lp in c("adjust", "remove", "certainty", "average")) {
    for (adj in c(TRUE, FALSE)) {
      old <- options(survey.lonely.psu = lp, survey.adjust.domain.lonely = adj)
      for (nm in names(designs)) {
        des <- suppressWarnings(designs[[nm]]())
        fd <- as_fast_design(des)
        a <- suppressWarnings(survey_mean_se(des, d$y, d$dom))
        b <- engine_mean_se(fd, d$y, d$dom)
        expect_equal(b, a, tolerance = 1e-10, info = paste(lp, adj, nm))
      }
      options(old)
    }
  }
})

test_that("ratio SE matches survey::svyratio on a domain", {
  skip_if_not_installed("survey")
  set.seed(1)
  n <- 200
  d <- data.frame(s = rep(1:8, each = 25), w = stats::runif(n, 1, 4),
                  num = stats::rpois(n, 2), den = stats::rpois(n, 5) + 1,
                  dom = rep(c(TRUE, FALSE), 100))
  des <- survey::svydesign(ids = ~1, strata = ~s, weights = ~w, data = d)
  r <- survey::svyratio(~num, ~den, subset(des, dom))
  fd <- as_fast_design(des)
  gs <- fastmsna:::build_group_spec(fd, NA)
  V <- matrix(as.integer(d$dom), ncol = 1)
  k <- fastmsna:::kernel_ratio(matrix(as.numeric(d$num)), matrix(as.numeric(d$den)), V, gs, fd,
                                fastmsna:::lonely_options(fd))
  expect_equal(k$est[1, 1], unname(stats::coef(r)), tolerance = 1e-12)
  expect_equal(sqrt(k$var[1, 1]), unname(survey::SE(r))[1], tolerance = 1e-10)
})

test_that("medians and Woodruff CIs match survey::svyquantile(qrule = 'school')", {
  skip_if_not_installed("survey")
  set.seed(9)
  n <- 300
  d <- data.frame(s = rep(1:6, each = 50), w = rep(c(1, 1.5, 2, 1, 3, 0.7), each = 50),
                  y = round(stats::rgamma(n, 2, 0.2)), g = sample(c("a", "b", "c"), n, TRUE))
  d$y[sample.int(n, 20)] <- NA
  des <- survey::svydesign(ids = ~1, strata = ~s, weights = ~w, data = d)
  fd <- as_fast_design(des)
  gs <- fastmsna:::build_group_spec(fd, "g")
  k <- fastmsna:::kernel_median(d$y, gs, fd, fastmsna:::lonely_options(fd), levels = 0.9)
  for (gi in seq_len(gs$G)) {
    sub <- subset(des, g == gs$values[gi] & !is.na(y))
    q <- survey::svyquantile(~y, sub, quantiles = 0.5, qrule = "school", ci = TRUE, alpha = 0.1)
    expect_equal(k$est[gi], unname(q$y[1, 1]))
    expect_equal(c(k$ci[["0.9"]]$low[gi], k$ci[["0.9"]]$upp[gi]), unname(q$y[1, 2:3]))
  }
})

test_that("median CIs match svyquantile on heavily tied data (median = maximum knife-edge)", {
  skip_if_not_installed("survey")
  set.seed(1)
  n <- 1200
  d <- data.frame(s = sample(sprintf("s%02d", 1:20), n, TRUE), g = sample(sprintf("g%02d", 1:40), n, TRUE))
  wst <- stats::setNames(round(stats::runif(20, 0.3, 50), sample(0:3, 20, TRUE)), sprintf("s%02d", 1:20))
  d$w <- unname(wst[d$s])
  d$y <- sample(c(0, 1, 2, 3, 3, 7, 7, 7, 7, 7), n, TRUE)
  d$y[stats::runif(n) < .1] <- NA
  with_survey_options("adjust", TRUE, {
    des <- survey::svydesign(ids = ~1, strata = ~s, weights = ~w, data = d)
    fd <- as_fast_design(des)
    gs <- fastmsna:::build_group_spec(fd, "g")
    k <- fastmsna:::kernel_median(d$y, gs, fd, fastmsna:::lonely_options(fd), levels = 0.9)
    for (gi in seq_len(gs$G)) {
      sub <- des[d$g == gs$values[gi] & !is.na(d$y), ]
      q <- suppressWarnings(survey::svyquantile(~y, sub, 0.5, qrule = "school", ci = TRUE, alpha = 0.1))$y[1, 1:3]
      expect_equal(c(k$est[gi], k$ci[["0.9"]]$low[gi], k$ci[["0.9"]]$upp[gi]), unname(q),
                   info = gs$values[gi])
    }
  })
})

test_that("a data frame design gives the same results as the srvyr design", {
  skip_if_not_installed("srvyr")
  syn <- make_synthetic_msna(300, n_select_one = 2, n_select_multiple = 1, n_numeric = 1, n_filler = 0)
  loa <- syn$loa
  loa$group_var <- "admin1"
  with_survey_options("adjust", TRUE, {
    a <- create_analysis_fast(srvyr::as_survey_design(syn$data, weights = weight, strata = strata), loa = loa)
    b <- create_analysis_fast(syn$data, loa = loa, weights = "weight", strata = "strata")
  })
  expect_true(compare_analysis_results(a, b, tolerance = 1e-12)$equivalent)
})
