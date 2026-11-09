# The port reproduces the frozen pipeline (iHEA_2027, freeze-ihea2027-abstract-v3).
# Fixtures: data-raw/pipeline_fixtures.R (pipeline fits on a 400-row synthetic
# train, refitted here under the same seed; plain numbers only).

.fx <- function(id) readRDS(test_path("fixtures", paste0("pipeline_", id, ".rds")))
.refit <- function(m, fx) withr::with_seed(fx$fit_seed, m$fit(fx$train, fx$set, fx$variant, NULL))

test_that("R06 matches the pipeline (fit, mean and draws)", {
  skip_if_not_installed("ordinal"); skip_if_not_installed("pbivnorm"); skip_if_not_installed("Matrix")
  fx <- .fx("R06"); m <- .p2e_method("R06"); vs <- p2e_value_set(fx$vt)
  fit <- .refit(m, fx)
  expect_equal(lapply(fit$dims, function(d) as.numeric(d$alpha)), fx$core$alpha, tolerance = 1e-10)
  expect_equal(lapply(fit$dims, function(d) as.numeric(d$beta)), fx$core$beta, tolerance = 1e-10)
  expect_equal(unclass(fit$coef), fx$core$coef, tolerance = 1e-10)
  expect_equal(unclass(fit$R), fx$core$R, tolerance = 1e-10)
  expect_equal(m$mean_utility(fit, fx$newdata, vs), fx$mean, tolerance = 1e-10)
  n <- nrow(fx$newdata)
  prof <- withr::with_seed(fx$seed, m$simulate(fit, fx$newdata, fx$M, seq_len(n)))
  expect_equal(unname(.score(prof, vs)), fx$draws$u, tolerance = 1e-12)
})

test_that("R09 matches the pipeline (draws exactly, Monte Carlo mean within 0.01)", {
  skip_if_not_installed("xgboost")
  fx <- .fx("R09"); m <- .p2e_method("R09"); vs <- p2e_value_set(fx$vt)
  fit <- .refit(m, fx)
  n <- nrow(fx$newdata)
  prof <- withr::with_seed(fx$seed, m$simulate(fit, fx$newdata, fx$M, seq_len(n)))
  expect_equal(unname(.score(prof, vs)), fx$draws$u, tolerance = 1e-12)
  expect_equal(m$mean_utility(fit, fx$newdata, vs), fx$mean, tolerance = 0.01)
})

test_that("D03 matches the pipeline (beta, out-of-fold pool, point and draws)", {
  skip_if_not_installed("quantreg")
  fx <- .fx("D03"); m <- .p2e_method("D03"); vs <- p2e_value_set(fx$vt)
  fit <- .refit(m, fx)
  expect_equal(fit$beta, fx$core$beta, tolerance = 1e-10)
  expect_equal(fit$oof_mu, fx$core$oof_mu, tolerance = 1e-10)
  expect_equal(fit$oof_resid, fx$core$oof_resid, tolerance = 1e-10)
  n <- nrow(fx$newdata); lo <- .vs_min(vs)
  expect_equal(.clamp(m$point(fit, fx$newdata), vs), fx$mean, tolerance = 1e-12)
  d <- withr::with_seed(fx$seed, m$draw(fit, fx$newdata, fx$M, seq_len(n), "as_published", lo))
  expect_equal(unname(.clamp(d, vs)), fx$draws$u, tolerance = 1e-12)
})

test_that("D12 matches the pipeline (terms, beta, sigma, point and draws)", {
  fx <- .fx("D12"); m <- .p2e_method("D12"); vs <- p2e_value_set(fx$vt)
  fit <- .refit(m, fx)
  expect_identical(fit$terms_kept, fx$core$terms_kept)
  expect_identical(names(fit$beta), fx$core$beta_names)
  expect_equal(as.numeric(fit$beta), fx$core$beta, tolerance = 1e-10)
  expect_equal(as.numeric(fit$sigma), fx$core$sigma, tolerance = 1e-10)
  n <- nrow(fx$newdata); lo <- .vs_min(vs)
  expect_equal(.clamp(m$point(fit, fx$newdata), vs), fx$mean, tolerance = 1e-12)
  d <- withr::with_seed(fx$seed, m$draw(fit, fx$newdata, fx$M, seq_len(n), "as_published", lo))
  expect_equal(unname(.clamp(d, vs)), fx$draws$u, tolerance = 1e-12)
})
