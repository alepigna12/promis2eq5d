test_that("gap rounding empties the gap and preserves the mean", {
  withr::local_seed(1); vs <- p2e_value_set("UK"); b <- 1 - .vs_gap(vs)
  u <- matrix(runif(2e5, b, 1), 1000); V <- matrix(runif(2e5), 1000)
  r <- .round_gap(u, vs, V)
  expect_true(all(r == 1 | abs(r - b) < 1e-12))
  expect_equal(mean(r), mean(u), tolerance = 1e-3)
  outside <- matrix(c(0.5, 1, -0.2, b), 2); expect_identical(.round_gap(outside, vs, matrix(0.3, 2, 2)), outside)
})

test_that("censored-normal location makes the draws' mean equal the target", {
  lo <- -0.567; s <- 0.12
  for (t in c(0.2, 0.85, 0.97, 0.999)) {
    m <- .censored_normal_location(t, s, lo, 1)
    z <- withr::with_seed(1, rnorm(4e5)); expect_lt(abs(mean(pmin(pmax(m + s * z, lo), 1)) - t), 1e-3)
  }
  expect_equal(.censored_normal_location(1, s, lo, 1), Inf)
})
test_that("censored-normal location returns -Inf at the minimum and handles vector sigma", {
  lo <- -0.567; hi <- 1
  expect_equal(.censored_normal_location(lo, 0.12, lo, hi), -Inf)
  expect_equal(.censored_normal_location(c(lo - 0.3, 1.2, 0.5), 0.12, lo, hi)[c(1, 2)], c(-Inf, Inf))
  v <- .censored_normal_location(c(0.3, 0.8), c(0.1, 0.2), lo, hi)
  expect_equal(v, c(.censored_normal_location(0.3, 0.1, lo, hi), .censored_normal_location(0.8, 0.2, lo, hi)))
  expect_true(all(is.finite(.censored_normal_location(c(0.3, 0.8), 0.15, lo, hi))))
})
