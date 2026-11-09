test_that("gap rounding empties the gap and preserves the mean", {
  withr::local_seed(1); vs <- p2e_value_set("UK"); b <- 1 - .vs_gap(vs)
  u <- matrix(runif(2e5, b, 1), 1000); V <- matrix(runif(2e5), 1000)
  r <- .round_gap(u, vs, V)
  expect_true(all(r == 1 | abs(r - b) < 1e-12))
  expect_equal(mean(r), mean(u), tolerance = 1e-3)
  outside <- matrix(c(0.5, 1, -0.2, b), 2); expect_identical(.round_gap(outside, vs, matrix(0.3, 2, 2)), outside)
})
