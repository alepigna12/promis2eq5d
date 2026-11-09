# Tests for the isotonic crosswalk between value sets.

test_that("crosswalk to itself is the identity on the table", {
  uk <- p2e_value_set("UK"); k <- .fit_crosswalk(uk, uk)
  expect_equal(.apply_crosswalk(uk$u, k, uk), uk$u)
})
test_that("a value set with ties at the minimum is anchored", {
  t <- cbind(.all_states_df(), utility = p2e_value_set("DE")$u)
  t$utility[t$utility < -0.6] <- min(t$utility)            # several states share the minimum
  v <- p2e_value_set(t); k <- .fit_crosswalk(p2e_value_set("UK"), v)
  expect_equal(k$y[1], min(v$u)); expect_equal(tail(k$y, 1), 1)
  expect_true(all(diff(k$y) >= 0))
})
test_that("source choice picks the lowest weighted RMSE", {
  us <- p2e_value_set("USA"); w <- rep(1, 3125)
  rmse <- sapply(c("FR", "DE", "UK"), function(s) { f <- p2e_value_set(s)
    sqrt(mean((.apply_crosswalk(f$u, .fit_crosswalk(f, us, w), us) - us$u)^2)) })
  expect_equal(.choose_crosswalk_source(us, weights = w), names(which.min(rmse)))
})
test_that("a source value set chooses itself, and NULL weights mean equal weights", {
  for (s in c("FR", "DE", "UK")) expect_identical(.choose_crosswalk_source(p2e_value_set(s)), s)
  us <- p2e_value_set("USA")
  expect_identical(.choose_crosswalk_source(us), .choose_crosswalk_source(us, weights = rep(2, 3125)))
})
test_that("crosswalk knots are isotonic, cover both ranges and pool tied values", {
  codes <- c("FR", "DE", "UK")
  for (from in codes) for (to in setdiff(codes, from)) {
    f <- p2e_value_set(from); t <- p2e_value_set(to)
    k <- .fit_crosswalk(f, t)
    expect_identical(names(k), c("x", "y"))
    expect_identical(attr(k, "from"), from); expect_identical(attr(k, "to"), to)
    expect_identical(k$x, sort(unique(f$u)))                 # one knot per distinct from value
    expect_true(all(diff(k$y) >= 0))
    expect_identical(c(k$x[1], tail(k$x, 1)), c(.vs_min(f), 1))
    expect_equal(c(k$y[1], tail(k$y, 1)), c(.vs_min(t), 1))
    expect_equal(.apply_crosswalk(.vs_min(f), k, t), .vs_min(t))
  }
  # Ties are pooled before the monotone fit, so knots do not depend on state order
  f <- p2e_value_set("DE"); t <- p2e_value_set("FR"); o <- rev(seq_len(3125))
  expect_equal(.fit_knots(f$u[o], t$u[o], rep(1, 3125), range(f$u), range(t$u)),
               { k <- .fit_crosswalk(f, t); data.frame(x = k$x, y = k$y) }, tolerance = 1e-12)
})
test_that("isotonic fit is the weighted least-squares monotone fit", {
  fit <- function(x, y, w = rep(1, length(x)), xr = range(x), yr = range(y)) .fit_knots(x, y, w, xr, yr)
  expect_equal(fit(c(1, 2, 2, 3), c(3, 1, 2, 0)), data.frame(x = c(1, 2, 3), y = rep(1.5, 3)))
  # x = 3 dropped by zero weight, then anchored at (max x, max y)
  expect_equal(fit(c(1, 2, 2, 3), c(3, 1, 2, 0), c(1, 1, 1, 0)),
               data.frame(x = c(1, 2, 3), y = c(2, 2, 3)))
  expect_equal(fit(1:4, c(1, 3, 2, 4))$y, c(1, 2.5, 2.5, 4))
  expect_equal(fit(1:4, c(1, 3, 2, 4), c(1, 1, 3, 1))$y, c(1, 2.25, 2.25, 4))   # weighted pool
})
test_that("weights are validated", {
  f <- p2e_value_set("DE"); t <- p2e_value_set("FR")
  expect_error(.fit_crosswalk(f, t, rep(1, 10)), "weights")
  expect_error(.fit_crosswalk(f, t, c(Inf, rep(1, 3124))), "weights")
  expect_error(.fit_crosswalk(f, t, c(-1, rep(1, 3124))), "weights")
  expect_error(.fit_crosswalk(f, t, rep(0, 3125)), "zero")
  expect_error(.choose_crosswalk_source(t, weights = rep(1, 10)), "weights")
})
test_that("dropping an end state with zero weight still anchors both ends", {
  f <- p2e_value_set("DE"); t <- p2e_value_set("FR")
  w <- rep(1, 3125); w[3125] <- 0; w[1] <- 0                # 55555 and 11111
  k <- .fit_crosswalk(f, t, w)
  expect_equal(c(k$x[1], tail(k$x, 1)), c(.vs_min(f), 1))
  expect_equal(c(k$y[1], tail(k$y, 1)), c(.vs_min(t), 1))
  expect_true(all(diff(k$y) >= 0))
})
test_that(".apply_crosswalk keeps shape, clamps and validates", {
  f <- p2e_value_set("DE"); t <- p2e_value_set("FR"); k <- .fit_crosswalk(f, t)
  u <- matrix(c(-3, .2, .9, 3, 0, 1), 2)
  out <- .apply_crosswalk(u, k, t)
  expect_identical(dim(out), dim(u))
  expect_true(all(out >= .vs_min(t) & out <= 1))
  expect_identical(out[, 2], .apply_crosswalk(u[, 2], k, t))
  expect_true(all(diff(.apply_crosswalk(seq(.vs_min(f), 1, length.out = 200), k, t)) >= 0))
  expect_error(.apply_crosswalk(c(.5, NA), k, t), "missing")
  expect_error(.apply_crosswalk(.5, data.frame(x = c(1, 0), y = c(.5, 1)), t), "knots")
  expect_error(.apply_crosswalk(.5, data.frame(x = c(0, 1), y = c(1, .5)), t), "knots")
})
