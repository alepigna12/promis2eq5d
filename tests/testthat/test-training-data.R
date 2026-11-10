test_that("p2e_training has the documented shape and no T-scores", {
  d <- .training_data()
  expect_identical(names(d), c(ITEMS, "P9", "AGE", "female", "country", DIMS))
  expect_equal(nrow(d), 4512L); expect_false(anyNA(d))
  expect_true(all(unlist(d[ITEMS]) %in% 1:5)); expect_true(all(d$P9 %in% 0:10))
  expect_true(all(unlist(d[DIMS]) %in% 1:5)); expect_setequal(unique(d$country), c("FR", "DE", "UK"))
  expect_true(all(d$female %in% 0:1)); expect_false(any(grepl("^T_", names(d))))
})
test_that("p2e_training stores plain columns of the digested types", {
  d <- .training_data()
  expect_identical(attributes(d)$names, names(d))
  for (v in names(d)) expect_null(attributes(d[[v]]), info = v)
  expect_true(all(vapply(d[c(ITEMS, "P9", "female", DIMS)], is.integer, logical(1))))
  expect_type(d$AGE, "double"); expect_type(d$country, "character")
})
test_that("util_band follows the pipeline rule", {
  d <- make_synthetic_data(400, 1, countries = c("FR", "DE", "UK"))
  P <- as.matrix(d[DIMS])
  ub <- .util_band(P, d$country)
  fh <- rowSums(P == 1) == 5
  expect_true(all(ub[fh] == "FH")); expect_setequal(unique(ub[!fh]), paste0("Q", 1:4))
  vs <- c(FR = "France", DE = "Germany", UK = "UK")
  for (cc in names(vs)) {
    i <- which(d$country == cc & !fh)
    u <- .score(P[i, , drop = FALSE], p2e_value_set(vs[[cc]]))
    q <- quantile(u, probs = 0:4 / 4, type = 7)
    expect_equal(ub[i], paste0("Q", cut(u, q, include.lowest = TRUE, labels = FALSE)))
  }
})
test_that("the folds of p2e_training are valid for K = 10", {
  d <- .training_data(); d$util_band <- .util_band(as.matrix(d[DIMS]), d$country)
  f <- make_cv_folds(d, v = 10, seed = 20261110L)$fold
  expect_setequal(f, 1:10); expect_lte(diff(range(tabulate(f, 10))), 1)
})
