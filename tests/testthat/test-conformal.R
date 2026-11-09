# Ported from iHEA_2027 tests/testthat/test-metrics.R at freeze-ihea2027-abstract-v3 (a79e102); changes: only the cvplus_interval tests are ported, renamed to .cvplus_interval; the conformal tests are new

test_that("cvplus_interval returns -Inf / +Inf when n is too small", {
  ci <- .cvplus_interval(matrix(.8, 2, 2), c(.1, .2), 1:2, alpha = .1)
  expect_equal(ci$lower, c(-Inf, -Inf)); expect_equal(ci$upper, c(Inf, Inf))
})

test_that("cvplus_interval uses each OOF row's own-fold model, per test person", {
  # two folds; fold 1 model predicts .2 / .6, fold 2 predicts .4 / .8 for two people
  tp <- cbind(c(.2, .6), c(.4, .8))
  ci <- .cvplus_interval(tp, c(.1, .1, .3, .3), c(1, 1, 2, 2), alpha = .5)  # n = 4: 2nd and 3rd
  # person 1 lower: .1 .1 .1 .1 (fold 1 -> .2-.1) ; fold 2 -> .4-.3 = .1 : all .1 -> .1
  # person 1 upper: .3 .3 .7 .7 -> 3rd = .7
  expect_equal(ci$lower, c(.1, .5)); expect_equal(ci$upper, c(.7, 1.1))
})

test_that("cvplus_interval one-row and single-vector inputs, signed residuals, bad folds", {
  ci <- .cvplus_interval(matrix(.8, 1, 4), c(.1, .2, .3, .4), 1:4, alpha = .5)
  expect_equal(c(ci$lower, ci$upper), c(.5, 1.1))
  a <- .cvplus_interval(matrix(.8, 1, 4), c(.1, .2, .3, .4), 1:4, alpha = .5)
  b <- .cvplus_interval(c(.8, .8, .8, .8), -c(.1, .2, .3, .4), 1:4, alpha = .5)  # vector = one person; signed
  expect_identical(a, b)
  expect_error(.cvplus_interval(matrix(.8, 1, 4), c(.1, .2), 1:4, .5), "one fold per residual")
  expect_error(.cvplus_interval(matrix(.8, 1, 2), c(.1, .2, .3), c(1, 2, 3), .5), "index the columns")
})

test_that("conformal intervals are ordered, capped and cover about level on held-out data", {
  b <- stand_in_bundle("D12_klapproth", "S3K", "UK"); te <- make_synthetic_data(300, 7, countries = "UK")
  vs <- p2e_value_set("UK"); ci <- .conformal(b, te, vs, 0.9)
  expect_true(all(ci$lower <= ci$upper & ci$lower >= .vs_min(vs) & ci$upper <= 1))
  y <- .score(as.matrix(te[DIMS]), vs)
  expect_gt(mean(y >= ci$lower & y <= ci$upper), 0.8)
})

test_that("response conformal works under a value set never trained on", {
  b <- stand_in_bundle("R06_mvop_dim", "S3"); te <- make_synthetic_data(20, 7)
  usa <- p2e_value_set("USA")
  ci <- .conformal(b, te, usa, 0.9)
  expect_length(ci$lower, 20)
  # hand computation: residuals and fold predictions under USA, then CV+ and the cap
  mu <- .p2e_method("R06_mvop_dim")$mean_utility
  resid <- abs(.score(b$y_profiles, usa) - rowMeans(.score(b$oof$profiles, usa)))
  tp <- vapply(b$folds, function(f) mu(f, te, usa), numeric(nrow(te)))
  ref <- .cvplus_interval(tp, resid, b$fold_id, alpha = 0.1)
  expect_equal(ci$lower, pmax(ref$lower, .vs_min(usa)))
  expect_equal(ci$upper, pmin(ref$upper, 1))
})

test_that("crosswalked D03 intervals are ordered, capped and centred on the crosswalked point", {
  b <- stand_in_bundle("D03_clad", "S3", "UK"); te <- make_synthetic_data(60, 5)
  usa <- p2e_value_set("USA")
  k <- .fit_crosswalk(p2e_value_set("UK"), usa, b$state_freq)
  ci <- .conformal(b, te, usa, 0.9, crosswalk = k)
  expect_true(all(ci$lower <= ci$upper))
  expect_true(all(ci$lower >= .vs_min(usa) & ci$upper <= 1))
  pt <- .apply_crosswalk(.clamp(.p2e_method("D03_clad")$point(b$point, te), b$vs), k, usa)
  expect_gt(mean(pt >= ci$lower & pt <= ci$upper), 0.5)
})

test_that("direct conformal needs a crosswalk for another value set, and a bundle value set", {
  b <- stand_in_bundle("D03_clad", "S3", "UK"); te <- make_synthetic_data(5, 2)
  expect_error(.conformal(b, te, p2e_value_set("USA"), 0.9), "crosswalk is needed")
  b$vs <- NULL
  expect_error(.conformal(b, te, p2e_value_set("UK"), 0.9), "no value set")
})

test_that("level must be one number strictly between 0 and 1", {
  b <- stand_in_bundle("D12_klapproth", "S3K", "UK"); te <- make_synthetic_data(5, 2, countries = "UK")
  vs <- p2e_value_set("UK")
  for (l in list(0, 1, 1.5, -0.1, c(.9, .9), NA_real_, "0.9")) {
    expect_error(.conformal(b, te, vs, l), "level")
  }
})

test_that("blocked and unblocked conformal intervals are identical", {
  b <- stand_in_bundle("D12_klapproth", "S3K", "UK"); te <- make_synthetic_data(10, 3, countries = "UK")
  vs <- p2e_value_set("UK")
  full <- .conformal(b, te, vs, 0.9)
  testthat::local_mocked_bindings(CONFORMAL_BLOCK = 3L)
  expect_identical(.conformal(b, te, vs, 0.9), full)
})

test_that("a bundle with no fold fits cannot give conformal intervals", {
  b <- stand_in_bundle("D12_klapproth", "S3K", "UK"); te <- make_synthetic_data(5, 2, countries = "UK")
  b$folds <- list()
  expect_error(.conformal(b, te, p2e_value_set("UK"), 0.9), "no fold fits.*interval = TRUE")
})
