# Ported in part from iHEA_2027 tests/testthat/test-design.R at freeze-ihea2027-abstract-v3 (a79e102); changes: only tests of build_design, .merge_levels and gaussian lasso_select; dropped ordinal/binomial lasso, state_weights, target_utility, folds, cohorts and model-data tests; added gaussian-only and synthetic-data tests
test_that("design ref reproduces columns on new data", {
  tr <- make_synthetic_data(400, 1); te <- make_synthetic_data(60, 2)
  for (f in c("numeric","dummies")) { a <- build_design(tr, "S1", f); b <- build_design(te, "S1", f, a$ref)
    expect_identical(colnames(a$X), colnames(b$X)) }
  for (f in c("linear","quadratic","spline3","spline3_int")) { a <- build_design(tr, "S3", f); b <- build_design(te, "S3", f, a$ref)
    expect_identical(colnames(a$X), colnames(b$X)) }
})
test_that("sparse item levels are merged", {
  tr <- make_synthetic_data(400, 1)
  tr$P1_1 <- c(rep(1L, 5), rep(2L, 95), rep(3L, 100), rep(4L, 100), rep(5L, 100))  # level 1 has 5 < 20 respondents
  a <- build_design(tr, "S2", "dummies")
  expect_equal(sum(grepl("^P1_1_", colnames(a$X))), 3)          # levels {1,2} merged as reference; dummies for 3, 4, 5
  te <- tr[1:10, ]; te$P1_1 <- 1L
  expect_true(all(build_design(te, "S2", "dummies", a$ref)$X[, grep("^P1_1_", colnames(a$X))] == 0))  # level 1 maps to the merged reference
})
test_that("merge rule: toward the middle, repeated until every level has 20", {
  m <- function(counts) .merge_levels(counts, 20)
  expect_identical(m(c(5, 95, 100, 100, 100)), c(1L, 1L, 3L, 4L, 5L))
  expect_identical(m(c(100, 100, 100, 100, 5)), c(1L, 2L, 3L, 4L, 4L))
  expect_identical(m(c(2, 3, 100, 100, 100)), c(1L, 1L, 1L, 4L, 5L))    # 1 -> 2, then {1,2} (5) -> 3
  expect_identical(m(c(100, 100, 5, 100, 100)), c(1L, 2L, 2L, 4L, 5L))  # sparse middle: nearest neighbour
  expect_identical(m(c(0, 100, 100, 100, 100)), c(1L, 1L, 3L, 4L, 5L))  # absent level merged
  expect_identical(m(rep(100, 5)), 1:5)
})
test_that("design is learned from ref data only (new levels, extreme T-scores)", {
  A <- make_synthetic_data(400, 1)
  A$P1_1 <- sample(c(2L, 3L, 4L, 5L), 400, TRUE)   # A lacks level 1
  B <- make_synthetic_data(80, 2)
  B$P1_1 <- 1L                                      # B has a level A lacks
  B$T_PF[1] <- 500; B$T_PI[2] <- -300; B$AGE[3] <- 400L
  for (f in c("dummies")) {
    a <- build_design(A, "S2", f); b <- build_design(B, "S2", f, a$ref)
    expect_identical(colnames(a$X), colnames(b$X)); expect_true(all(is.finite(b$X)))
  }
  for (f in c("linear","quadratic","spline3","spline3_int")) {
    a <- build_design(A, "S3", f); b <- build_design(B, "S3", f, a$ref)
    expect_identical(colnames(a$X), colnames(b$X)); expect_identical(ncol(a$X), ncol(b$X))
    expect_true(all(is.finite(b$X)))
    # applying ref to A reproduces A's own matrix exactly
    expect_equal(build_design(A, "S3", f, a$ref)$X, a$X)
  }
  # ref decides the columns: a sparse level in B must not change them
  B2 <- B; B2$P1_1[1:70] <- 5L
  expect_identical(colnames(build_design(B2, "S2", "dummies", build_design(A, "S2", "dummies")$ref)$X),
                   colnames(build_design(A, "S2", "dummies")$X))
  expect_error(build_design(A, "S2", "spline3_int"))
})
test_that("spline knots survive a T-score ceiling with heavy ties", {
  A <- make_synthetic_data(300, 8)
  A$T_PF <- ifelse(seq_len(300) <= 200, 57, 40 + (seq_len(300) %% 15))  # two thirds at the maximum
  a <- build_design(A, "S3", "spline3")
  expect_true(all(is.finite(a$X)))
  b <- build_design(make_synthetic_data(50, 9), "S3", "spline3", a$ref)
  expect_identical(colnames(a$X), colnames(b$X))
})

test_that("lasso_select picks the signal variable", {
  skip_if_not_installed("glmnet")
  set.seed(1); n <- 200
  X <- matrix(rnorm(n * 6), n, 6, dimnames = list(NULL, paste0("v", 1:6)))
  fold <- rep(1:5, length.out = n)
  yg <- 3 * X[, 1] + rnorm(n)
  expect_true("v1" %in% lasso_select(X, yg, "gaussian", fold))
})
test_that("lasso_select accepts character fold ids", {
  skip_if_not_installed("glmnet")
  set.seed(3); n <- 150
  X <- matrix(rnorm(n * 5), n, 5, dimnames = list(NULL, paste0("v", 1:5)))
  y <- 3 * X[, 1] + rnorm(n)
  fi <- rep(1:5, length.out = n); fc <- c("FR", "DE", "UK", "XX", "YY")[fi]
  # factor() orders character ids alphabetically, so compare with that integer coding
  fi2 <- as.integer(factor(fc))
  set.seed(9); a <- lasso_select(X, y, "gaussian", fc)
  set.seed(9); b <- lasso_select(X, y, "gaussian", fi2)
  expect_identical(a, b)
})
test_that("lasso_select supports the gaussian family only", {
  X <- matrix(rnorm(40), 20, 2, dimnames = list(NULL, c("a", "b")))
  expect_error(lasso_select(X, rnorm(20), "binomial", rep(1:5, 4)), "gaussian")
  expect_error(lasso_select(X, rnorm(20), "ordinal", rep(1:5, 4)), "gaussian")
})
test_that("synthetic data has the pipeline columns and scores U_own on the value sets", {
  d <- make_synthetic_data(60, 1, c("FR", "DE", "UK"))
  expect_identical(names(d), .SYN_COLS)
  for (cc in c("FR", "DE", "UK")) {
    i <- d$country == cc
    expect_equal(d$U_own[i], .score(as.matrix(d[i, DIMS]), p2e_value_set(cc)))
  }
  expect_error(make_synthetic_data(10, 1, "XX"), "subset")
})
